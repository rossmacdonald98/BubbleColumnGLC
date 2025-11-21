import gmsh
import dolfinx
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
from dolfinx.io import VTXWriter, XDMFFile
import festim as F
from dolfinx import fem
import ufl
from basix.ufl import element
import numpy as np
from dolfinx import cpp as _cpp
from dolfinx.log import set_log_level, LogLevel
import h_transport_materials as htm

from adv_flux import AdvectionFlux


def evaluate_stabilisation_term(mesh, u, delta):
    """See more at https://www.comsol.com/blogs/understanding-stabilization-methods"""

    # evaluate Cell size
    tdim = mesh.topology.dim
    num_cells = mesh.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)
    mesh_ = _cpp.mesh.Mesh_float64(
        mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
    )
    h = _cpp.mesh.h(mesh_, tdim, cells)
    V0 = fem.functionspace(mesh, ("DG", 0))
    h_as_function = fem.Function(V0)
    h_as_function.x.array[:] = h

    # Compute magnitude of velocity
    v_mag = ufl.sqrt(ufl.dot(u, u))

    D_art = delta * v_mag * h_as_function

    return D_art


def generate_mesh(r_inner=0.1, r_tube=0.11, length=0.4):
    assert r_tube > r_inner, "Tube radius must be larger than fluid radius"
    gmsh.initialize()
    gmsh.model.add("mwe")

    factory = gmsh.model.occ

    fluid = factory.addRectangle(0, 0, 0, length, r_inner)
    tube = factory.addRectangle(0, 0, 0, length, r_tube)
    walls, interface = factory.cut(
        [(2, tube)], [(2, fluid)], removeObject=True, removeTool=False
    )

    factory.synchronize()
    interface_curves = [s[1] for s in interface[0]]
    interface_tag = gmsh.model.addPhysicalGroup(1, interface_curves)
    gmsh.model.setPhysicalName(1, interface_tag, "interface")

    fluid_tag = gmsh.model.addPhysicalGroup(2, [fluid])
    gmsh.model.setPhysicalName(2, fluid_tag, "fluid")

    wall_tag = gmsh.model.addPhysicalGroup(2, [walls[0][1]])
    gmsh.model.setPhysicalName(2, wall_tag, "wall")

    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0002)

    surfaces = gmsh.model.getEntities(dim=1)

    inlet_surfaces = []
    outlet_surfaces = []

    for s in surfaces:
        com = gmsh.model.occ.getCenterOfMass(s[0], s[1])
        # inlet at x=0, outlet at x=length
        if abs(com[0]) < 1e-6:
            if abs(com[1]) < r_inner + 1e-6:
                inlet_surfaces.append(s[1])
        elif abs(com[0]) > length - 1e-6:
            if abs(com[1]) < r_inner + 1e-6:
                outlet_surfaces.append(s[1])

    # inlet_tag = gmsh.model.addPhysicalGroup(1, inlet_surfaces)
    # gmsh.model.setPhysicalName(1, inlet_tag, "inlet")

    # outlet_tag = gmsh.model.addPhysicalGroup(1, outlet_surfaces)
    # gmsh.model.setPhysicalName(1, outlet_tag, "outlet")

    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write("mwe.msh")
    mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=2)
    my_mesh = mesh_data.mesh
    gmsh.finalize()

    return my_mesh, mesh_data, fluid_tag, wall_tag


def calculate_mesh_peclet(mesh, velocity, D):
    # evaluate Cell size
    tdim = mesh.topology.dim
    num_cells = mesh.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)
    mesh_ = _cpp.mesh.Mesh_float64(
        mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
    )
    h = _cpp.mesh.h(mesh_, tdim, cells)

    V0 = fem.functionspace(mesh, ("DG", 0))
    V1 = fem.functionspace(mesh, ("CG", 1))
    h_as_function = fem.Function(V0)
    h_as_function.x.array[:] = h

    # Compute magnitude of velocity
    v_mag = ufl.sqrt(ufl.dot(velocity, velocity))

    # evaluate Peclet number
    Pe_expr = v_mag * h_as_function / D

    Pe_local = fem.Function(V1)
    Pe_local.interpolate(fem.Expression(Pe_expr, V1.element.interpolation_points))

    return Pe_local


r_inner = 0.01  # fluid radius (m)
r_tube = r_inner + 0.001  # tube wall radius (m)
length = 0.1  # cylinder length (m)

scaling_factor = 1e-15

temperature = 800  # K

my_mesh, mesh_data, fluid_tag, wall_tag = generate_mesh(r_inner, r_tube, length)

my_model = F.HydrogenTransportProblemDiscontinuous()
my_model.mesh = F.Mesh(my_mesh)

el = element("Lagrange", my_mesh.topology.cell_name(), 2, shape=(my_mesh.geometry.dim,))


V = dolfinx.fem.functionspace(my_model.mesh.mesh, el)

velocity = dolfinx.fem.Function(V)

entities = mesh_data.cell_tags.find(2)
# flow_rate = 0.0001  # m3/s
average_velocity = 1e-3  # flow_rate / (np.pi * r_inner**2)
velocity.interpolate(
    lambda x: (
        average_velocity * (1.0 - (x[1] / r_inner) ** 2),
        np.full_like(x[0], 0.0),
    ),
    cells0=entities,
)
fluid_mat = htm.LIPB
tube_mat = htm.VANADIUM
D_fluid_diff = htm.diffusivities.filter(material=fluid_mat).mean()
K_S_fluid = htm.solubilities.filter(material=fluid_mat).mean()
D_solid_diff = htm.diffusivities.filter(material=tube_mat).mean()
K_S_solid = htm.solubilities.filter(material=tube_mat).mean()

# add stabilization term for diffusion
D_art = evaluate_stabilisation_term(mesh=my_mesh, u=velocity, delta=0.25)

D_expr = D_fluid_diff.value(temperature).magnitude + D_art
V = fem.functionspace(my_mesh, ("CG", 1))
D_fluid = fem.Function(V)
D_fluid.interpolate(fem.Expression(D_expr, V.element.interpolation_points))


my_writer = VTXWriter(MPI.COMM_WORLD, "velocity_field.bp", velocity, "BP5")
my_writer.write(t=0)

my_writer_2 = VTXWriter(MPI.COMM_WORLD, "D_field.bp", D_fluid, "BP5")
my_writer_2.write(t=0)


dummy_fluid = F.Material(
    D=D_fluid,
    K_S_0=K_S_fluid.pre_exp.magnitude * scaling_factor,
    E_K_S=K_S_fluid.act_energy.magnitude,
)
dummy_tube = F.Material(
    D_0=D_solid_diff.pre_exp.magnitude,
    E_D=D_solid_diff.act_energy.magnitude,
    K_S_0=K_S_solid.pre_exp.magnitude * scaling_factor,
    E_K_S=K_S_solid.act_energy.magnitude,
)

Pe_local = calculate_mesh_peclet(mesh=my_mesh, velocity=velocity, D=D_fluid)
writer = dolfinx.io.VTXWriter(MPI.COMM_WORLD, "Pe_local.bp", Pe_local, "BP5")
writer.write(t=0.0)


inlet = F.SurfaceSubdomain(
    id=4, locator=lambda x: np.logical_and(np.isclose(x[0], 0.0), x[1] < r_inner)
)
outlet = F.SurfaceSubdomain(
    id=5,
    locator=lambda x: np.logical_and(np.isclose(x[0], length), x[1] < r_inner + 1e-6),
)
top = F.SurfaceSubdomain(id=6, locator=lambda x: np.isclose(x[1], r_tube))
fluid = F.VolumeSubdomain(
    id=fluid_tag, material=dummy_fluid, locator=lambda x: x[1] < r_inner + 1e-6
)
tube = F.VolumeSubdomain(
    id=wall_tag, material=dummy_tube, locator=lambda x: x[1] > r_inner - 1e-6
)

my_model.subdomains = [inlet, outlet, top, fluid, tube]

my_model.surface_to_volume = {inlet: fluid, outlet: fluid, top: tube}

H = F.Species("H", subdomains=my_model.volume_subdomains)
my_model.species = [H]
my_model.interfaces = [F.Interface(id=1, subdomains=[fluid, tube], penalty_term=100)]
my_model.method_interface = F.InterfaceMethod.nitsche
my_model.temperature = temperature

advection_term = F.AdvectionTerm(
    velocity=velocity,
    subdomain=fluid,
    species=H,
)

my_model.advection_terms = [advection_term]

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1e15 * scaling_factor, species=H),
    F.FixedConcentrationBC(subdomain=top, value=0, species=H),
]
my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

concentration_field_fluid = F.VTXSpeciesExport(
    filename=f"H_fluid.bp", field=H, subdomain=fluid
)
concentration_field_tube = F.VTXSpeciesExport(
    filename=f"H_tube.bp", field=H, subdomain=tube
)

permeated_flux = F.SurfaceFlux(field=H, surface=top)
outlet_flux = AdvectionFlux(field=H, surface=outlet, velocity_field=velocity)
outlet_flux_diff = F.SurfaceFlux(field=H, surface=outlet)
inlet_flux = AdvectionFlux(field=H, surface=inlet, velocity_field=velocity)
inlet_flux_diff = F.SurfaceFlux(field=H, surface=inlet)
my_model.exports = [
    concentration_field_fluid,
    concentration_field_tube,
    permeated_flux,
    outlet_flux,
    outlet_flux_diff,
    inlet_flux,
    inlet_flux_diff,
]


my_model.initialise()


with XDMFFile(MPI.COMM_WORLD, "facet_tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(my_mesh)
    xdmf.write_meshtags(my_model.facet_meshtags, my_mesh.geometry)

with XDMFFile(MPI.COMM_WORLD, "volume_tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(my_mesh)
    xdmf.write_meshtags(my_model.volume_meshtags, my_mesh.geometry)
set_log_level(LogLevel.INFO)
my_model.run()


print(f"Permeated flux: {permeated_flux.value:.2e} atoms/m2/s")
print(f"Outlet flux: {outlet_flux.value:.2e} atoms/m2/s")
print(f"Outlet diffusive flux: {outlet_flux_diff.value:.2e} atoms/m2/s")
print(f"Inlet flux: {inlet_flux.value:.2e} atoms/m2/s")
print(f"Inlet diffusive flux: {inlet_flux_diff.value:.2e} atoms/m2/s")

print(f"Total out: {outlet_flux.value + permeated_flux.value:.2e} atoms/m2/s")
