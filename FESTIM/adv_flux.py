import festim as F
import ufl
from scifem import assemble_scalar
from dolfinx import fem


class AdvectionFlux(F.SurfaceFlux):
    def __init__(self, field, surface, velocity_field, filename=None):
        super().__init__(field=field, surface=surface, filename=filename)

        self.velocity_field = velocity_field

    @property
    def title(self):
        return f"{self.field.name} advective flux surface {self.surface.id}"

    def compute(self, u, ds: ufl.Measure, entity_maps=None):
        """Computes the value of the flux at the surface

        Args:
            ds (ufl.Measure): surface measure of the model
        """

        # obtain mesh normal from field
        # if case multispecies, solution is an index, use sub_function_space
        if isinstance(u, ufl.indexed.Indexed):
            mesh = self.field.sub_function_space.mesh
        else:
            mesh = u.function_space.mesh
        n = ufl.FacetNormal(mesh)

        surface_flux = assemble_scalar(
            fem.form(
                -self.D * ufl.dot(ufl.grad(u), n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )
        advective_flux = assemble_scalar(
            fem.form(
                u * ufl.inner(self.velocity_field, n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        self.value = surface_flux + advective_flux
        self.data.append(self.value)
