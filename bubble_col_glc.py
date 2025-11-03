import numpy as np
from scipy.integrate import solve_bvp
import scipy.constants as const
from scipy.optimize import root_scalar


def solve(params):

    def solve_tritium_extraction(dimensionless_params, y_T2_in, elements):
        """
        Solves the BVP for tritium extraction in a bubble column.

        Args:
            params (dict): A dictionary containing the dimensionless parameters:
                        Bo_l, phi_l, Bo_g, phi_g, psi, nu.
            y_T2_in (float): Inlet tritium molar fraction in the gas phase, y_T2(0-).

        Returns:
            sol: The solution object from scipy.integrate.solve_bvp.
        """

        Bo_l = dimensionless_params["Bo_l"]
        phi_l = dimensionless_params["phi_l"]
        Bo_g = dimensionless_params["Bo_g"]
        phi_g = dimensionless_params["phi_g"]
        psi = dimensionless_params["psi"]
        nu = dimensionless_params["nu"]

        def ode_system(xi, S):
            """
            Defines the system of 4 first-order ODEs.
            S[0] = x_T  (dimensionless liquid concentration)
            S[1] = dx_T/d(xi)
            S[2] = y_T2 (dimensionless gas concentration)
            S[3] = dy_T2/d(xi)
            """
            x_T, dx_T_dxi, y_T2, dy_T2_dxi = S

            # Dimensionless driving force theta. Eq. 8.8
            theta = x_T - np.sqrt(np.maximum(0, (1 - psi * xi) * y_T2 / nu))

            # Equation for d(S[0])/d(xi) = d(x_T)/d(xi)
            dS0_dxi = dx_T_dxi

            # Equation for d(S[1])/d(xi) = d^2(x_T)/d(xi)^2
            dS1_dxi = Bo_l * (phi_l * theta - dx_T_dxi)  # Eq. 9.1.4

            # Equation for d(S[2])/d(xi) = d(y_T2)/d(xi)
            dS2_dxi = dy_T2_dxi

            # Equation for d(S[3])/d(xi) = d^2(y_T2)/d(xi)^2 from eq (11)
            # Avoid division by zero if (1 - psi * xi) is close to zero at xi=1
            denominator = 1 - psi * xi
            denominator = np.where(np.isclose(denominator, 0), 1e-9, denominator)

            term1 = (1 + 2 * psi / Bo_g) * dy_T2_dxi
            term2 = phi_g * theta
            dS3_dxi = (Bo_g / denominator) * (term1 - term2)

            return np.vstack((dS0_dxi, dS1_dxi, dS2_dxi, dS3_dxi))

        def boundary_conditions(Sa, Sb):
            """
            Defines the boundary conditions for the BVP.
            Sa: solution at xi = 0 (liquid outlet)
            Sb: solution at xi = 1 (liquid inlet)
            """
            # Residuals that should be zero for a valid solution.
            # Based on equations (16) and (17) in the paper.

            # At xi = 0: dx_T/d(xi) = 0
            res1 = Sa[1]  # Eq. 10.1

            # At xi = 1: x_T(1) = 1 - (1/Bo_l) * dx_T/d(xi)|_1
            res2 = Sb[0] - (1 - (1 / Bo_l) * Sb[1]) # Eq. 10.2

            # At xi = 0: y_T2(0) = y_T2(0-) + (1/Bo_g) * dy_T2/d(xi)|_0
            res3 = Sa[2] - y_T2_in - (1 / Bo_g) * Sa[3]  # Eq. 10.3

            # At xi = 1: dy_T2/d(xi) = 0
            res4 = Sb[3]  # Eq. 10.4

            return np.array([res1, res2, res3, res4])

        # Set up the mesh and an initial guess for the solver.
        xi = np.linspace(0, 1, elements + 1)

        # An initial guess that is physically plausible can significantly help convergence.
        # We expect liquid concentration (x_T) to decrease from inlet (xi=1) to outlet (xi=0).
        # We expect gas concentration (y_T2) to increase from inlet (xi=0) to outlet (xi=1).
        y_guess = np.zeros((4, xi.size))
        y_guess[0, :] = np.linspace(0.5, 1.0, xi.size)  # Guess for x_T (linear decrease)
        y_guess[1, :] = -0.5 # Guess for dx_T/dxi
        y_guess[2, :] = np.linspace(y_T2_in, y_T2_in + 1e-4, xi.size)  # Guess for y_T2 (linear increase)
        y_guess[3, :] = 1e-4 # Guess for dy_T2/dxi

        # Run the BVP solver
        sol = solve_bvp(ode_system, boundary_conditions, xi, y_guess, tol=1e-5, max_nodes=10000)

        return sol





    # --- Unpack User Input Parameters ---
    # Bubble Column parameters
    c_T_inlet = params["c_T_inlet"]  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid just before inlet
    y_T2_in = params["y_T2_in"]  # Inlet tritium molar fraction in gas (0 = pure purge gas)
    P_0 = params["P_0"]  # Pa, Gas total pressure at (gas) inlet)

    L = params["L"]  # m, Height of the bubble column
    D = params["D"]  # m, Column diameter

    Flow_l = params["Flow_l"] # kg/s, Liquid mass flow rate
    Flow_g = params["Flow_g"] # mol/s, Gas molar flow rate
    
    T = params["T"]  # K, Temperature

    # Solver parameters
    elements = params["elements"] # Number of mesh elements for solver

    # --- Constants ---
    g = 9.81  # m/s^2, Gravitational acceleration
    R = const.R  # J/(mol·K), Universal gas constant
    N_A = const.N_A  # 1/mol, Avogadro's number
    M_LiPb = 2.875E-25 # Kg/molecule, Lipb molecular mass

    # --- Calculated parameters ---

    # Calculate empirical correlations
    ρ_l = 10.45e3 * (1 - 1.61e-4 * T)  # kg/m^3, Liquid (LiPb) density
    σ_l = 0.52 - 0.11e-3 * T  # N/m, Surface tension, liquid (LiPb) - gas (He) interface
    μ_l = 1.87e-4 * np.exp(11640 / (R*T))  # Pa.s, Dynamic viscosity of liquid LiPb
    ν_l = μ_l / ρ_l # m^2/s, Kinematic viscosity of liquid LiPb
    D_T = 2.5E-7 * np.exp(-27000/(R*T))  # m^2/s, Tritium diffusion coefficient in liquid LiPb

    K_s = 2.32E-8 * np.exp(-1350 / (R*T))  # atfrac*Pa^0.5, Sievert's constant for tritium in liquid LiPb
    K_s = K_s * (ρ_l / (M_LiPb * N_A))  # mol/(m^3·Pa^0.5)


    # Calculate the volumetric flow rates
    Q_l = Flow_l / ρ_l  # m^3/s, Volumetric flow rate of liquid phase
    Q_g = (Flow_g * R * T) / P_0  # m^3/s, Volumetric flow rate of gas phase at inlet

    # Calculate the superficial flow velocities
    A = np.pi * (D / 2) ** 2  # m^2, Cross-sectional area of the column
    u_g0 = Q_g / A  # m/s, superficial gas inlet velocity
    u_l = Q_l / A  # m/s, superficial liquid inlet velocity

    # Calculate Bond, Galilei, Schmidt and Froude numbers
    Bn = (g * D**2 * ρ_l) / σ_l  # Bond number
    Ga = (g * D**3) / ν_l**2  # Galilei number
    Sc = ν_l / D_T # Schmidt number
    Fr = u_g0 / (g * D)**0.5  # Froude number

    # Calculate dispersion coefficients
    E_l = (D * u_g0) / ((13 * Fr) / (1 + 6.5*(Fr**0.8)) ) # m^2/s, Effective axial dispersion coefficient, liquid phase
    E_g = (0.2 * D**2) * u_g0 # m^2/s, Effective axial dispersion coefficient, gas phase

    # Calculate gas hold-up (phase fraction) & mass transfer coefficient
    C = 0.2 * (Bn**(1/8)) * (Ga**(1/12)) * Fr # C = ε_g / (1 - ε_g)^4

    print(C)

    def solveEqn(ε_g, C):
        # Define the equation to solve
        eqn = ε_g / (1 - ε_g)**4 - C
        return eqn
    ε_g_initial_guess = 0.1
    try:
        # bracket=[0.0001, 0.9999] tells it to *only* look in this range
        sol = root_scalar(solveEqn, args=(C,), bracket=[0.00001, 0.99999])

    except ValueError as e:
        print(f"Solver failed. This can happen if C is so large that no solution exists between 0 and 1.")
        print(f"Error: {e}")

    ε_g = sol.root # Gas phase fraction

    # Calculate outlet pressure hydrostatically & check non-negative
    P_outlet = P_0 - ( ρ_l * (1 - ε_g) * g * L)

    if P_outlet <= 0:
        raise ValueError(
            f"Calculated gas outlet pressure P_outlet must be positive, but got {P_outlet:.2e} Pa. Check P_0, rho_l, g, and L are realistic."
        )

    # Calculate mean bubble diameter
    d_b = (26 * (Bn**-0.5) * (Ga**-0.12) * (Fr**-0.12)) * D  # m, Mean bubble diameter AGREES WITH PAPER

    # Calculate interfacial area
    a = 6 * ε_g / d_b  # m^-1, Specific interfacial area, assuming spherical bubbles

    # Calculate volumetric mass transfer coefficient, liquid-gas
    h_l_a = D_T * (0.6 * Sc**0.5 * Bn**0.62 * Ga**0.31 * ε_g**1.1) / (D**2)  # Volumetric mass transfer coefficient, liquid-gas

    # Calculate mass transfer coefficient
    h_l = h_l_a / a # Mass transfer coefficient

    # Calculate dimensionless values
    ε_l = 1 - ε_g  # Liquid phase fraction

    ψ = (ρ_l * g * ε_l * L) / P_0  # Hydrostatic pressure ratio (Eq. 8.3)
    ν = (c_T_inlet / K_s) ** 2 / P_0  # Tritium partial pressure ratio (Eq. 8.5)
    Bo_l = u_l * L / (ε_l * E_l)  # Bodenstein number, liquid phase (Eq. 8.9)
    phi_l = a * h_l * L / u_l  # Transfer units parameter, liquid phase (Eq. 8.11)
    Bo_g = u_g0 * L / (ε_g * E_g)  # Bodenstein number, gas phase (Eq. 8.10)
    phi_g = (
        0.5 * (R * T * c_T_inlet / P_0) * (a * h_l * L / u_g0)
    )  # Transfer units parameter, gas phase (Eq. 8.12)

    dimensionless_params = {
        "Bo_l": Bo_l,
        "phi_l": phi_l,
        "Bo_g": Bo_g,
        "phi_g": phi_g,
        "psi": ψ,
        "nu": ν,
    }

    # Solve the model
    solution = solve_tritium_extraction(dimensionless_params, y_T2_in, elements)

    # --- Results ---
    if solution.success:
        # --- Dimensionless Results ---
        x_T_outlet_dimless = solution.y[0, 0]
        efficiency = 1 - x_T_outlet_dimless
        y_T2_outlet_gas = solution.y[2, -1]  # y_T2 at xi=1

        # --- Dimensional Results ---
        # Liquid concentration at outlet (xi=0)
        c_T_outlet = x_T_outlet_dimless * c_T_inlet

        # Gas partial pressure at outlet (xi=1)
        P_outlet = P_0 * (
            1 - dimensionless_params["psi"]
        )  # Derived from Eq. 8.4 at xi=1
        P_T2_out = y_T2_outlet_gas * P_outlet

        # Mass transfer consistency check
        N_A = const.N_A  # Avogadro's number, 1/mol
        # Tritium molar flow rate into the column via liquid
        n_T_in_liquid = c_T_inlet * Q_l * N_A  # Triton/s

        # Tritium molar flow rate out of the column via liquid
        n_T_out_liquid = c_T_outlet * Q_l * N_A  # Tritons/s

        # Tritium molar flow rate into the column via gas
        P_T2_in = y_T2_in * P_0 # [Pa]
        n_T2_in_gas = (P_T2_in * Q_g / (R * T)) * N_A  # T2/s
        n_T_in_gas = n_T2_in_gas * 2  # Triton/s

        # Calculate outlet gas volumetric flow rate (gas expands as pressure drops)
        Q_g_out = (P_0 * Q_g) / P_outlet
        # Tritium molar flow rate out of the column via gas
        n_T2_out_gas = (P_T2_out * Q_g_out / (R * T)) * N_A  # T2/s
        n_T_out_gas = n_T2_out_gas * 2  # Triton/s

        T_in = n_T_in_liquid + n_T_in_gas
        T_out = n_T_out_liquid + n_T_out_gas

        results = {
            "Total tritium in [T/s]": T_in,
            "Total tritium out [T/s]": T_out,
            "extraction_efficiency [%]": efficiency * 100,
            "c_T_inlet [mol/m^3]": c_T_inlet,
            "c_T_outlet [mol/m^3]": c_T_outlet,
            "liquid_vol_flow [m^3/s]": Q_l,
            "P_T2_inlet_gas [Pa]": P_T2_in,
            "P_T2_outlet_gas [Pa]": P_T2_out,
            "gas_vol_flow [m^3/s]": Q_g,
            "total_gas_P_outlet [Pa]": P_outlet,
        }

    else:
        print("BVP solver failed to converge.")

    return results
