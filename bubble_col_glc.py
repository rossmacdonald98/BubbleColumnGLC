import numpy as np
from scipy.integrate import solve_bvp


def solve(params):

    def solve_tritium_extraction(dimensionless_params, y_T2_in):
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

            # Dimensionless driving force theta. This is equation (12) in the paper.
            theta = x_T - np.sqrt(np.maximum(0, (1 - psi * xi) * y_T2 / nu))

            # Equation for d(S[0])/d(xi) = d(x_T)/d(xi)
            dS0_dxi = dx_T_dxi

            # Equation for d(S[1])/d(xi) = d^2(x_T)/d(xi)^2 from eq (10)
            dS1_dxi = Bo_l * (phi_l * theta - dx_T_dxi)

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
            res1 = Sa[1]

            # At xi = 1: x_T(1) = 1 - (1/Bo_l) * dx_T/d(xi)|_1
            res2 = Sb[0] - (1 - (1 / Bo_l) * Sb[1])

            # At xi = 0: y_T2(0) = y_T2(0-) + (1/Bo_g) * dy_T2/d(xi)|_0
            res3 = Sa[2] - (y_T2_in + (1 / Bo_g) * Sa[3])

            # At xi = 1: dy_T2/d(xi) = 0
            res4 = Sb[3]

            return np.array([res1, res2, res3, res4])

        # Set up the mesh and an initial guess for the solver.
        xi = np.linspace(0, 1, 51)

        # A flat initial guess is often good enough to get the solver started.
        y_guess = np.zeros((4, xi.size))
        y_guess[0, :] = 0.8  # Guess for x_T
        y_guess[2, :] = 1e-5  # Guess for y_T2

        # Run the BVP solver
        sol = solve_bvp(ode_system, boundary_conditions, xi, y_guess, tol=1e-5)

        return sol

    # Unpack parameters
    c_T_inlet = params["c_T_inlet"]
    y_T2_in = params["y_T2_in"]
    P_0 = params["P_0"]
    L = params["L"]
    u_l = params["u_l"]
    u_g0 = params["u_g0"]
    ε_g = params["ε_g"]
    ε_l = params["ε_l"]
    E_g = params["E_g"]
    E_l = params["E_l"]
    a = params["a"]
    h_l = params["h_l"]
    ρ_l = params["ρ_l"]
    K_s = params["K_s"]
    g = params["g"]
    D = params["D"]

    # Calculate some other parameters

    A = np.pi * (D / 2) ** 2  # m^2, Cross-sectional area of the column
    A_l = A * ε_l  # m^2, Cross-sectional area of the liquid phase
    A_g = A * ε_g  # m^2, Cross-sectional area of the gas phase

    Q_l = u_l * A_l  # m^3/s, Volumetric flow rate of liquid phase

    # Calculate dimensionless values
    # Eq 13 in the paper
    Bo_l = u_l * L / ((1 - ε_g) * E_l)  # Bodenstein number, liquid phase
    Φ_l = a * h_l * L / u_l  # Mass transfer parameter, liquid phase
    # Eq 14 in the paper
    Bo_g = (u_g0 * L / (ε_g * E_g))  # Bodenstein number, gas phase (assumed near plug-flow condition)
    Φ_g = a * h_l * L / u_l  # Mass transfer parameter, gas phase
    # Eq 15 in the paper
    ψ = (ρ_l * g * (1 - ε_g) * L) / P_0  # Hydrostatic pressure ratio
    ν = (c_T_inlet / K_s) ** 2 / P_0  # Tritium partial pressure ratio

    dimensionless_params = {
        "Bo_l": Bo_l,
        "phi_l": Φ_l,
        "Bo_g": Bo_g,
        "phi_g": Φ_g,
        "psi": ψ,
        "nu": ν,
    }

    # Solve the model
    solution = solve_tritium_extraction(dimensionless_params, y_T2_in)

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
        P_outlet_gas = P_0 * (1 - dimensionless_params["psi"])
        P_T2_outlet_gas = y_T2_outlet_gas * P_outlet_gas

        # Tritium extraction rate (mol/s)
        extraction_rate = Q_l * (c_T_inlet - c_T_outlet)  # mol/s

        # Gas velocity at outlet (xi=1)
        u_g_outlet = u_g0 / (1 - dimensionless_params["psi"])

        results = {
            "extraction_efficiency [%]": efficiency * 100,
            "extraction_rate [mol/s]": extraction_rate,
            "c_T_outlet [mol/m^3]": c_T_outlet,
            "liquid_vol_flow [m^3/s]": Q_l,
            "total_gas_P_outlet [Pa]": P_outlet_gas,
            "P_T2_outlet_gas [Pa]": P_T2_outlet_gas,
            "gas_vol_flow_outlet [m^3/s]": u_g_outlet * A_g,

        }

    else:
        print("BVP solver failed to converge.")

    return results
