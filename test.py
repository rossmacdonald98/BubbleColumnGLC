import bubble_col_glc
import numpy as np

# These parameters are estimated from the paper's reference case in Table 2
# for L=3m, D=0.25m, and c_T(L+) = 3.09e-2 mol/m^3.

# --- Parameters ---
# Input parameters
c_T_inlet = 1.96e-2  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid just before inlet
y_T2_in = 1e-20  # Inlet tritium molar fraction in gas (0 = pure purge gas)
P_0 = 5e5  # Pa, Gas total pressure at gas inlet / liquid outlet (bottom of column)
BCs = "O-C" # Boundary conditions type: "O-C" (Open-Closed) or "C-C" (Closed-Closed)
L = 3  # m, Height of the bubble column
D = 0.5  # m, Column diameter
Flow_l = 560  # kg/s, Liquid mass flow rate
Flow_g = 0.0195  # mol/s, Gas molar flow rate
T = 623  # K, Temperature
elements = 50  # Number of mesh elements for solver

params = {
    "c_T_inlet": c_T_inlet,
    "y_T2_in": y_T2_in,
    "P_0": P_0,
    "BCs": BCs,
    "L": L,
    "D": D,
    "Flow_l": Flow_l,
    "Flow_g": Flow_g,
    "T": T,
    "elements": elements,
}

result = bubble_col_glc.solve(params)

for key, value in result.items():
    print(f"{key}: {value}")
