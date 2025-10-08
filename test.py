import bubble_col_glc
import numpy as np

# These parameters are estimated from the paper's reference case in Table 2
# for L=3m, D=0.25m, and c_T(L+) = 3.09e-2 mol/m^3.

# --- Parameters ---
# Input parameters
c_T_inlet = 3.09e-2  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid just before inlet 
y_T2_in = 0.0        # Inlet tritium molar fraction in gas (0 = pure purge gas)
P_0 = 5e5            # Pa, Gas total pressure at inlet
ρ_l = 9000           # kg/m^3, Liquid density
K_s = 2e-6           # Tritium Sievert's constant in liquid

L = 3.0              # m, Height of the bubble column
D = 0.5              # m, Column diameter
ε_g = 0.04           # Gas phase fraction

Q_l = 0.01           # m^3/s, Volumetric flow rate of liquid phase
Q_g = 0.0005         # m^3/s, Volumetric flow rate of gas phase at inlet

a = 20               # m^-1, Specific liquid-gas interfacial area

E_g = 0.05           # m^2/s, Effective axial dispersion coefficient, gas phase
E_l = 0.01           # m^2/s, Effective axial dispersion coefficient, liquid phase

h_l = 1e-4           # m/s, Mass transfer coefficient, tritium liquid - gas

g = 9.81             # m/s^2, Gravitational acceleration

# Calculated parameters
ε_l = 1 - ε_g        # Liquid phase fraction

A = np.pi * (D/2)**2 # m^2, Cross-sectional area of the column
A_l = A * ε_l        # m^2, Cross-sectional area of the liquid phase
A_g = A * ε_g        # m^2, Cross-sectional area of the gas phase

u_l = Q_l / A_l      # m/s, superficial liquid inlet velocity (positive for downward flow)
u_g0 = Q_g / A_g      # m/s, gas inlet velocity (positive for upward flow)

params = {
    "c_T_inlet": c_T_inlet,
    "y_T2_in": y_T2_in,
    "P_0": P_0,
    "L": L,
    "u_l": u_l,
    "u_g0": u_g0,
    "ε_g": ε_g,
    "ε_l": ε_l,
    "E_g": E_g,
    "E_l": E_l,
    "a": a,
    "h_l": h_l,
    "ρ_l": ρ_l,
    "K_s": K_s,
    "g": g,
    "D": D,
}

result = bubble_col_glc.solve(params)

print(result)
