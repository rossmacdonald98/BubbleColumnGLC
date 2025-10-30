import bubble_col_glc
import numpy as np

# These parameters are estimated from the paper's reference case in Table 2
# for L=3m, D=0.25m, and c_T(L+) = 3.09e-2 mol/m^3.

# --- Parameters ---
# Input parameters
c_T_inlet = 3.09e-2  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid just before inlet
y_T2_in = 0.0  # Inlet tritium molar fraction in gas (0 = pure purge gas)
P_outlet = 5e5  # Pa, Gas total pressure at outlet
ρ_l = 9000  # kg/m^3, Liquid density
K_s = 2e-6  # Tritium Sievert's constant in liquid

L = 100.0  # m, Height of the bubble column
D = 0.5  # m, Column diameter
ε_g = 0.04  # Gas phase fraction

Q_l = 0.01  # m^3/s, Volumetric flow rate of liquid phase
Q_g = 0.0005  # m^3/s, Volumetric flow rate of gas phase at inlet

a = 20  # m^-1, Specific liquid-gas interfacial area

E_g = 0.05  # m^2/s, Effective axial dispersion coefficient, gas phase
E_l = 0.01  # m^2/s, Effective axial dispersion coefficient, liquid phase

h_l = 1e-4  # m/s, Mass transfer coefficient, tritium liquid - gas

g = 9.81  # m/s^2, Gravitational acceleration
T = 630  # K, Temperature

elements = 50  # Number of mesh elements for solver

params = {
    "c_T_inlet": c_T_inlet,
    "y_T2_in": y_T2_in,
    "P_outlet": P_outlet,
    "elements": elements,
    "ρ_l": ρ_l,
    "K_s": K_s,
    "L": L,
    "D": D,
    "ε_g": ε_g,
    "Q_l": Q_l,
    "Q_g": Q_g,
    "a": a,
    "E_g": E_g,
    "E_l": E_l,
    "h_l": h_l,
    "g": g,
    "T": T,
}

result = bubble_col_glc.solve(params)

print(result)
