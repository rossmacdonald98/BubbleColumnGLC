import bubble_col_glc

# These parameters are estimated from the paper's reference case in Table 2
# for L=3m, D=0.25m, and c_T(L+) = 3.09e-2 mol/m^3.

# --- Parameters ---
c_T_inlet = 3.09e-2  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid just before inlet
y_T2_in = 0.0  # Inlet tritium molar fraction in gas (pure purge gas)
P_0 = 5e5  # Pa, Gas total pressure at inlet
L = 3.0  # m, Height of the bubble column
u_l = 0.02  # m/s, superficial liquid inlet velocity (positive for downward flow)
u_g = 0.2  # m/s, gas inlet velocity (positive for upward flow)
ε_g = 0.04  # Gas phase fraction
ε_l = 1 - ε_g  # Liquid phase fraction
E_g = 0.05  # m^2/s, Effective axial dispersion coefficient, gas phase
E_l = 0.01  # m^2/s, Effective axial dispersion coefficient, liquid phase
a = 20  # m^-1, Specific liquid-gas interfacial area
h_l = 1e-4  # m/s, Mass transfer coefficient, tritium liquid - gas
ρ_l = 9000  # kg/m^3, Liquid density
K_s = 1e-6  # Tritium Sievert's constant in liquid
g = 9.81  # m/s^2, Gravitational acceleration

D = 0.25  # m, Column diameter

params = {
    "c_T_inlet": c_T_inlet,
    "y_T2_in": y_T2_in,
    "P_0": P_0,
    "L": L,
    "u_l": u_l,
    "u_g": u_g,
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
