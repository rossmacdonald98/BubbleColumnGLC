"""
Bubble Column Gas-Liquid Contactor Model Runner

This script configures and runs a bubble column gas-liquid contactor simulation
using specified operating parameters and boundary conditions.
"""

import bubble_col_glc

def main():
    # Operating Parameters
    # ------------------
    # Input parameters for liquid phase
    c_T_inlet = 1.96e-2  # mol/m^3 (c_T(L+)), Inlet tritium concentration in liquid
    
    # Input parameters for gas phase
    y_T2_in = 0  # Inlet tritium molar fraction in gas (0 = pure purge gas)
    
    # Avoid model instability with zero concentration
    if y_T2_in == 0:
        y_T2_in = 1e-20 
    
    # Physical parameters
    P_0 = 5e5        # Pa, Total pressure at gas inlet / liquid outlet
    L = 3            # m, Height of the bubble column
    D = 0.5          # m, Column diameter
    Flow_l = 560     # kg/s, Liquid mass flow rate
    Flow_g = 0.19    # mol/s, Gas molar flow rate
    T = 623          # K, Temperature
    
    # Solver parameters
    BCs = "O-C"      # Boundary conditions: "O-C" (Open-Closed) or "C-C" (Closed-Closed)
    elements = 50    # Number of initial mesh elements for solver
    
    # Package parameters for solver
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
    
    # Run simulation
    result = bubble_col_glc.solve(params)
    
    # Print results
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
