# Mathematical Model of a Gas-Liquid Contactor Bubble Column

This document explains the mathematical equations used to model the behavior of tritium in a gas-liquid contactor bubble column, based on the paper "Tritium Extraction from Pb-17Li by Bubble Columns" by Carlo Malara. The model consists of mass balance equations for both the liquid and gas phases, boundary conditions, and several supporting equations.

The original paper is available at: https://doi.org/10.13182/FST95-A30485

## Nomenclature

* **$a$**: liquid-gas interfacial area per unit of volume [$m^{-1}$]
* **$Bo_g$**: Bodenstein number for the gas phase [-]
* **$Bo_l$**: Bodenstein number for the liquid phase [-]
* **$c_T$**: tritium concentration in liquid phase [$mol \cdot m^{-3}$]
* **$c_{T,eq}$**: tritium concentration in liquid phase at equilibrium with the gas phase [$mol \cdot m^{-3}$]
* **$E_l$**: dispersion coefficient in the liquid phase [$m^2 \cdot s^{-1}$]
* **$E_g$**: dispersion coefficient in the gas phase [$m^2 \cdot s^{-1}$]
* **$g$**: gravitational acceleration [$m \cdot s^{-2}$]
* **$h_l$**: tritium mass transfer coefficient from the liquid to the gas phase [$m \cdot s^{-1}$]
* **$J_T$**: tritium flux from liquid to gas phase [$mol \cdot m^{-2} \cdot s^{-1}$]
* **$K_s$**: tritium Sievert's constant in liquid [$mol \cdot m^{-3} \cdot Pa^{-0.5}$]
* **$L$**: extractor height [m]
* **$P$**: total gas pressure [Pa]
* **$P_0$**: total gas pressure at the extractor inlet [Pa]
* **$R$**: gas constant [$J \cdot mol^{-1} \cdot K^{-1}$]
* **$T$**: temperature [K]
* **$u_g$**: superficial gas velocity [$m \cdot s^{-1}$]
* **$u_{g0}$**: superficial gas velocity at the extractor inlet [$m \cdot s^{-1}$]
* **$u_l$**: superficial liquid velocity [$m \cdot s^{-1}$]
* **$x_T$**: dimensionless liquid concentration [-]
* **$y_{T_2}$**: tritium molar fraction in the gas phase [-]
* **$z$**: axial coordinate in the extractor [m]

## Greek Letters

* **$\epsilon_g$**: gas phase fraction [-]
* **$\epsilon_l$**: liquid phase fraction [-]
* **$\xi$**: dimensionless axial position [-]
* **$\theta$**: dimensionless driving force ratio [-]
* **$\phi_l$**: dimensionless 'number of transfer units' parameter for the liquid phase [-]
* **$\phi_g$**: dimensionless 'number of transfer units' parameter for the gas phase parameter [-]
* **$\psi$**: dimensionless pressure ratio [-]
* **$\nu$**: dimensionless equilibrium ratio [-]
* **$\rho_l$**: Liquid density [$kg \cdot m^{-3}$]

---


## 1. Liquid phase Tritium mass balance

$$
\underbrace{\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}}}_{\text{axial dispersion}} + \underbrace{u_{l}\frac{dc_{T}}{dz}}_{\text{convection}} - \underbrace{J_{T}a}_{\text{mass transfer}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{1.1}
$$
This is equation (1) in the paper.

**Axial dispersion** acts to flatten the tritium concentration profile due to the random net movement of particles from areas of high concentration to low concentration. If the local concentration profile is a valley, the curvature and dispersion terms are positive; if it's a peak, dispersion is negative.

**Convection** either locally accumulates or removes tritium as a result of the fluid velocity and concentration gradient. If the upstream concentration is greater than downstream (positive gradient), there is more tritium flowing in than out and tritium accumulates, and vice-versa. If the concentration profile is flat then the flow in equals the flow out and the net accumulation/removal is 0.

**Mass transfer** is the flux of tritium atoms from the liquid to the gas.


## 2. Gas phase Tritium mass balance

$$
\frac{1}{RT} \left( \underbrace{-\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}}}_{\text{axial dispersion}} + \underbrace{\frac{d(u_{g}Py_{T_{2}})}{dz}}_{\text{advection}} \right) - \underbrace{\frac{1}{2}J_{T}a}_{\text{mass transfer}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{2.1}
$$

This is equation (2) in the paper.
-   The term $\frac{1}{RT}$ is used to convert partial pressures to number of molecules.
-   The factor of $\frac{1}{2}$ is to convert from tritium atoms $T$ to molecules $T_{2}$.

**Axial dispersion** acts to flatten the tritium partial pressure profile due to the random net movement and mixing of gas bubbles from areas of high concentration to low concentration. If the local concentration profile is a valley, the curvature and dispersion terms are positive; if it's a peak, dispersion is negative.

**Advection** transports tritium due to the movement of the gas phase through the liquid. As the gas velocity, pressure, and molar fraction of tritium change with $z$, they are all included in the derivative.

**Mass transfer** is the flux of $T_{2}$ molecules from the liquid to the gas.


## 3. Liquid Phase Boundary Conditions
This is equation (3) in the paper.
-   **At liquid outlet ($z=0$):** Tritium concentration profile is flat. Only convection, no dispersion.
    $$
    \frac{dc_{T}}{dz}\bigg|_{z=0} = 0 \quad [mol \cdot m^{-4}] \tag{3.1}
    $$

-   **At liquid inlet ($z=L$):** Tritium entering the column via convection is equal to the sum entering/leaving via dispersion.
    $$
    \underbrace{u_{l}[c_{T}(L) - c_{T}(L^{+})]}_{\text{convective flux}} = \underbrace{-\epsilon_{l}E_{l}\frac{dc_{T}}{dz}\bigg|_{z=L}}_{\text{dispersive flux}} \quad [mol \cdot m^{-2} \cdot s^{-1}] \tag{3.2}
    $$


## 4. Gas Phase Boundary Conditions
This is equation (4) in the paper.
-   **At gas inlet ($z=0$):** Tritium entering the column via advection is equal to the sum entering/leaving via dispersion.
    $$
    \underbrace{u_{g}[y_{T_{2}}(0^{-}) - y_{T_{2}}(0)]}_{\text{advective flux}} = \underbrace{-\epsilon_{g}E_{g}\frac{dy_{T_{2}}}{dz}\bigg|_{z=0}}_{\text{dispersive flux}} \quad [m \cdot s^{-1}] \tag{4.1}
    $$

    The units for this equation are $m \cdot s^{-1}$ because it represents a balance of species velocities, not molar fluxes. On the left, the gas velocity `u_g` ($m \cdot s^{-1}$) is multiplied by the dimensionless molar fraction `y_T2`. On the right, the dispersion coefficient `E_g` ($m^{2} \cdot s^{-1}$) is multiplied by the molar fraction gradient `dy/dz` ($m^{-1}$). Both sides resolve to $m \cdot s^{-1}$.

-   **At gas outlet ($z=L$):** Tritium concentration profile is flat. Only advection, no dispersion.
    $$
    \frac{dy_{T_{2}}}{dz}\bigg|_{z=L} = 0 \quad [m^{-1}] \tag{4.2}
    $$


## 5. Gas Phase Total Mass Balance
This is equation (5) in the paper.

Enforces that the only change in the total gas mass flow along the length of the column is due to the source term of tritium flux from the liquid to the gas.

$$
\underbrace{\frac{1}{RT}\frac{d(Pu_{g})}{dz}}_{\text{Change in total gas flow}} - \underbrace{\frac{1}{2}J_{T}a}_{\text{Source term}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{5.1}
$$


## 6. Gas Phase Velocity

The pressure in the column is assumed to be purely hydrostatic, based on the height $z$ and pressure $P_{0}$ at $z=0$.

$$
P = P_{0} - \rho_{l}g(1-\epsilon_{g})z \quad [Pa] \tag{6.1}
$$
This is equation (6) in the paper.

By neglecting the change in mass flow due to the source term (assuming it is negligible), and assuming the gas phase fraction remains constant along the length of the column, the equation for the gas velocity profile can be obtained:

$$
\frac{d(Pu_{g})}{dz} = 0 \quad \rightarrow \quad Pu_{g} = \text{constant} \quad [Pa \cdot m \cdot s^{-1}] \tag{6.2}
$$

$$
P(z) \cdot u_{g}(z) = P_{0} \cdot u_{g0}   \quad \rightarrow \quad    u_{g}(z) = \frac{P_{0} \cdot u_{g0}}{P(z)} \quad [m \cdot s^{-1}] \tag{6.3}
$$

Substituting (6.1) in (6.3) gives:

$$
u_{g}(z) = u_{g0}\frac{P_{0}}{P_{0} - \rho_{l}g(1-\epsilon_{g})z} \quad [m \cdot s^{-1}] \tag{6.4}
$$

This is equation (7) in the paper.


## 7. Liquid – Gas Flux

A simple linear driving force model is used to specify the rate of tritium transport across the liquid-gas interface. This is some mass transfer coefficient $h_{l}$ multiplied by the difference between the liquid concentration and concentration when it is in equilibrium with the gas phase partial pressure.

From Sievert's law:

$$
c_{T,eq} = K_{S}(P_{T_{2}})^{0.5} \quad [mol \cdot m^{-3}] \tag{7.1}
$$

The flux $J_{T}$ is:

$$
J_{T}(z) = h_{l}(c_{T} - K_{s}(Py_{T_{2}})^{0.5}) \quad [mol \cdot m^{-2} \cdot s^{-1}] \tag{7.2}
$$

This is equation (8) in the paper.


## 8. Dimensionless Equations

Nondimensionalization of the equations makes the model easier to solve. The following nondimensional variables are used:

-   **Axial position $\xi$** from $z=0$ to $z=L$, equation (9) in the paper:
    $$
    \xi = z/L \tag{8.1}
    $$

-   **Liquid tritium concentration $x_{T}$** scaled relative to the inlet concentration $c_{T}(L^{+})$, equation (9) in the paper:
    $$
    x_{T} = c_{T}/c_{T}(L^{+}) \tag{8.2}
    $$

-   **Pressure ratio $\psi$** is the ratio of the hydrostatic pressure head to the gas inlet pressure $P_{0}$:
    $$
    \psi = \frac{\rho_{l}g(1-\epsilon_{g})L}{P_{0}} \tag{8.3}
    $$
    This is equation (15) in the paper

    By substituting `z = ξL` (from eq. 8.1) into the dimensional pressure equation (6.1) and then factoring out `P_0`:
    $$
    P = P_{0} - \rho_{l}g(1-\epsilon_{g})\xi L = P_{0} \left( 1 - \frac{\rho_{l}g(1-\epsilon_{g})L}{P_0}\xi \right)
    $$
    The fractional term is the definition of `ψ` (eq. 8.3), which simplifies the expression to:

    $$
    P = P_{0}(1-\xi\psi) \tag{8.4}
    $$


-   **Equilibrium ratio $\nu$** is the ratio of tritium partial pressure that would be in equilibrium with the inlet liquid to the total gas pressure at the inlet, which is equation (15) in the paper:
    $$
    v = \frac{[c_{T}(L^{+})/K_{s}]^{2}}{P_{0}} \tag{8.5}
    $$

-   **Dimensionless driving force $\theta$**: The liquid-gas tritium flux $J_{T}$ is driven by the term $(c_{T} - c_{T,eq})$. It is made dimensionless by dividing by the liquid concentration at the inlet, $c_{T}(L^{+})$:
    $$
    \frac{(c_{T}-c_{T,eq})}{c_{T}(L^{+})} = \frac{c_{T}}{c_{T}(L^{+})} - \frac{c_{T,eq}}{c_{T}(L^{+})} = x_{T} - \frac{K_{s}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})} \tag{8.6}
    $$
    Bringing all terms inside the square root:
    $$
    \frac{K_{S}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})} \quad \rightarrow \quad \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5} \tag{8.7}
    $$
    Using the definition for $v$ gives $\frac{K_{s}^{2}}{(c_{T}(L^{+}))^{2}} = \frac{1}{vP_{0}}$, and with $P = P_{0}(1-\xi\psi)$:
    $$
    \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5} = \left[\left(\frac{1}{vP_{0}}\right) \cdot P_{0}(1-\xi\psi) \cdot y_{T_{2}}\right]^{0.5} = \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5}
    $$
    Finally, this gives the dimensionless driving force $\theta$ as:
    $$
    \theta = x_{T} - \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5} \tag{8.8}
    $$
    This is equation (12) in the paper.

-   **Bodenstein numbers**: $Bo_{l}$ and $Bo_{g}$ represent the ratio of advective (bulk flow) transport to dispersive (axial mixing) transport.
    $$
    Bo_{l} = \frac{u_{l}L}{(1-\epsilon_{g})E_{l}} \tag{8.9}
    $$
    Equation (13) in the paper.
    $$
    Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}} \tag{8.10}
    $$
    Equation (14) in the paper.
    (Note: The original paper has $u_g$, but it should be $u_{g0}$ to be consistent with the dimensionless group definition).

-   **Number of Transfer Units**: $\phi_{l}$ and $\phi_{g}$ compare the residence time of the fluid in the column to the characteristic time required for mass transfer.
    $$
    \phi_{l} = \frac{ah_{l}L}{u_{l}} \tag{8.11}
    $$
    Equation (13) in the paper.
    $$
    \phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \cdot \frac{ah_{l}L}{u_{g0}} \tag{8.12}
    $$
    Equation (14) in the paper.

## 9. Dimensionless ODE System

To solve the model, the original partial differential equations are transformed into a system of dimensionless ordinary differential equations (ODEs) using the variables and parameters defined in Section 8. This process simplifies the equations and groups physical parameters into meaningful dimensionless numbers.

### 9.1. Liquid Phase ODE

We start with the liquid phase mass balance (Equation 1.1):
$$
\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}} + u_{l}\frac{dc_{T}}{dz} - J_{T}a = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{1.1}
$$
First, we substitute the dimensionless variables for concentration $c_T = x_T c_T(L^{+})$ and axial position $z = \xi L$:
$$
\frac{d}{dz} = \frac{1}{L}\frac{d}{d\xi} \quad \text{and} \quad \frac{d^2}{dz^2} = \frac{1}{L^2}\frac{d^2}{d\xi^2}
$$
Equation 1.1 becomes:
$$
\frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}(x_T c_T(L^{+}))}{d\xi^{2}} + \frac{u_{l}}{L}\frac{d(x_T c_T(L^{+}))}{d\xi} - J_{T}a = 0 \tag{9.1.1}
$$
Since $c_T(L^{+})$ is a constant, it can be factored out. We also substitute the expression for the flux, $J_T = h_l(c_T - c_{T,eq})$:
$$
c_T(L^{+})\left( \frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}x_T}{d\xi^{2}} + \frac{u_{l}}{L}\frac{dx_T}{d\xi} \right) - ah_l(c_T - c_{T,eq}) = 0 \tag{9.1.2}
$$
Now, we divide the entire equation by $\frac{u_l c_T(L^{+})}{L}$:
$$
\frac{\epsilon_{l}E_{l}}{u_l L}\frac{d^{2}x_T}{d\xi^{2}} + \frac{dx_T}{d\xi} - \frac{ah_l L}{u_l} \frac{(c_T - c_{T,eq})}{c_T(L^{+})} = 0 \tag{9.1.3}
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the liquid phase: $Bo_{l} = \frac{u_{l}L}{\epsilon_{l}E_{l}}$ (assuming $\epsilon_l = 1-\epsilon_g$)
-   Number of transfer units for the liquid phase: $\phi_{l} = \frac{ah_{l}L}{u_{l}}$
-   Dimensionless driving force: $\theta = \frac{c_T - c_{T,eq}}{c_T(L^{+})}$

We arrive at the final dimensionless ODE for the liquid phase:
$$
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0 \tag{9.1.4}
$$
This is equation (10) in the paper.

### 9.2. Gas Phase ODE

The derivation for the gas phase is more involved. We start with the gas phase mass balance (Equation 2.1):
$$
\frac{1}{RT} \left( -\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}} + \frac{d(u_{g}Py_{T_{2}})}{dz} \right) - \frac{1}{2}J_{T}a = 0 \tag{2.1}
$$
We substitute the dimensionless variables and expressions for $P(\xi)$ and $u_g(\xi)$:
-   $P(\xi) = P_0(1-\xi\psi)$
-   $u_g(\xi) = \frac{u_{g0}}{1-\xi\psi}$
-   $J_T = h_l c_T(L^{+}) \theta$

The advection term $u_g P y_{T_2}$ simplifies nicely:
$$
u_g P y_{T_2} = \left(\frac{u_{g0}}{1-\xi\psi}\right) \left(P_0(1-\xi\psi)\right) y_{T_2} = u_{g0}P_0 y_{T_2} \tag{9.2.1}
$$
The derivative of the advection term with respect to $z$ is:
$$
\frac{d(u_{g}Py_{T_{2}})}{dz} = \frac{1}{L}\frac{d(u_{g0}P_0 y_{T_2})}{d\xi} = \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi} \tag{9.2.2}
$$

This simplification is the primary reason for using the molar fraction $y_{T_{2}}$ as the state variable for the gas phase instead of the partial pressure $P_{T_{2}}$. By choosing $y_{T_{2}}$, the position-dependent terms for gas velocity $u_g(z)$ and total pressure $P(z)$ cancel out in the advection term, leading to a much simpler derivative.

For the dispersion term, we must differentiate the full product $P(\xi)y_{T_2}(\xi)$ twice with respect to $\xi$. Let's denote $y_{T_2}$ as $y$ for brevity.
The first derivative is:
$$
\frac{d(Py)}{d\xi} = \frac{d}{d\xi}[P_0(1-\xi\psi)y] = P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right] \tag{9.2.3}
$$
The second derivative is:
$$
\frac{d^2(Py)}{d\xi^2} = \frac{d}{d\xi} \left( P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right] \right) = P_0 \left[ (1-\xi\psi)\frac{d^2y}{d\xi^2} - 2\psi \frac{dy}{d\xi} \right] \tag{9.2.4}
$$
The full dispersion term in equation 2.1 is therefore:
$$
-\epsilon_{g}E_{g}\frac{d^{2}(Py_{T_{2}})}{dz^{2}} = -\frac{\epsilon_{g}E_{g}}{L^2} \frac{d^2(Py_{T_2})}{d\xi^2} = -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] \tag{9.2.5}
$$
Substituting the dimensionless dispersion and advection terms back into the gas phase balance:
$$
-\frac{1}{RT} \left( -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] + \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi} \right) - \frac{1}{2}ah_l c_T(L^{+}) \theta = 0 \tag{9.2.6}
$$
Now, we divide the entire equation by $\frac{u_{g0}P_0}{RTL}$:
$$
\frac{\epsilon_{g}E_{g}}{u_{g0}L} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] - \frac{dy_{T_2}}{d\xi} + \frac{RTL}{2u_{g0}P_0} ah_l c_T(L^{+}) \theta = 0 \tag{9.2.7}
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the gas phase: $Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}}$
-   Number of transfer units for the gas phase: $\phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \frac{ah_{l}L}{u_{g0}}$

Substituting these groups into equation 9.2.7 gives:
$$
\frac{1}{Bo_g} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] - \frac{dy_{T_2}}{d\xi} + \phi_g \theta = 0 \tag{9.2.8}
$$
Next, we expand the first term by distributing the $\frac{1}{Bo_g}$ factor:
$$
\frac{(1-\xi\psi)}{Bo_g}\frac{d^2y_{T_2}}{d\xi^2} - \frac{2\psi}{Bo_g}\frac{dy_{T_2}}{d\xi} - \frac{dy_{T_2}}{d\xi} + \phi_g \theta = 0 \tag{9.2.9}
$$
Finally, we group the terms containing the first derivative $\frac{dy_{T_2}}{d\xi}$:

We arrive at the final dimensionless ODE for the gas phase:
$$
\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} - \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} + \phi_{g}\theta = 0 \tag{9.2.10}
$$
In the paper they write this equation with opposite signs:
$$
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0 \tag{9.2.11}
$$
This is equation (11) in the paper.
### 9.3. Summary of the System

The complete model is a system of two second-order ODEs for the dimensionless concentrations $x_T$ and $y_{T_2}$ as a function of the dimensionless position $\xi$.

The two governing equations are the liquid phase balance (9.1.4) and the gas phase balance (9.2.10):
$$
\begin{cases}
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0 \quad \text{(Liquid Phase)} \\
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0 \quad \text{(Gas Phase)}
\end{cases} \tag{9.3.1}
$$
Where:
$$
\theta = x_{T} - \sqrt{\frac{(1-\xi\psi)y_{T_{2}}}{\nu}} \tag{8.8}
$$

This system is:
-   **Coupled:** The two equations are linked through the dimensionless driving force, $\theta$, which depends on both $x_T$ and $y_{T_2}$.
-   **Non-linear:** The presence of the square root in the definition of $\theta$ makes the system non-linear.
-   A **Boundary Value Problem (BVP):** The conditions are specified at the two ends of the domain ($\xi=0$ and $\xi=1$), not just at the start.

To solve this system numerically, it is standard practice to convert the two second-order ODEs into a system of four first-order ODEs. We define a state vector $S$ as:
$$
S = \begin{bmatrix} S_0 \\ S_1 \\ S_2 \\ S_3 \end{bmatrix} = \begin{bmatrix} x_T \\ dx_T/d\xi \\ y_{T_2} \\ dy_{T_2}/d\xi \end{bmatrix}
$$
The system can then be written as $\frac{dS}{d\xi}$ by rearranging the equations in (9.3.1) to solve for the second derivatives:
$$
\frac{dS}{d\xi} =
\begin{cases}
dS_0/d\xi = S_1 \\
dS_1/d\xi = Bo_{l} (\phi_{l}\theta - S_1) \\
dS_2/d\xi = S_3 \\
dS_3/d\xi = \frac{Bo_{g}}{1-\xi\psi} \left[ \left(1 + \frac{2\psi}{Bo_{g}}\right)S_3 - \phi_{g}\theta \right]
\end{cases} \tag{9.3.3}
$$
This system of four first-order ODEs, along with the four boundary conditions from Section 10, can be solved using a numerical BVP solver, such as `scipy.integrate.solve_bvp` in Python.

## 10. Dimensionless Boundary Conditions

-   **Liquid outlet ($\xi=0$):** Zero-gradient condition states that tritium leaves the bottom of the column only through the bulk flow (advection). This is equation (16) in the paper.
    $$
    \frac{dx_{T}}{d\xi}\bigg|_{\xi=0} = 0 \tag{10.1}
    $$

-   **Liquid inlet ($\xi=1$):** The concentration of the liquid feed before it enters is different from the concentration just inside the column because of back-mixing. The size of this jump is determined by the Bodenstein number $Bo_{l}$. If there were no dispersion ($Bo_{l} \rightarrow \infty$), the concentration profile would be continuous. This is equation (16) in the paper.
    $$
    x_{T}(1) = 1 - \frac{1}{Bo_{l}}\frac{dx_{T}}{d\xi}\bigg|_{\xi=1} \tag{10.2}
    $$

-   **Gas inlet ($\xi=0$):** The concentration of the incoming gas is immediately altered to a slightly higher value just inside the column due to tritium from higher up mixing back down via dispersion. This is equation (17) in the paper.
    $$
    y_{T_{2}}(0) = y_{T_{2}}(0^{-}) + \frac{1}{Bo_{g}}\frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=0} \tag{10.3}
    $$

-   **Gas outlet ($\xi=1$):** Zero-gradient condition enforces that tritium leaves the top of the column only through the bulk upward flow of the gas bubbles (advection). This is equation (17) in the paper.
    $$
    \frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=1} = 0 \tag{10.4}
    $$


## 11. Parameter Acquisition

The core model requires several parameters that describe the hydrodynamic and mass transfer behaviour of the bubble column (e.g., $\epsilon_g, a, h_l, E_l, E_g$). These parameters are not fundamental material properties but depend on the column's geometry, operating conditions, and the physical properties of the fluids. They are typically estimated using well-established empirical correlations derived from experimental data.

### 11.1. Bond, Galilei, Schmidt & Froude numbers ($Bn$, $Ga$, $Sc$, & $Fr$)

$$
Fr = \frac{u_{g0}}{\sqrt{gD}}
$$

$$
Bn = \frac{g D^{2} \rho_{l}}{\sigma_{l}}
$$
$$
Ga = \frac{g D^{3}}{\nu_{l}^{2}}
$$
$$
Sc = \frac{\nu_{l}}{D_{T,l}}
$$

### 11.2. Liquid Phase Dispersion Coefficient ($E_l$)

Correlation for the liquid-phase dispersion coefficient in terms of the column diameter and liquid properties:

$$
E_{l} = \frac{D \cdot u_{g0}}{((13Fr)/(1+6.5 \cdot (Fr)^{0.8}))} 
$$

### 11.3. Gas Phase Dispersion Coefficient ($E_g$)

A commonly used correlation for the gas-phase dispersion coefficient:

$$
E_{g} = 0.2\,D^{2}u_{g0} \tag{19}
$$


### 11.4. Gas Holdup ($\epsilon_g$)

A correlation that has been used to estimate the gas holdup in similar systems:

$$
\frac{\epsilon_{g}}{(1 - \epsilon_{g})^{4}} = 0.2\,Bn^{1/8}\,Ga^{1/12}\,Fr \tag{21}
$$

(Iterate to solve for $\epsilon_g$.)

### 11.5. Liquid Phase Fraction ($\epsilon_l$)

The liquid phase fraction is simply related to the gas phase fraction:

$$
\epsilon_{l} = 1 - \epsilon_{g} \tag{22}
$$

### 11.6. Mean Bubble Diameter ($d_b$)

$$
d_b = (26\,Bn^{-0.5}\,Ga^{-0.12}\,Fr^{-0.12}) \cdot D 
$$

### 11.7. Specific Interfacial Area ($a$)

????????????????????

### 11.8. Volumetric Mass Transfer Coefficient ($a*h_l$)

An empirical correlation for the volumetric mass transfer coefficient used in the reference paper:

$$
\frac{a h_{l} D^{2}}{D_{T,l}} = 0.6\,Sc^{0.5}\,Bn^{0.62}\,Ga^{0.31}\,\epsilon_{g}^{1.1}
$$
$$
a h_{l}= \frac{(0.6\,Sc^{0.5}\,Bn^{0.62}\,Ga^{0.31}\,\epsilon_{g}^{1.1}) \cdot D_{T,l}}{D^2}
$$

### 11.9. Mass Transfer Coefficient ($h_l$)

When the interfacial area $a$ is known, the mass transfer coefficient can be calculated from the volumetric mass transfer coefficient:

$$
h_l = \frac{a h_{l}}{a}
$$


