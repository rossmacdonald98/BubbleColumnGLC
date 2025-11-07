# Mathematical Model of a Gas-Liquid Contactor Bubble Column Reactor

This document explains the mathematical equations used to model the behavior of a gas-liquid contactor bubble column reactor (BCR) for extraction of Tritium from liquid Pb-17Li.

This model was originally defined in the 1995 paper by Carlos Malara "Tritium Extraction from Pb-17Li by Bubble Columns". The paper is available at: https://doi.org/10.13182/FST95-A30485

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

## Model Assumptions

* The BCR is a vertically sparged, single stage reactor.
* The liquid phase enters the BCR from the top of the column, and the gas phase flows in a countercurrent direction from the bottom.
* The bubble shape is spherical. 

## Model Theory
We need to simulate the transport of tritium in both the liquid and gas phases, and the transfer between phases. This includes:
* Advective and Dispersive transport in the gas phase.
* Advective and Dispersive transport in the liquid phase.
* Mass transfer between the liquid and gas phases.

### Advection
Advection is the transport of a solute due to the bulk movement of a fluid. 

The rate of change of Tritium concentration due to advection in the liquid phase at a position $z$ in a BCR is:

$$ \frac{dc_{T}}{dt}\bigg|_{z, Advection} = u_{l}\frac{dc_{T}}{dz} $$

This implys that if the concentration gradient and liquid velocity are positive (recalling that the positive velocity direction for the liquid is downwards), advection will act to increase the local concentration as the upstream concentration is higher, and vice-versa.

For the gas phase the direction of flow is bottom to top, so a positive partial pressure gradient and gas velocity will cause advection to reduce the local partial pressure of Tritium:

$$ \frac{dPy_{T_{2}}}{dt}\bigg|_{z, Advection} = -u_{g}\frac{d(u_{g}Py_{T_{2}})}{dz} $$

### Dispersion
Dispersion is the spreading or mixing that occurs due to a combination of molecular diffusion and non-ideal flow patterns.
It acts to flatten the concentration profile because it moves molecules away from areas of high concentration and into areas of low concentration.

If the local concentration profile is a valley, the curvature and dispersion terms are positive; if it's a peak, dispersion is negative.

The rate of change of Tritium concentration due to dispersion in the liquid phase at a position $z$ in a BCR is:

$$ \frac{dc_{T}}{dt}\bigg|_{z, Dispersion} = \epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}} $$

Similarly, the rate of change of Tritium partial pressure due to dispersion in the gas phase at a position $z$ in a BCR is:

$$ \frac{dPy_{T_{2}}}{dt}\bigg|_{z, Dispersion} = \epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}} $$


### Mass Transfer
Gas dissolved in the liquid phase can transfer into the gas phase by moving across the liquid-gas interface.

For a tritium & liquid metal system, the mechanisms involved in mass transfer over the liquid-gas interface include [A]:
*  Diffusion through the liquid boundary layer.
*  Recombination/dissociation reactions at the liquid-gas interface.
*  Diffusion through the gas boundary layer.

Generally the rate limiting step is transfer through the liquid boundary layer. Therefore, the liquid-gas interface transfer flux is defined as:

$$
J_{T} = h_{l}(c_{T} - c_{T}^{eq})
$$

Where $c_{T}$ is the bulk liquid phase tritium concentration and $c_{T}^{eq}$ is the concentration at the interface that is in equilibrium with the gas phase.

The equilbrium concentration $c_{T}^{eq}$ is defined by Sievert's law:

$$
c_{T}^{eq} = K_{S}(Py_{T_{2}})^{0.5}
$$

Where $K_{S}$ is the Sievert's constant for Tritium solubility in the metal, $Py_{T_{2}}$ is the partial pressure of tritium in the gas phase, and $h_l$ is the liquid-gas mass transfer coefficient.



---


## 1. Liquid phase Tritium mass balance

The mass of Tritium must be conserved over the entire column. Therefore we can combine the dispersion, advection and mass transfer into a mass balance equation:

$$
\underbrace{\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}}}_{\text{dispersion}} + \underbrace{u_{l}\frac{dc_{T}}{dz}}_{\text{advection}} - \underbrace{J_{T}a}_{\text{mass transfer}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{1.1}
$$
This is saying that the rate of Tritium mass transfer from the liquid to the gas is balanced by the advective and dispersive transport within the liquid phase. This is equation (1) in the paper.



## 2. Gas phase Tritium mass balance

$$
\frac{1}{RT} \left( \underbrace{\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}}}_{\text{dispersion}} - \underbrace{\frac{d(u_{g}Py_{T_{2}})}{dz}}_{\text{advection}} \right) + \underbrace{\frac{1}{2}J_{T}a}_{\text{mass transfer}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}]
$$
* The term $\frac{1}{RT}$ is used to convert partial pressures to number of molecules.
* The factor of $\frac{1}{2}$ is to convert from tritium atoms $T$ to molecules $T_{2}$.

This is saying that the rate of mass transfer of $T_2$ molecules from the liquid to the gas is balanced by the advective and dispersive transport within the gas phase.

In the paper they use opposite signs:
$$
\frac{1}{RT} \left(-\underbrace{\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}}}_{\text{dispersion}} + \underbrace{\frac{d(u_{g}Py_{T_{2}})}{dz}}_{\text{advection}} \right) - \underbrace{\frac{1}{2}J_{T}a}_{\text{mass transfer}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{2.1}
$$

The $1/RT$ term converts from partial pressures to number of molecules, and the $1/2$ converts from Tritium atoms to $T_2$ molecules. This is equation (2) in the paper.


## 3. Boundary Conditions
There are different types of BCR design that impose different boundary conditions. The original paper defines a "closed-closed" set of boundary conditions, and later work [B] defined an "open-closed" condition. 

The "open" and "closed" terminology for boundary conditions in BCR models has different physical meanings at the inlet versus the outlet:
* At the Inlet:
    * Closed: Models a "messy" boundary with back-mixing, causing a concentration discontinuity (jump).
    * Open: Models a "simple" boundary, ignoring back-mixing and assuming a continuous concentration profile.
* At the Outlet:
    * Closed: Models a "simple" boundary where back-mixing stops, defined by a continuous concentration profile.
    * Open: Models a "messy" boundary where the internal back-mixing process continues into the exit stream.

### 'Closed' BCs

-   **'Closed' liquid inlet ($z=L$):** Tritium entering the column via advection is equal to the sum entering & leaving via dispersion.
    $$
    \underbrace{u_{l}[c_{T}(L) - c_{T}(L^{+})]}_{\text{convective flux}} = \underbrace{-\epsilon_{l}E_{l}\frac{dc_{T}}{dz}\bigg|_{z=L}}_{\text{dispersive flux}} \tag{3.1}
    $$
-   **'Closed' gas inlet ($z=0$):** Tritium entering the column via advection is equal to the sum entering & leaving via dispersion.
    $$
    \underbrace{u_{g}[y_{T_{2}}(0^{-}) - y_{T_{2}}(0)]}_{\text{advective flux}} = \underbrace{-\epsilon_{g}E_{g}\frac{dy_{T_{2}}}{dz}\bigg|_{z=0}}_{\text{dispersive flux}} \tag{3.2}
    $$
-   **'Closed' liquid outlet ($z=0$):** Tritium concentration profile is flat. Only advection, no dispersion.
    $$
    \frac{dc_{T}}{dz}\bigg|_{z=0} = 0  \tag{3.3}
    $$
-   **'Closed' gas outlet ($z=L$):** Tritium partial pressure profile is flat. Only advection, no dispersion.
    $$
    \frac{dy_{T_{2}}}{dz}\bigg|_{z=L} = 0 \tag{3.4}
    $$


### 'Open' BCs

-   **'Open' liquid inlet ($z=L$):** Tritium concentration profile is flat. Only advection, no dispersion. Concentration at inlet is equal to concentration just upstream from inlet.
    $$
    c_{T}(L^+) = c_{T}(L)  \tag{3.5}
    $$

-   **'Open' gas inlet ($z=0$):** Tritium partial pressure profile is flat. Only advection, no dispersion. Partial pressure at inlet is equal to partial pressure just upstream from inlet.
    $$
    y_{T_2}(0) = y_{T_2}(0^-)  \tag{3.6}
    $$

This model currently implements two options for boundary conditions: 
-   **"Closed-Closed":** Closed inlets and closed outlets, as in the Malara paper [A].
-   **"Open-Closed":** Open inlets and closed outlets, as in the Mohan paper [B].


## 4. Gas Phase Total Mass Balance
This is equation (5) in the Malara paper.

Enforces that the only change in the total gas mass flow along the length of the column is due to the source term of tritium flux from the liquid to the gas.

$$
\underbrace{\frac{1}{RT}\frac{d(Pu_{g})}{dz}}_{\text{Change in total gas flow}} - \underbrace{\frac{1}{2}J_{T}a}_{\text{Source term}} = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}] \tag{4.1}
$$


## 5. Gas Phase Velocity

The pressure in the column is assumed to be purely hydrostatic, based on the height $z$ and pressure $P_{0}$ at $z=0$.

$$
P = P_{0} - \rho_{l}g(1-\epsilon_{g})z \quad [Pa] \tag{5.1}
$$
This is equation (6) in the paper.

By neglecting the change in mass flow due to the source term (assuming it is negligible), and assuming the gas phase fraction remains constant along the length of the column, the equation for the gas velocity profile can be obtained:

$$
\frac{d(Pu_{g})}{dz} = 0 \quad \rightarrow \quad Pu_{g} = \text{constant}  \tag{5.2}
$$

$$
P(z) \cdot u_{g}(z) = P_{0} \cdot u_{g0}   \quad \rightarrow \quad    u_{g}(z) = \frac{P_{0} \cdot u_{g0}}{P(z)} \quad [m \cdot s^{-1}] \tag{5.3}
$$

Substituting (5.1) in (5.3) gives:

$$
u_{g}(z) = u_{g0}\frac{P_{0}}{P_{0} - \rho_{l}g(1-\epsilon_{g})z} \quad [m \cdot s^{-1}] \tag{5.4}
$$

This is equation (7) in the paper.


## 6. Liquid – Gas Mass Transfer Flux

A simple linear driving force model is used to specify the rate of tritium transport across the liquid-gas interface. This is some mass transfer coefficient $h_{l}$ multiplied by the difference between the liquid concentration and concentration when it is in equilibrium with the gas phase partial pressure.

From Sievert's law:

$$
c_{T,eq} = K_{S}(P_{T_{2}})^{0.5} \quad [mol \cdot m^{-3}] \tag{6.1}
$$

The flux $J_{T}$ is:

$$
J_{T} = h_{l}(c_{T} - K_{s}(Py_{T_{2}})^{0.5}) \quad [mol \cdot m^{-2} \cdot s^{-1}] \tag{6.2}
$$

This is equation (8) in the paper.


## 7. Dimensionless Equations

Nondimensionalization of the equations makes the model easier to solve. The following nondimensional variables are used:

-   **Axial position $\xi$** from $z=0$ to $z=L$, equation (9) in the paper:
    $$
    \xi = z/L \tag{7.1}
    $$

-   **Liquid tritium concentration $x_{T}$** scaled relative to the inlet concentration $c_{T}(L^{+})$, equation (9) in the paper:
    $$
    x_{T} = c_{T}/c_{T}(L^{+}) \tag{7.2}
    $$

-   **Pressure ratio $\psi$** is the ratio of the hydrostatic pressure head to the gas inlet pressure $P_{0}$:
    $$
    \psi = \frac{\rho_{l}g(1-\epsilon_{g})L}{P_{0}} \tag{7.3}
    $$
    This is equation (15) in the paper

    By substituting `z = ξL` (from eq. 7.1) into the dimensional pressure equation (5.1) and then factoring out `P_0`:
    $$
    P = P_{0} - \rho_{l}g(1-\epsilon_{g})\xi L = P_{0} \left( 1 - \frac{\rho_{l}g(1-\epsilon_{g})L}{P_0}\xi \right)
    $$
    The fractional term is the definition of `ψ` (eq. 7.3), which simplifies the expression to:

    $$
    P = P_{0}(1-\xi\psi) \tag{7.4}
    $$


-   **Equilibrium ratio $\nu$** is the ratio of tritium partial pressure that would be in equilibrium with the inlet liquid to the total gas pressure at the inlet, which is equation (15) in the paper:
    $$
    v = \frac{[c_{T}(L^{+})/K_{s}]^{2}}{P_{0}} \tag{7.5}
    $$

-   **Dimensionless driving force $\theta$**: The liquid-gas tritium flux $J_{T}$ is driven by the term $(c_{T} - c_{T,eq})$. It is made dimensionless by dividing by the liquid concentration at the inlet, $c_{T}(L^{+})$:
    $$
    \frac{(c_{T}-c_{T,eq})}{c_{T}(L^{+})} = \frac{c_{T}}{c_{T}(L^{+})} - \frac{c_{T,eq}}{c_{T}(L^{+})} = x_{T} - \frac{K_{s}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})} \tag{7.6}
    $$
    Bringing all terms inside the square root:
    $$
    \frac{K_{S}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})} \quad \rightarrow \quad \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5} \tag{7.7}
    $$
    Using the definition for $v$ gives $\frac{K_{s}^{2}}{(c_{T}(L^{+}))^{2}} = \frac{1}{vP_{0}}$, and with $P = P_{0}(1-\xi\psi)$:
    $$
    \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5} = \left[\left(\frac{1}{vP_{0}}\right) \cdot P_{0}(1-\xi\psi) \cdot y_{T_{2}}\right]^{0.5} = \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5}
    $$
    Finally, this gives the dimensionless driving force $\theta$ as:
    $$
    \theta = x_{T} - \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5} \tag{7.8}
    $$
    This is equation (12) in the paper.

-   **Bodenstein numbers**: $Bo_{l}$ and $Bo_{g}$ represent the ratio of advective (bulk flow) transport to dispersive (axial mixing) transport.
    $$
    Bo_{l} = \frac{u_{l}L}{(1-\epsilon_{g})E_{l}} \tag{7.9}
    $$
    Equation (13) in the paper.
    $$
    Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}} \tag{7.10}
    $$
    Equation (14) in the paper.
    (Note: The original paper has $u_g$, but it should be $u_{g0}$ to be consistent with the dimensionless group definition).

-   **Number of Transfer Units**: $\phi_{l}$ and $\phi_{g}$ compare the residence time of the fluid in the column to the characteristic time required for mass transfer.
    $$
    \phi_{l} = \frac{ah_{l}L}{u_{l}} \tag{7.11}
    $$
    Equation (13) in the paper.
    $$
    \phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \cdot \frac{ah_{l}L}{u_{g0}} \tag{7.12}
    $$
    Equation (14) in the paper.

## 8. Dimensionless ODE System

To solve the model, the original partial differential equations are transformed into a system of dimensionless ordinary differential equations (ODEs) using the variables and parameters defined in Section 8. This process simplifies the equations and groups physical parameters into meaningful dimensionless numbers.

### 8.1. Liquid Phase Dimensionless ODE

We start with the liquid phase mass balance (Equation 1.1):
$$
\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}} + u_{l}\frac{dc_{T}}{dz} - J_{T}a = 0 \quad [mol \cdot m^{-3} \cdot s^{-1}]
$$
First, we substitute the dimensionless variables for concentration $c_T = x_T c_T(L^{+})$ and axial position $z = \xi L$:
$$
\frac{d}{dz} = \frac{1}{L}\frac{d}{d\xi} \quad \text{and} \quad \frac{d^2}{dz^2} = \frac{1}{L^2}\frac{d^2}{d\xi^2}
$$
Equation 1.1 becomes:
$$
\frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}(x_T c_T(L^{+}))}{d\xi^{2}} + \frac{u_{l}}{L}\frac{d(x_T c_T(L^{+}))}{d\xi} - J_{T}a = 0 \tag{8.1.1}
$$
Since $c_T(L^{+})$ is a constant, it can be factored out. We also substitute the expression for the flux, $J_T = h_l(c_T - c_{T,eq})$:
$$
c_T(L^{+})\left( \frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}x_T}{d\xi^{2}} + \frac{u_{l}}{L}\frac{dx_T}{d\xi} \right) - ah_l(c_T - c_{T,eq}) = 0 \tag{8.1.2}
$$
Now, we divide the entire equation by $\frac{u_l c_T(L^{+})}{L}$:
$$
\frac{\epsilon_{l}E_{l}}{u_l L}\frac{d^{2}x_T}{d\xi^{2}} + \frac{dx_T}{d\xi} - \frac{ah_l L}{u_l} \frac{(c_T - c_{T,eq})}{c_T(L^{+})} = 0 \tag{8.1.3}
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the liquid phase: $Bo_{l} = \frac{u_{l}L}{\epsilon_{l}E_{l}}$ (assuming $\epsilon_l = 1-\epsilon_g$)
-   Number of transfer units for the liquid phase: $\phi_{l} = \frac{ah_{l}L}{u_{l}}$
-   Dimensionless driving force: $\theta = \frac{c_T - c_{T,eq}}{c_T(L^{+})}$

We arrive at the final dimensionless ODE for the liquid phase:
$$
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0 \tag{8.1.4}
$$
This is equation (10) in the paper.

### 8.2. Gas Phase Dimensionless ODE

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
u_g P y_{T_2} = \left(\frac{u_{g0}}{1-\xi\psi}\right) \left(P_0(1-\xi\psi)\right) y_{T_2} = u_{g0}P_0 y_{T_2} \tag{8.2.1}
$$
The derivative of the advection term with respect to $z$ is:
$$
\frac{d(u_{g}Py_{T_{2}})}{dz} = \frac{1}{L}\frac{d(u_{g0}P_0 y_{T_2})}{d\xi} = \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi} \tag{8.2.2}
$$

This simplification is the primary reason for using the molar fraction $y_{T_{2}}$ as the state variable for the gas phase instead of the partial pressure $P_{T_{2}}$. By choosing $y_{T_{2}}$, the position-dependent terms for gas velocity $u_g(z)$ and total pressure $P(z)$ cancel out in the advection term, leading to a much simpler derivative.

For the dispersion term, we must differentiate the full product $P(\xi)y_{T_2}(\xi)$ twice with respect to $\xi$. Let's denote $y_{T_2}$ as $y$ for brevity.
The first derivative is:
$$
\frac{d(Py)}{d\xi} = \frac{d}{d\xi}[P_0(1-\xi\psi)y] = P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right] \tag{8.2.3}
$$
The second derivative is:
$$
\frac{d^2(Py)}{d\xi^2} = \frac{d}{d\xi} \left( P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right] \right) = P_0 \left[ (1-\xi\psi)\frac{d^2y}{d\xi^2} - 2\psi \frac{dy}{d\xi} \right] \tag{8.2.4}
$$
The full dispersion term in equation 2.1 is therefore:
$$
-\epsilon_{g}E_{g}\frac{d^{2}(Py_{T_{2}})}{dz^{2}} = -\frac{\epsilon_{g}E_{g}}{L^2} \frac{d^2(Py_{T_2})}{d\xi^2} = -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] \tag{8.2.5}
$$
Substituting the dimensionless dispersion and advection terms back into the gas phase balance:
$$
-\frac{1}{RT} \left( -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] + \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi} \right) - \frac{1}{2}ah_l c_T(L^{+}) \theta = 0 \tag{8.2.6}
$$
Now, we divide the entire equation by $\frac{u_{g0}P_0}{RTL}$:
$$
\frac{\epsilon_{g}E_{g}}{u_{g0}L} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] - \frac{dy_{T_2}}{d\xi} + \frac{RTL}{2u_{g0}P_0} ah_l c_T(L^{+}) \theta = 0 \tag{8.2.7}
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the gas phase: $Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}}$
-   Number of transfer units for the gas phase: $\phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \frac{ah_{l}L}{u_{g0}}$

Substituting these groups into equation 9.2.7 gives:
$$
\frac{1}{Bo_g} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] - \frac{dy_{T_2}}{d\xi} + \phi_g \theta = 0 \tag{8.2.8}
$$
Next, we expand the first term by distributing the $\frac{1}{Bo_g}$ factor:
$$
\frac{(1-\xi\psi)}{Bo_g}\frac{d^2y_{T_2}}{d\xi^2} - \frac{2\psi}{Bo_g}\frac{dy_{T_2}}{d\xi} - \frac{dy_{T_2}}{d\xi} + \phi_g \theta = 0 \tag{8.2.9}
$$
Finally, we group the terms containing the first derivative $\frac{dy_{T_2}}{d\xi}$:

We arrive at the final dimensionless ODE for the gas phase:
$$
\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} - \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} + \phi_{g}\theta = 0 \tag{8.2.10}
$$
In the paper they write this equation with opposite signs:
$$
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0 \tag{8.2.11}
$$
This is equation (11) in the paper.
### 8.3. Summary of the System

The complete model is a system of two second-order ODEs for the dimensionless concentrations $x_T$ and $y_{T_2}$ as a function of the dimensionless position $\xi$.

The two governing equations are the liquid phase balance (8.1.4) and the gas phase balance (8.2.10):
$$
\begin{cases}
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0 \quad \text{(Liquid Phase)} \\
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0 \quad \text{(Gas Phase)}
\end{cases} \tag{8.3.1}
$$
Where:
$$
\theta = x_{T} - \sqrt{\frac{(1-\xi\psi)y_{T_{2}}}{\nu}} \tag{8.3.2}
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
\end{cases} \tag{8.3.3}
$$
This system of four first-order ODEs, along with the four boundary conditions from Section 10, can be solved using a numerical BVP solver, such as `scipy.integrate.solve_bvp` in Python.

## 9. Dimensionless Boundary Conditions

### 9.1 "Closed-Closed" BCs case
-   **Liquid outlet ($\xi=0$):** Zero-gradient condition states that tritium leaves the bottom of the column only through the bulk flow (advection). This is equation (16) in the paper.
    $$
    \frac{dx_{T}}{d\xi}\bigg|_{\xi=0} = 0 \tag{9.1.1}
    $$

-   **Liquid inlet ($\xi=1$):** The concentration of the liquid feed before it enters is different from the concentration just inside the column because of back-mixing. The size of this jump is determined by the Bodenstein number $Bo_{l}$. If there were no dispersion ($Bo_{l} \rightarrow \infty$), the concentration profile would be continuous. This is equation (16) in the paper.
    $$
    x_{T}(1) = 1 - \frac{1}{Bo_{l}}\frac{dx_{T}}{d\xi}\bigg|_{\xi=1} \tag{9.1.2}
    $$

-   **Gas inlet ($\xi=0$):** The concentration of the incoming gas is immediately altered to a slightly higher value just inside the column due to tritium from higher up mixing back down via dispersion. This is equation (17) in the paper.
    $$
    y_{T_{2}}(0) = y_{T_{2}}(0^{-}) + \frac{1}{Bo_{g}}\frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=0} \tag{9.1.3}
    $$

-   **Gas outlet ($\xi=1$):** Zero-gradient condition enforces that tritium leaves the top of the column only through the bulk upward flow of the gas bubbles (advection). This is equation (17) in the paper.
    $$
    \frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=1} = 0 \tag{9.1.4}
    $$

### 9.2 "Open-Closed" BCs case
-   **Liquid outlet ($\xi=0$):** 
    $$
    \frac{dx_{T}}{d\xi}\bigg|_{\xi=0} = 0 \tag{9.2.1}
    $$

-   **Liquid inlet ($\xi=1$):** 
    $$
    x_{T}(1) = 1  \tag{9.2.2}
    $$

-   **Gas inlet ($\xi=0$):** 
    $$
    y_{T_{2}}(0) = y_{T_{2}}(0^{-}) \tag{9.2.3}
    $$

-   **Gas outlet ($\xi=1$):**
    $$
    \frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=1} = 0 \tag{9.2.4}
    $$

## 10. Parameter Acquisition

The core model requires several parameters that describe the hydrodynamic and mass transfer behaviour of the bubble column (e.g., $\epsilon_g, a, h_l, E_l, E_g$). These parameters are not fundamental material properties but depend on the column's geometry, operating conditions, and the physical properties of the fluids. They are typically estimated using well-established empirical correlations derived from experimental data.

### 10.1. Bond, Galilei, Schmidt & Froude numbers ($Bn$, $Ga$, $Sc$, & $Fr$)

$$
Fr = \frac{u_{g0}}{\sqrt{gD}} \tag{10.1.1}
$$

$$
Bn = \frac{g D^{2} \rho_{l}}{\sigma_{l}} \tag{10.1.2}
$$
$$
Ga = \frac{g D^{3}}{\nu_{l}^{2}} \tag{10.1.3}
$$
$$
Sc = \frac{\nu_{l}}{D_{T,l}} \tag{10.1.4}
$$

### 10.2. Liquid Phase Dispersion Coefficient ($E_l$)

Correlation for the liquid-phase dispersion coefficient in terms of the column diameter and liquid properties:

$$
E_{l} = \frac{D \cdot u_{g0}}{((13Fr)/(1+6.5 \cdot (Fr)^{0.8}))} \tag{10.2.1}
$$

### 10.3. Gas Phase Dispersion Coefficient ($E_g$)

A commonly used correlation for the gas-phase dispersion coefficient:

$$
E_{g} = 0.2\,D^{2}u_{g0} \tag{10.3.1}
$$


### 10.4. Gas Holdup ($\epsilon_g$)

A correlation that has been used to estimate the gas holdup in similar systems:

$$
\frac{\epsilon_{g}}{(1 - \epsilon_{g})^{4}} = 0.2\,Bn^{1/8}\,Ga^{1/12}\,Fr \tag{10.4.1}
$$


### 10.5. Liquid Phase Fraction ($\epsilon_l$)

The liquid phase fraction is simply related to the gas phase fraction:

$$
\epsilon_{l} = 1 - \epsilon_{g} \tag{10.5.1}
$$

### 10.6. Mean Bubble Diameter ($d_b$)

$$
d_b = (26\,Bn^{-0.5}\,Ga^{-0.12}\,Fr^{-0.12}) \cdot D \tag{10.6.1}
$$

### 10.7. Specific Interfacial Area ($a$)

In this model we assume the bubbles are spherical, have a mean diameter of $d_b$, and are evenly distributed across the column.  Knowing the gas phase fraction $\epsilon_{g}$, we can then determine the specific interfacial area $a$:

Bubble volume $v_b$: 
$$\frac{4}{3} \cdot \pi \cdot r_b^3$$

Bubble area $a_b$:
$$a_b = 4 \cdot \pi \cdot r_b^2$$

Number of bubbles per unit volume of two-phase mixture, $n_b$:
$$n_b=\frac{\epsilon_{g}}{v_b}$$

Interfacial area per unit volume $a$:
$$a = n_b \cdot a_b = \frac{\epsilon_{g}}{\frac{4}{3} \cdot \pi \cdot r_b^3} \cdot (4 \cdot \pi \cdot r_b^2) = \frac{\epsilon_{g}}{\frac{1}{3} \cdot \frac{d_b}{2}}$$

$$a =  \frac{6\epsilon_{g}}{d_b} \quad [{m^{-1}}] \tag{10.7.1}$$

### 10.8. Volumetric Mass Transfer Coefficient ($a \cdot h_l$)

An empirical correlation for the volumetric mass transfer coefficient used in the reference paper:

$$
\frac{a h_{l} D^{2}}{D_{T,l}} = 0.6\,Sc^{0.5}\,Bn^{0.62}\,Ga^{0.31}\,\epsilon_{g}^{1.1}
$$
$$
a h_{l}= \frac{(0.6\,Sc^{0.5}\,Bn^{0.62}\,Ga^{0.31}\,\epsilon_{g}^{1.1}) \cdot D_{T,l}}{D^2} \tag{10.8.1}
$$

### 10.9. Mass Transfer Coefficient ($h_l$)

When the interfacial area $a$ is known, the mass transfer coefficient can be calculated from the volumetric mass transfer coefficient:

$$
h_l = \frac{a h_{l}}{a} \tag{10.9.1}
$$

### 11. Physical Properties Correlations

The Malara paper provides several correlations for the physical properties of LiPb, which are implemented in this model. 

-   **LiPb Density :** $\rho_{l} = 10.45 \cdot 10^{3} \cdot (1-1.64 \cdot 10^{-4} \cdot T) \quad [{kg/m^3}] $
-   **LiPb Dynamic Viscosity :** $ \mu_{l} = 1.87 \cdot 10^{-4} \cdot e^{\frac{11640}{RT}} \quad [{Pa \cdot s}] $
-   **LiPb Tritium Diffusivity :** $D_{T,l} = 2.7 \cdot 10^{-7} \cdot e^{\frac{-27000}{RT}} \quad [{m^2/s}] $
-   **Solubility of Tritium in LiPb, Sieverts constant :** $K_{l, atfrac} = 2.32 \cdot 10^{-8} \cdot e^{\frac{-1350}{RT}} \quad [at.frac \cdot Pa^{0.5}]$

To get the Sieverts constant in SI units:
-   **Solubility of Tritium in LiPb, Sieverts constant :** $K_{l} = K_{l, atfrac} \cdot \frac{\rho_{l}}{M \cdot N_{Av}} \quad [mol \cdot m^{3} \cdot Pa^{0.5}]$

Where $M$ is the molecular weight of LiPb $(M=2.875 \cdot 10^{-25}) \quad [{kg/mol}]$ and $N_{Av}$ is the Avogadro constant.

To obtain the kinematic viscosity $\nu_{l}$ we divide the dynamic viscosity by the density:
-   **LiPb Kinematic Viscosity :** $\nu_{l} = \frac{\mu_{l}}{\rho_{l}} \quad [{m^2/s}] $

The Malara paper does not provide a correlation for surface tension of the LiPb - Helium interface, however this property is needed for computing the Bond number. We use the correlation from [C] in this model:

-   **LiPb - Helium Surface Tension :** $\sigma_{l} = 0.52 - 0.11 \cdot 10^{-3} \cdot T \quad [{N/m}]$


 



### References
[A] - TRITIUM EXTRACTION FROM LIQUID Pb-l6li: A CRJTICAL REVIEW OF CANDIDATE TECHNOLOGIES
FOR ITER AND DEMO APPLICATIONS

[B] - Experimental design of tritium extraction loop from lead lithium eutectic

[C] - Lead–lithium eutectic material database for nuclear fusion technology

