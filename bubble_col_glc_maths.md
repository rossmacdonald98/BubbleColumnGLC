# Mathematical Model of a Gas-Liquid Contactor Bubble Column

This document explains the mathematical equations used to model the behavior of tritium in a gas-liquid contactor bubble column, based on the paper "Tritium Extraction from Pb-17Li by Bubble Columns" by Carlo Malara. The model consists of mass balance equations for both the liquid and gas phases, boundary conditions, and several supporting equations.

The original paper is available at: https://doi.org/10.13182/FST95-A30485

## 1. Liquid phase Tritium mass balance

$$
\underbrace{\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}}}_{\text{axial dispersion}} + \underbrace{u_{l}\frac{dc_{T}}{dz}}_{\text{convection}} - \underbrace{J_{T}a}_{\text{mass transfer}} = 0
$$
This is equation (1) in the paper.

**Axial dispersion** acts to flatten the tritium concentration profile due to the random net movement of particles from areas of high concentration to low concentration. If the local concentration profile is a valley, the curvature and dispersion terms are positive; if it's a peak, dispersion is negative.

**Convection** either locally accumulates or removes tritium as a result of the fluid velocity and concentration gradient. If the upstream concentration is greater than downstream (positive gradient), there is more tritium flowing in than out and tritium accumulates, and vice-versa. If the concentration profile is flat then the flow in equals the flow out and the net accumulation/removal is 0.

**Mass transfer** is the flux of tritium atoms from the liquid to the gas.


## 2. Gas phase Tritium mass balance

$$
\frac{1}{RT} \left( \underbrace{-\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}}}_{\text{axial dispersion}} + \underbrace{\frac{d(u_{g}Py_{T_{2}})}{dz}}_{\text{advection}} \right) - \underbrace{\frac{1}{2}J_{T}a}_{\text{mass transfer}} = 0
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
    \frac{dc_{T}}{dz}\bigg|_{z=0} = 0
    $$

-   **At liquid inlet ($z=L$):** Tritium entering the column via convection is equal to the sum entering/leaving via dispersion.
    $$
    \underbrace{u_{l}[c_{T}(L) - c_{T}(L^{+})]}_{\text{convective flux}} = \underbrace{-\epsilon_{l}E_{l}\frac{dc_{T}}{dz}\bigg|_{z=L}}_{\text{dispersive flux}}
    $$


## 4. Gas Phase Boundary Conditions
This is equation (4) in the paper.
-   **At gas inlet ($z=0$):** Tritium entering the column via advection is equal to the sum entering/leaving via dispersion.
    $$
    \underbrace{u_{g}[y_{T_{2}}(0^{-}) - y_{T_{2}}(0)]}_{\text{advective flux}} = \underbrace{-\epsilon_{g}E_{g}\frac{dy_{T_{2}}}{dz}\bigg|_{z=0}}_{\text{dispersive flux}}
    $$

-   **At gas outlet ($z=L$):** Tritium concentration profile is flat. Only advection, no dispersion.
    $$
    \frac{dy_{T_{2}}}{dz}\bigg|_{z=L} = 0
    $$


## 5. Gas Phase Mass Balance
This is equation (5) in the paper.

Enforces that the only change in the total gas mass flow along the length of the column is due to the source term of tritium flux from the liquid to the gas.

$$
\underbrace{\frac{1}{RT}\frac{d(Pu_{g})}{dz}}_{\text{Change in total gas flow}} - \underbrace{\frac{1}{2}J_{T}a}_{\text{Source term}} = 0
$$


## 6. Gas Phase Velocity

The pressure in the column is assumed to be purely hydrostatic, based on the height $z$ and pressure $P_{0}$ at $z=0$.

$$
P = P_{0} - \rho_{l}g(1-\epsilon_{g})z
$$
This is equation (6) in the paper.

By neglecting the change in mass flow due to the source term (assuming it is negligible), equation 5 becomes much simpler:

$$
\frac{d(Pu_{g})}{dz} = 0 \rightarrow Pu_{g} = \text{constant}
$$

$$
P(z) \cdot u_{g}(z) = P_{0} \cdot u_{g0} \rightarrow u_{g}(z) = \frac{P_{0} \cdot u_{g0}}{P(z)}
$$

Substituting in the pressure profile gives:

$$
u_{g}(z) = u_{g0}\frac{P_{0}}{P_{0} - \rho_{l}g(1-\epsilon_{g})z}
$$

This is equation (7) in the paper.


## 7. Liquid – Gas Flux

A simple linear driving force model is used to specify the rate of tritium transport across the liquid-gas interface. This is some mass transfer coefficient $h_{l}$ multiplied by the difference between the liquid concentration and concentration when it is in equilibrium with the gas phase partial pressure.

From Sievert's law:

$$
c_{T,eq} = K_{S}(P_{T_{2}})^{0.5}
$$

The flux $J_{T}$ is:

$$
J_{T}(z) = h_{l}(c_{T} - K_{s}(Py_{T_{2}})^{0.5})
$$

This is equation (8) in the paper.


## 8. Dimensionless Equations

Nondimensionalization of the equations makes the model easier to solve. The following nondimensional variables are used:

-   **Axial position $\xi$** from $z=0$ to $z=L$, equation (9) in the paper:
    $$
    \xi = z/L
    $$

-   **Liquid tritium concentration $x_{T}$** scaled relative to the inlet concentration $c_{T}(L^{+})$, equation (9) in the paper:
    $$
    x_{T} = c_{T}/c_{T}(L^{+})
    $$

-   **Pressure ratio $\psi$** is the ratio of the hydrostatic pressure head to the gas inlet pressure $P_{0}$:
    $$
    \psi = \frac{\rho_{l}g(1-\epsilon_{g})L}{P_{0}}
    $$
    From this, the pressure profile becomes:
    $$
    P = P_{0}(1-\xi\psi)
    $$
    This is equation (15) in the paper

-   **Equilibrium ratio $\nu$** is the ratio of tritium partial pressure that would be in equilibrium with the inlet liquid to the total gas pressure at the inlet, equation (15) in the paper:
    $$
    v = \frac{[c_{T}(L^{+})/K_{s}]^{2}}{P_{0}}
    $$

-   **Dimensionless driving force $\theta$**: The liquid-gas tritium flux $J_{T}$ is driven by the term $(c_{T} - c_{T,eq})$. It is made dimensionless by dividing by the liquid concentration at the inlet, $c_{T}(L^{+})$:
    $$
    \frac{(c_{T}-c_{T,eq})}{c_{T}(L^{+})} = \frac{c_{T}}{c_{T}(L^{+})} - \frac{c_{T,eq}}{c_{T}(L^{+})} = x_{T} - \frac{K_{s}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})}
    $$
    Bringing all terms inside the square root:
    $$
    \frac{K_{S}(Py_{T_{2}})^{0.5}}{c_{T}(L^{+})} = \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5}
    $$
    Using the definition for $v$ gives $\frac{K_{s}^{2}}{(c_{T}(L^{+}))^{2}} = \frac{1}{vP_{0}}$, and with $P = P_{0}(1-\xi\psi)$:
    $$
    \left[\frac{K_{S}^{2} \cdot P \cdot y_{T_{2}}}{(c_{T}(L^{+}))^{2}}\right]^{0.5} = \left[\left(\frac{1}{vP_{0}}\right) \cdot P_{0}(1-\xi\psi) \cdot y_{T_{2}}\right]^{0.5} = \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5}
    $$
    Finally, this gives the dimensionless driving force $\theta$ as:
    $$
    \theta = x_{T} - \left(\frac{(1-\xi\psi) \cdot y_{T_{2}}}{v}\right)^{0.5}
    $$
    This is equation (12) in the paper.

-   **Bodenstein numbers**: $Bo_{l}$ and $Bo_{g}$ represent the ratio of advective (bulk flow) transport to dispersive (axial mixing) transport.
    $$
    Bo_{l} = \frac{u_{l}L}{(1-\epsilon_{g})E_{l}}
    $$
    Equation (13) in the paper.
    $$
    Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}}
    $$
    Equation (14) in the paper.
    (Note: The original paper has $u_g$, but it should be $u_{g0}$ to be consistent with the dimensionless group definition).

-   **Number of Transfer Units**: $\phi_{l}$ and $\phi_{g}$ compare the residence time of the fluid in the column to the characteristic time required for mass transfer.
    $$
    \phi_{l} = \frac{ah_{l}L}{u_{l}}
    $$
    Equation (13) in the paper.
    $$
    \phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \cdot \frac{ah_{l}L}{u_{g0}}
    $$
    Equation (14) in the paper.

## 9. Dimensionless ODE System

To solve the model, the original partial differential equations are transformed into a system of dimensionless ordinary differential equations (ODEs) using the variables and parameters defined in Section 8. This process simplifies the equations and groups physical parameters into meaningful dimensionless numbers.

### 9.1. Liquid Phase ODE

We start with the liquid phase mass balance (Equation 1):
$$
\epsilon_{l}E_{l}\frac{d^{2}c_{T}}{dz^{2}} + u_{l}\frac{dc_{T}}{dz} - J_{T}a = 0
$$
First, we substitute the dimensionless variables for concentration $c_T = x_T c_T(L^{+})$ and axial position $z = \xi L$:
$$
\frac{d}{dz} = \frac{1}{L}\frac{d}{d\xi} \quad \text{and} \quad \frac{d^2}{dz^2} = \frac{1}{L^2}\frac{d^2}{d\xi^2}
$$
The equation becomes:
$$
\frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}(x_T c_T(L^{+}))}{d\xi^{2}} + \frac{u_{l}}{L}\frac{d(x_T c_T(L^{+}))}{d\xi} - J_{T}a = 0
$$
Since $c_T(L^{+})$ is a constant, it can be factored out. We also substitute the expression for the flux, $J_T = h_l(c_T - c_{T,eq})$:
$$
c_T(L^{+})\left( \frac{\epsilon_{l}E_{l}}{L^2}\frac{d^{2}x_T}{d\xi^{2}} + \frac{u_{l}}{L}\frac{dx_T}{d\xi} \right) - ah_l(c_T - c_{T,eq}) = 0
$$
Now, we divide the entire equation by $\frac{u_l c_T(L^{+})}{L}$:
$$
\frac{\epsilon_{l}E_{l}}{u_l L}\frac{d^{2}x_T}{d\xi^{2}} + \frac{dx_T}{d\xi} - \frac{ah_l L}{u_l} \frac{(c_T - c_{T,eq})}{c_T(L^{+})} = 0
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the liquid phase: $Bo_{l} = \frac{u_{l}L}{\epsilon_{l}E_{l}}$ (assuming $\epsilon_l = 1-\epsilon_g$)
-   Number of transfer units for the liquid phase: $\phi_{l} = \frac{ah_{l}L}{u_{l}}$
-   Dimensionless driving force: $\theta = \frac{c_T - c_{T,eq}}{c_T(L^{+})}$

We arrive at the final dimensionless ODE for the liquid phase:
$$
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0
$$
This is equation (10) in the paper.

### 9.2. Gas Phase ODE

The derivation for the gas phase is more involved. We start with the gas phase mass balance (Equation 2):
$$
\frac{1}{RT} \left( -\epsilon_{g}E_{g}\frac{d^{2}Py_{T_{2}}}{dz^{2}} + \frac{d(u_{g}Py_{T_{2}})}{dz} \right) - \frac{1}{2}J_{T}a = 0
$$
We substitute the dimensionless variables and expressions for $P(\xi)$ and $u_g(\xi)$:
-   $P(\xi) = P_0(1-\xi\psi)$
-   $u_g(\xi) = \frac{u_{g0}}{1-\xi\psi}$
-   $J_T = h_l c_T(L^{+}) \theta$

The advection term $u_g P y_{T_2}$ simplifies nicely:
$$
u_g P y_{T_2} = \left(\frac{u_{g0}}{1-\xi\psi}\right) \left(P_0(1-\xi\psi)\right) y_{T_2} = u_{g0}P_0 y_{T_2}
$$
The derivative of the advection term with respect to $z$ is:
$$
\frac{d(u_{g}Py_{T_{2}})}{dz} = \frac{1}{L}\frac{d(u_{g0}P_0 y_{T_2})}{d\xi} = \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi}
$$
For the dispersion term, we must differentiate the full product $P(\xi)y_{T_2}(\xi)$ twice with respect to $\xi$. Let's denote $y_{T_2}$ as $y$ for brevity.
The first derivative is:
$$
\frac{d(Py)}{d\xi} = \frac{d}{d\xi}[P_0(1-\xi\psi)y] = P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right]
$$
The second derivative is:
$$
\frac{d^2(Py)}{d\xi^2} = \frac{d}{d\xi} \left( P_0 \left[ -\psi y + (1-\xi\psi)\frac{dy}{d\xi} \right] \right) = P_0 \left[ (1-\xi\psi)\frac{d^2y}{d\xi^2} - 2\psi \frac{dy}{d\xi} \right]
$$
The full dispersion term in the original equation is therefore:
$$
-\epsilon_{g}E_{g}\frac{d^{2}(Py_{T_{2}})}{dz^{2}} = -\frac{\epsilon_{g}E_{g}}{L^2} \frac{d^2(Py_{T_2})}{d\xi^2} = -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right]
$$
Substituting the full dispersion and advection terms back into the gas phase balance:
$$
\frac{1}{RT} \left( -\frac{\epsilon_{g}E_{g}P_0}{L^2} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] + \frac{u_{g0}P_0}{L}\frac{dy_{T_2}}{d\xi} \right) - \frac{1}{2}ah_l c_T(L^{+}) \theta = 0
$$
Now, we divide the entire equation by $\frac{u_{g0}P_0}{RTL}$:
$$
-\frac{\epsilon_{g}E_{g}}{u_{g0}L} \left[ (1-\xi\psi)\frac{d^2y_{T_2}}{d\xi^2} - 2\psi \frac{dy_{T_2}}{d\xi} \right] + \frac{dy_{T_2}}{d\xi} - \frac{RTL}{2u_{g0}P_0} ah_l c_T(L^{+}) \theta = 0
$$
Recognizing the dimensionless groups:
-   Bodenstein number for the gas phase: $Bo_{g} = \frac{u_{g0}L}{\epsilon_{g}E_{g}}$
-   Number of transfer units for the gas phase: $\phi_{g} = \frac{1}{2}\frac{RTc_{T}(L^{+})}{P_{0}} \frac{ah_{l}L}{u_{g0}}$

We arrive at the final dimensionless ODE for the gas phase:
$$
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0
$$
This is equation (11) in the paper.
### 9.3. Summary of the System

The complete model is a system of two second-order ODEs for the dimensionless concentrations $x_T$ and $y_{T_2}$ as a function of the dimensionless position $\xi$.

$$
\begin{cases}
\frac{1}{Bo_{l}}\frac{d^{2}x_{T}}{d\xi^{2}} + \frac{dx_{T}}{d\xi} - \phi_{l}\theta = 0 \\
-\frac{(1-\xi\psi)}{Bo_{g}}\frac{d^{2}y_{T_{2}}}{d\xi^{2}} + \left(1 + \frac{2\psi}{Bo_{g}}\right)\frac{dy_{T_{2}}}{d\xi} - \phi_{g}\theta = 0
\end{cases}
$$
where the coupling term $\theta$ is given by:
$$
\theta = x_{T} - \sqrt{\frac{(1-\xi\psi)y_{T_{2}}}{\nu}}
$$
This system is solved over the domain $\xi \in [0, 1]$ subject to the four dimensionless boundary conditions specified in Section 10.

## 10. Dimensionless Boundary Conditions

-   **Liquid outlet ($\xi=0$):** Zero-gradient condition states that tritium leaves the bottom of the column only through the bulk flow (advection). This is equation (16) in the paper.
    $$
    \frac{dx_{T}}{d\xi}\bigg|_{\xi=0} = 0
    $$

-   **Liquid inlet ($\xi=1$):** The concentration of the liquid feed before it enters is different from the concentration just inside the column because of back-mixing. The size of this jump is determined by the Bodenstein number $Bo_{l}$. If there were no dispersion ($Bo_{l} \rightarrow \infty$), the concentration profile would be continuous. This is equation (16) in the paper.
    $$
    x_{T}(1) = 1 - \frac{1}{Bo_{l}}\frac{dx_{T}}{d\xi}\bigg|_{\xi=1}
    $$

-   **Gas inlet ($\xi=0$):** The concentration of the incoming gas is immediately altered to a slightly higher value just inside the column due to tritium from higher up mixing back down via dispersion. This is equation (17) in the paper.
    $$
    y_{T_{2}}(0) = y_{T_{2}}(0^{-}) + \frac{1}{Bo_{g}}\frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=0}
    $$

-   **Gas outlet ($\xi=1$):** States that tritium leaves the top of the column only through the bulk upward flow of the gas bubbles. This is equation (17) in the paper.
    $$
    \frac{dy_{T_{2}}}{d\xi}\bigg|_{\xi=1} = 0
    $$

---



## Nomenclature

* **$a$**: liquid-gas interfacial area per unit of volume [$m^{-1}$]
* **$Bo_g$**: Bodenstein number for the gas phase (eq. (14)) [-]
* **$Bo_l$**: Bodenstein number for the liquid phase (eq. (13)) [-]
* **$c_T$**: tritium concentration in liquid phase [$mol \cdot m^{-3}$]
* **$c_{T,e}$**: tritium concentration in liquid phase at equilibrium with the gas phase [$mol \cdot m^{-3}$]
* **$D_{T,l}$**: tritium diffusion coefficient in liquid phase[$m^2 \cdot s^{-1}$]
* **$E_l$**: dispersion coefficient in the liquid phase [$m^2 \cdot s^{-1}$]
* **$E_g$**: dispersion coefficient in the gas phase [$m^2 \cdot s^{-1}$]
* **$g$**: gravitational acceleration [$m \cdot s^{-2}$]
* **$h_l$**: tritium mass transfer coefficient from the liquid to the gas phase [$m \cdot s^{-1}$]
* **$J_T$**: tritium flux from liquid to gas phase [$mol \cdot m^{-2} \cdot s^{-1}$]
* **$K_s$**: tritium Sievert's constant in liquid [$mol \cdot m^{-3} \cdot Pa^{1/2}$]
* **$L$**: extractor height [m]
* **$P$**: total gas pressure [Pa]
* **$P_0$**: total gas pressure at the extractor inlet [Pa]
* **$R$**: gas constant [$J \cdot mol^{-1} \cdot K^{-1}$]
* **$T$**: temperature [K]
* **$u_g$**: gas velocity [$m \cdot s^{-1}$]
* **$u_{g0}$**: gas velocity at the extractor inlet [$m \cdot s^{-1}$]
* **$u_l$**: liquid velocity [$m \cdot s^{-1}$]
* **$x_T$**: dimensionless variable defined by eq. (9) [-]
* **$y_{T_2}$**: tritium molar fraction in the gas phase [-]
* **$z$**: axial coordinate in the extractor [m]

## Greek Letters

* **$\epsilon_g$**: gas phase fraction [-]
* **$\epsilon_l$**: liquid phase fraction [-]
* **$\xi$**: dimensionless axial position defined by eq. (9) [-]
* **$\theta$**: dimensionless driving force ratio defined by eq. (12) [-]
* **$\phi_l$**: dimensionless 'number of transfer units' parameter for the liquid phase defined by eq. (13) [-]
* **$\phi_g$**: dimensionless 'number of transfer units' parameter for the gas phase parameter defined by eq. (14) [-]
* **$\psi$**: dimensionless pressure ratio defined by eq. (15) [-]
* **$\nu$**: dimensionless equilibrium ratio defined by eq. (15) [-]
* **$\rho_l$**: Liquid density [$kg \cdot m^{-3}$]
