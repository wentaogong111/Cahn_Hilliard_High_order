# 2D Cahn-Hilliard Equation with Pseudospectral Method and 4th order time Integration
## Introduction

We aim to solve the 2D Cahn-Hilliard equation:

$$
\frac{\partial c}{\partial t} = \nabla^2(-\nabla^2 c + f(c))
$$

The variable $c(x,y,t)$ represents the concentration.

## Spatial Discretization Using Pseudospectral Method

### Fourier Transform

The Fourier coefficients are given by 
$$
\widehat{c}_{\boldsymbol{k}}(t)=\mathcal {FT}{c(\boldsymbol{r},t)}=\int c(\boldsymbol{r},t)e^{i\boldsymbol{k}\boldsymbol{r}}d\boldsymbol{r}
$$
and $k_i= \{-\pi N_i/L_i, -\pi(N_i-1)/L_i, \ldots, \pi(N_i-1)/L_i,\pi N_i/L_i\}$, where $N_i = L_i/\Delta_i$  and $\Delta_i$ is the stepsize of the meshgrid on the i direction.

The Fourier transform of the dynamical equation is 
$$
\frac{\partial\widehat{c}_{\boldsymbol{k}}}{\partial t}= M[-\kappa k^4\widehat{c}_{\boldsymbol{k}}-k^2\mathcal {FT}\{\frac{\delta f}{\delta c}\}]
$$

Explanation of Terms

1. **$\frac{\partial \hat{c}(k, t)}{\partial t}$**: Rate of change of the Fourier component of $c$ with respect to time $t$.
  
2. **$-k^4 \hat{c}_{\boldsymbol{k}}$**: The Laplacian term in Fourier space. $k^4$ is a consequence of two spatial derivatives $\nabla^2$.

3. **$k^2 \hat{f'}$**: Fourier transform of the derivative of free energy density. $k^2$ arises from one spatial derivative $\nabla^2$.




### Implicit Euler Method (Existing Method)

The existing implicit Euler scheme for solving the equation is:

$$
\widehat{c}_{\boldsymbol{k}}^{n+1} = \frac{\widehat{c}^n_{\boldsymbol{k}}-\Delta tMk^2\mathcal{FT}\{f'(c^n)\} }{1+\Delta t\kappa M k^4}
$$
Solver for $\widehat{c}_{\boldsymbol{k}}^{n+1}$: isolate $\widehat{c}_{\boldsymbol{k}}^{n+1}$ to find its value at the next time step:
$$
\widehat{c}_{\boldsymbol{k}}^{n+1}=\frac{\widehat{c}^n_{\boldsymbol{k}}-\Delta tMk^2\mathcal{FT}\{f'(c^n)\} }{1+\Delta t\kappa M k^4}
$$
where $\kappa$ is the gradient coefficient, **k** is the wave vector, M is mobility, $\Delta t$ is the time step value

### Using RK4 for first four time steps (Explicit method)

You can adapt the RK4 scheme to solve the Cahn-Hilliard equation. For the RHS defined as:

$$
\text{RHS} =M[-\kappa k^4\widehat{c}_{\boldsymbol{k}}^n-k^2\mathcal {FT}\{\frac{\delta f}{\delta c}\}]
$$

The RK4 scheme becomes:

$$
k_1 = \Delta t * RHS(\hat{c}_{\boldsymbol{k}}^n)
$$
$$
k_2 = \Delta t * RHS(\hat{c}_{\boldsymbol{k}}^n + k_1 / 2)
$$
$$
k_3 = \Delta t * RHS(\hat{c}_{\boldsymbol{k}}^n + k_2 / 2)
$$
$$
k_4 = \Delta t * RHS(\hat{c}_{\boldsymbol{k}}^n + k_3)
$$
$$
\hat{c}_{\boldsymbol{k}}^{n+1} = \hat{c}_{\boldsymbol{k}}^n + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
$$


### The fourth-order Adams-Moulton method for subsequent time steps (Implict method)

To adapt this for the Cahn-Hilliard equation in Fourier space, we can use:

$$
\text{RHS}(\hat{c}) = M[-\kappa k^4 \hat{c}_{\boldsymbol{k}}^n - k^2 \mathcal{FT}\{f'(c)\}]
$$

And the adapted Adams-Moulton formula becomes:

$$
\hat{c}_{\boldsymbol{k}}^{n+4} = \hat{c}_{\boldsymbol{k}}^{n+3} + \frac{\Delta t}{720} \left[ 251 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+4}) + 646 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+3}) - 264 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+2}) + 106 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+1}) - 19 \text{RHS}(\hat{c}_{\boldsymbol{k}}^n) \right]
$$


We want to isolate $\hat{c}_{\boldsymbol{k}}^{n+4}$ from the equation 


### Isolating $\hat{c}_{\boldsymbol{k}}^{n+4}$ from the Adams-Moulton Equation



To isolate $\hat{c}_{\boldsymbol{k}}^{n+4}$, we can rewrite the equation as follows:
$$
\hat{c}_{\boldsymbol{k}}^{n+4} \left( 1 - \frac{251 \Delta t M \kappa k^4}{720} \right) = \hat{c}_{\boldsymbol{k}}^{n+3} + \frac{\Delta t}{720} \left( 646 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+3}) - 264 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+2}) + 106 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+1}) - 19 \text{RHS}(\hat{c}_{\boldsymbol{k}}^n) \right)
$$

Finally, isolating $\hat{c}_{\boldsymbol{k}}^{n+4}$ gives:

$$
\hat{c}_{\boldsymbol{k}}^{n+4} = \frac{\hat{c}_{\boldsymbol{k}}^{n+3} + \frac{\Delta t}{720} \left( 646 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+3}) - 264 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+2}) + 106 \text{RHS}(\hat{c}_{\boldsymbol{k}}^{n+1}) - 19 \text{RHS}(\hat{c}_{\boldsymbol{k}}^n) \right)}{1 - \frac{251 \Delta t M \kappa k^4}{720}}
$$

1. **Inverse Fourier Transform**:Take the inverse Fourier transform of $\widehat{c}_{\boldsymbol{k}}^{n+4}$ to find ${c}_{\boldsymbol{k}}^{n+4}$ the concentration at the next time step in real space.
2. **Repeat**: Go back to RK4 and fourth-order Adams-Moulton method and repeat for all desired time steps.



## Hyperparameter Guidelines

### Mobility (M)
A smaller value of $M$ can make the time evolution more stable but slower. A larger value can make the system evolve very quickly and could lead to numerical instabilities.

### Interfacial Parameter ($\kappa$)
This controls the width of the interface between the two phases. A smaller value would mean a sharper interface, while a larger value results in a more diffuse interface. The value of $\kappa$ must be chosen carefully to ensure that the numerical grid can resolve the interface adequately.

### Well-Depth ($W$)
Controls the depth of the double well potential. Larger values make the free energy wells deeper, causing the phase separation to be more pronounced. However, this could also make the system more sensitive to numerical errors.

## Stability Analysis

Stability analysis for complex nonlinear PDEs like the Cahn-Hilliard equation can be challenging. Here are some empirical strategies for ensuring numerical stability:

### CFL Condition
Make sure the Courant-Friedrichs-Lewy (CFL) condition is satisfied. This means ensuring that

$$
dt \leq \text{CFL} \times \left( \frac{dx^2}{M} \right)
$$

where CFL is a number less than or equal to 1.

### Grid Resolution
Make sure the grid is fine enough to resolve the smallest features in your simulation. This is particularly important for the interface width, which is often related to 

$$
\sqrt{\frac{\kappa}{W}} 
$$.

### Initial Test
It may be beneficial to start with a simpler system (e.g., 1D instead of 2D, or a smaller grid) and scale up as you verify the model's stability.

### Sensitivity Analysis
Perform a sensitivity analysis by varying one parameter at a time while keeping the others constant. This will help you understand how each parameter influences the stability of the solution.

Example
![PNG](https://github.com/wentaogong111/Cahn_Hilliard_High_order/blob/main/code/cahn-hilliard-c0-0.5.png)
![GIF](https://github.com/wentaogong111/Cahn_Hilliard_High_order/blob/main/code/cahn-hilliard_high_order2.gif)

