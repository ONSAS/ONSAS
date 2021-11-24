# Principle of Virtual Thermal Work

## Variational formulation of the heat flow equation


```math
    \int_{\Omega} k \nabla T \cdot \nabla \delta T dV
  + \int_{\Omega} \rho c \frac{\partial T}{\partial t} \delta T dV
  =
    \int_{\Omega} Q_h \delta T d V
  + \int_{\partial \Omega} k \nabla T \cdot \hat{\mathbf{n}} \delta T d S 
 \qquad \forall \delta T \in \tilde{\mathcal{T}}
```

Considering the boundary conditions

```math
\left\{
\begin{array}{lr}
T(\mathbf{x},t) = f_D(t) & \text{on} \,\Gamma_D \\
-k \frac{\partial T}{\partial n} (\mathbf{x},t) = f_N(\mathbf{x},t) & \text{on} \, \Gamma_N \\
-k \frac{\partial T}{\partial n} (\mathbf{x},t) = h \left( T(\mathbf{x},t)-T_\infty(t) \right)  & \text{on} \, \Gamma_R
\end{array}
\right.
```

where $h$ is the convection coefficient and $T_\infty(t)$ is the ambient temperature at time $t$.


```math
    \int_{\Omega} k \nabla T \cdot \nabla \delta T dV
  + \int_{\Omega} \rho c \frac{\partial T}{\partial t} \delta T dV
  + \int_{\Gamma_R} h T(\mathbf{x},t) \delta T d S 
  =
    \int_{\Omega} Q_h \delta T d V
  + \int_{\Gamma_N} q_{inp}(\mathbf{x},t)  \delta T d S 
  + \int_{\Gamma_R} h T_\infty(t) \delta T d S 

 \qquad \forall \delta T \in \tilde{\mathcal{T}}
```

where $q_{inp}$ is the input heat flow $q_{inp} = -f_N$.


## Finite Elements



```math
\mathbf{K}_{diff}^e = \frac{ k^e A^e}{\ell^e} 
\left[
\begin{matrix}
1 & -1 \\
-1 & 1
\end{matrix}
\right]
```

```math
\mathbf{C}_{intE}^e = \rho^e c^e A^e \ell^e \frac{1}{6} 
\left[
\begin{matrix}
2 & 1 \\
1 & 2
\end{matrix}
\right]
```

