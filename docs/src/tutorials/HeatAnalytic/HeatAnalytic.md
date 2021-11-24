# Heat Analytic
 In this example a heat transfer problem with known analytical solution is presented and solved. The FEM solution is verified using the analytical solution.

## problem description

 Let us consider the one-dimensional heat equation, $\partial_t T(x, t) = \alpha \partial^2_{xx}T(x, t)$ where $\alpha := k / \rho c$ is assumed uniform in the domain $[0,1]$ and constant. $Q_h = 0$ is also assumed. The boundary conditions are given by Dirichlet conditions at both boundaries for all times $T(0,t) = 0$ and $T(1,t)=0$. The initial condition is given by the following temperature distribution function:

```math
T(x, t=0) = \phi(x) := \sin \pi x + \frac{1}{2}\sin 3\pi x
```

```
using Plots, LaTeXStrings

ϕ(x) = sin(π*x) + sin(3π*x)/2

ne = 50 # number of elements
xdom = 0:1/ne:1
T0 = ϕ.(xdom)

plot(0:1e-3:1, ϕ, seriestype=:line, lab=L"\phi(x)",
     xlab=L"x", ylab=L"T", legend=:bottomright, title="Initial temperature profile")
plot!(xdom, T0, seriestype = :scatter, lab=L"T(x=x_e, 0)")
#savefig("plot_initial_temperature.svg")

nothing
```

The analytic solution in this case is

```math
T(x, t) = e^{- \pi^2 \alpha t}\sin \pi x + \frac{1}{2}e^{-(3\pi)^2 \alpha t}\sin 3\pi x,\qquad 0 \leq x \leq 1, t \geq 0.
```

## Numerical resolution

### 1D numerical solution
```
close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions
```

```
timeIncr  = 0.0001 ;   Tfinal    = 0.02   ;
rho       = 1.     ;   cSpHe     = 1.     ;
kCond     = 4      ;   L         = 1      ;
Area      = 0.25   ;   nelem     = 20     ;
Tdiri     = 0      ;   nPlots = 4 ;
```

```

diriDofs = [ 1 nelem+1 ];  robiDofs = [] ;

problemName     = 'Heat1DAnalytic' ;
initialTempFunc = 'myInitialTemp' ;
ambTempFuncName = 'myAmbTemp'     ;

ndivs     = [ nelem  ];
hConv = []; Tamb = [];
qInpLeft = []; qInpRight = [];
Tdiri = 0 ;

[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 L Area ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, ambTempFuncName, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], [], [], []);


```
### 3D numerical solution

```

timeIncr  = 0.0001 ; Tfinal    = 0.01   ;
Lx        = 1      ; Ly        = .5     ;
Lz        = .5     ; Tdiri     = 0      ;

assert( Area == (Ly*Lz) ) ;

diriDofs = []; robiDofs = [] ;

nPlots = inf ;

initialTempFunc = 'myInitialTemp' ;

ndivs     = [ 10 2 2 ];

hConv = []; Tamb = [];
qInpLeft = []; qInpRight = [];

Tdiri = 0 ; qInp = 0;

diriFaces = [ Tdiri 1 2 ];
neumFaces = [ qInp  3 4 5 6 ] ;
robiFaces = [  ] ;

[Ts3D, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], ...
  diriFaces, neumFaces, robiFaces  );






```
## solution verification

 analytical solution computation
```
% analytic solution computation
indplot = round(nelem/2) ;
xsAnly = NodesCoord(indplot) ;
alpha  = kCond / ( rho * cSpHe ) ;
tanali =  exp( -(  pi)^2 * alpha * times ) *       sin(     pi * xsAnly) ...
        + exp( -(3*pi)^2 * alpha * times ) * 0.5 * sin( 3 * pi * xsAnly ) ;
```
 plot settings
```
figure, grid on, hold on
MS = 18 ; LW = 1.5 ;

ntimes = length(times);
indsMarkers = floor( linspace(1,ntimes, 6)' ) ;

plot( times(indsMarkers(1)), Ts(indplot,indsMarkers(1)), 'b-x',  'markersize', MS,'linewidth',LW )
plot( times(indsMarkers(1)), tanali(indsMarkers(1)), 'r-o',  'markersize', MS,'linewidth',LW )

% numerical solution plot

plot( times, Ts(indplot,:), 'b',  'markersize', MS,'linewidth',LW )
plot( times(indsMarkers), Ts(indplot,indsMarkers), 'x',  'markersize', MS,'linewidth',LW )


% analytic solution plot
plot( times, tanali, 'r',  'markersize', MS,'linewidth',LW )
plot( times(indsMarkers), tanali(indsMarkers), 'o',  'markersize', MS,'linewidth',LW )
xlabel('t'), ylabel('Temp')

legend('FEM1D','analytic')

print( [ './pngs/' problemName '.png'],'-dpng' )
% ----------------------------------
