# Uniaxial Extension Solid example
---

In this tutorial example an elastic solid is submitted to a uniaxial extension test. The problem is inspired by Exercise 4 from section 6.5 in (Holzapfel,2000). The geometry and tension applied are shown in the figure, where the $Lx$, $Ly$ and $Lz$ are the dimensions and the tension $p$ is applied on the face $x=Lx$, as nominal traction (see (Holzapfel,2000)).

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/diagramSolidUniaxialHTML.svg" alt="structure diagram" width="500"/>
```

```@raw latex
\begin{center}
\def\svgwidth{0.7\textwidth}
\input{diagramSolidUniaxialPDF.pdf_tex}
\end{center}
```

## Analytic solution
---
Let us consider that a uniform deformation is produced, with a nonzero axial stretch $\alpha$ and nonzero transversal stretch $\beta$. The corresponding deformation gradient and Green-Lagrange strain tensor are given by:
```math
\textbf{F} = \left[ \begin{matrix} \alpha & 0 & 0 \\ 0 & \beta & 0 \\ 0 & 0 & \beta \end{matrix} \right]
\qquad
\textbf{E} = \left[  \begin{matrix} \frac{1}{2} \left(\alpha^2 -1\right) & 0 & 0 \\ 0 &  \frac{1}{2} \left(\beta^2 -1\right) & 0 \\ 0 & 0 &  \frac{1}{2} \left(\beta^2 -1\right) \end{matrix}
\right]
```
The second Piola-Kirchhoff tensor $\textbf{S}$ is given by
```math
\textbf{S}( \textbf{E} ) = p_1 tr(\textbf{E}) \textbf{I} + 2 p_2 \textbf{E}
```
 then, using the relation $\textbf{P}=\textbf{F}\textbf{S}$, the $P_{yy}$ component is computed and set to zero (using the boundary conditions)
```math
P_{yy}( \textbf{E} ) =
p_1 \beta \left(
             \frac{1}{2} \left(\alpha^2 -1 \right) + \left( \beta^2 -1\right)
 \right) + 2 p_2 \beta (\frac{1}{2} \left(\beta^2 -1 \right)) = 0
```
then, using that $\beta\neq0$ (since $\text{det}( \textbf{F} ) \neq0$), we obtain
```math
 p_1 \frac{1}{2} \left(\alpha^2 -1 \right)
 = - (p_1+p_2) \left(\beta^2 -1 \right)
```
 then using the relation between the Lamé parameters $p_2$ and $p_1$ and the Young modulus and Poisson ratio, we obtain:

```math
 \left(\beta^2 -1 \right) = -\nu \left(\alpha^2 -1 \right).
```

 The axial component of the nominal stress is
```math
P_{xx}( \textbf{E} ) =
p_1 \alpha \left(
             \frac{1}{2} \left(\alpha^2 -1 \right) + \left( \beta^2 -1\right)
 \right) + 2 p_2 \alpha (\frac{1}{2} \left(\alpha^2 -1 \right)) = 0
```
 and substituting we obtain

```math
P_{xx}( \alpha ) =
p_1 \alpha \frac{1-2\nu}{2} \left(\alpha^2 -1 \right) + p_2 \alpha \left(\alpha^2 -1 \right) =
 \left( \frac{E \nu}{(1+\nu)2}  + \frac{E}{(1+\nu)2} \right)  \alpha \left(\alpha^2 -1 \right)
```
 thus, considering the axial displacement $u$ and using the stretch definition $\alpha = (1+u/Lx)$, we obtain
```math
P_{xx}( u ) =
 \frac{E}{2}  \left( \left( 1+\frac{u}{Lx} \right)^3 - \left( 1+ \frac{u}{Lx} \right) \right)
```
 where $u$ is the $x$ displacement of the points located on face $x=Lx$.

## Numerical solution: case 1
---
Before defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined.
```
clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ;
% scalar parameters
E = 1 ; nu = 0.3 ; p = 3 ; Lx = 2 ; Ly = 1 ; Lz = 1 ;
```


### MEBI parameters

#### materials
 The material of the solid considered is the Saint-Venant-Kirchhoff with Lamé parameters computed as
```
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
```
 since only one material is considered, a scalar struct is defined as follows
```
materials.hyperElasModel = 'SVK' ;
materials.hyperElasParams = [ lambda mu ] ;
```

#### elements
 In this model two kinds of elements are used: `tetrahedron` for the solid and `triangle` for introducing the external loads. Since two kinds of elements are used, the struct have length 2:
```
elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;
```

#### boundaryConds
 in this case four BCs are considered, one corresponding to a load and three to displacements.
 the first BC introduced is a load, then the coordinate system, loadfactor time function and base load vector are defined
```
boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) t ;
boundaryConds(1).loadsBaseVals = [ p 0 0 0 0 0 ] ;
```
 the other BCs have imposed displacements
```
boundaryConds(2).imposDispDofs = [1] ;
boundaryConds(2).imposDispVals =  0  ;
%
boundaryConds(3).imposDispDofs = [3] ;
boundaryConds(3).imposDispVals =  0  ;
%
boundaryConds(4).imposDispDofs = [5] ;
boundaryConds(4).imposDispVals =  0  ;
%
```

#### initialConds
 since no initial non-homogeneous initial conditions are used, an empty struct is used .
```
initialConds = struct();
```

### Mesh
 A simple hand-made 8-node mesh, with 6 tetrahedrons is considered

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/solidCubeMeshHTML.svg" alt="mesh diagram" width="500"/>
```
```@raw latex
\begin{center}
\def\svgwidth{0.6\textwidth}
\input{solidCubeMeshPDF.pdf_tex}
\end{center}
```
 The node coordinates matrix is given by the following
```
mesh.nodesCoords = [ 0    0    0 ; ...
                     0    0   Lz ; ...
                     0   Ly   Lz ; ...
                     0   Ly    0 ; ...
                     Lx   0    0 ; ...
                     Lx   0   Lz ; ...
                     Lx  Ly   Lz ; ...
                     Lx  Ly    0 ] ;
```
 and the connectivity cell is defined as follows with the four MEBI parameters for each element followed by the indexes of the nodes of each element. All the eight triangle elements are considered with no material (since they are used only to include load) and the following six elements are solid SVK material tetrahedrons.
```
mesh.conecCell = {[ 0 1 1 0    5 8 6   ]; ... % loaded face
                  [ 0 1 1 0    6 8 7   ]; ... % loaded face
                  [ 0 1 2 0    4 1 2   ]; ... % x=0 supp face
                  [ 0 1 2 0    4 2 3   ]; ... % x=0 supp face
                  [ 0 1 3 0    6 2 1   ]; ... % y=0 supp face
                  [ 0 1 3 0    6 1 5   ]; ... % y=0 supp face
                  [ 0 1 4 0    1 4 5   ]; ... % z=0 supp face
                  [ 0 1 4 0    4 8 5   ]; ... % z=0 supp face
                  [ 1 2 0 0    1 4 2 6 ]; ... % tetrahedron
                  [ 1 2 0 0    6 2 3 4 ]; ... % tetrahedron
                  [ 1 2 0 0    4 3 6 7 ]; ... % tetrahedron
                  [ 1 2 0 0    4 1 5 6 ]; ... % tetrahedron
                  [ 1 2 0 0    4 6 5 8 ]; ... % tetrahedron
                  [ 1 2 0 0    4 7 6 8 ]  ... % tetrahedron
                } ;
```

### Analysis parameters

```
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30     ;
analysisSettings.stopTolDeltau = 1.0e-8 ;
analysisSettings.stopTolForces = 1.0e-8 ;
analysisSettings.finalTime      = 1      ;
analysisSettings.deltaT        = .125   ;
```

### Output parameters
```
otherParams.plotsFormat = 'vtk' ;
otherParams.problemName = 'uniaxialExtension_HandMadeMesh' ;
```

```
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
```

### Analytic solution computation
```
analyticFunc = @(w) 1/p *E * 0.5 * ( ( 1 + w/Lx ).^3 - ( 1 + w/Lx) ) ;
```

```
analyticCheckTolerance = 1e-6 ;
analyticFunc           = @(w) 1/p * E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) ) ;
controlDisps = matUs(6*6+1,:) ;
analyticVals = analyticFunc( controlDisps ) ;
controlDispsValsCase1         = controlDisps  ;
loadFactorAnalyticalValsCase1 = analyticVals  ;
loadFactorNumericalValsCase1  = loadFactorsMat ;
```

## Numerical solution: case 2
In this analysis case, the mesh information is read from a gmsh-generated
mesh file, the pressure is applied using local coordinates and the stiffness
 matrix is computed using the complex-step method.

```
otherParams.problemName = 'uniaxialExtension_GMSH_ComplexStep' ;

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;
boundaryConds(1).loadsCoordSys = 'local';
boundaryConds(1).loadsBaseVals = [0 0 0 0 p 0 ] ;
elements(2).elemTypeParams = [ 2 ] ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

controlDisps = matUs(6*6+1,:) ;
analyticVals = analyticFunc( controlDisps ) ;
controlDispsValsCase2         = controlDisps  ;
loadFactorAnalyticalValsCase2 = analyticVals  ;
loadFactorNumericalValsCase2  = loadFactorsMat ;

aux1 = loadFactorNumericalValsCase1' - loadFactorAnalyticalValsCase1 ;
aux2 = loadFactorNumericalValsCase2' - loadFactorAnalyticalValsCase1 ;

verifBoolean = ...
     ( norm( aux1 ) / norm( loadFactorNumericalValsCase1 ) < analyticCheckTolerance ) ...
  && ( norm( aux2 ) / norm( loadFactorNumericalValsCase1 ) < analyticCheckTolerance ) ;
```


## Plot
The numerical and analytic solutions are plotted.
```
lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure, hold on, grid on
plot( controlDispsValsCase1, loadFactorAnalyticalValsCase1, 'r-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsValsCase1, loadFactorNumericalValsCase1,  'k-o' , 'linewidth', lw,'markersize',ms )
plot( controlDispsValsCase2, loadFactorNumericalValsCase2,  'g-s' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('\lambda(t)') ;
legend( 'Analytic', 'Numeric-1', 'Numeric-2', 'location', 'North' )
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/verifUniaxial.png','-dpng')
```

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/verifUniaxial.png" alt="plot check" width="500"/>
```
