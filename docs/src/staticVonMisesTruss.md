# Static Von Mises Truss example
---

In this tutorial, the Static Von Mises Truss example and its resolution using ONSAS are described. The aim of this example is to validate the Newton-Raphson method implementation by comparing the results provided by ONSAS with analytic solutions. The structural model is formed by two truss elements with length $L$ as it is shown in the figure, with node $2$ submitted to a nodal load $P$ and restrained to move in the $x-z$ plane and nodes $1$ and $3$ fixed.

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/vonMisesTruss.svg" alt="structure diagram" width="500"/>
```
## Analytic solution
--------------------

The solutions are developed in [section 2.3 of (Bazzano and Pérez Zerpa, 2017)](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf#section.2.3). The expressions obtained for different strain measures are:
 * Rotated-Engineering: $P = \frac{EA_o(z_2+w)\left(\sqrt{(w+z_2)^2+x_2^2}-l_o\right)}{l_o\sqrt{(w+z_2)^2+x_2^2}}$
 * SVK: $P = \frac{EA_o (z_2+w)\left( 2 z_2 w + w^2 \right) }{ 2 l_o^3 }$
 where $x_2$ and $z_2$ are the coordinates of node 2 and $w$ is measured positive as $z$.

## Numerical solution
---------------------


The Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/staticVonMisesTruss/onsasExample_staticVonMisesTruss.m).

Before defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar auxiliar parameters are defined.
```
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ; nu = 0 ;
% x and z coordinates of node 2
x2 = cos( ang1*pi/180 ) * L ;
z2 = sin( ang1*pi/180 ) * L ;
```

### MEBI parameters

The modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.

#### materials
 Since both bars are formed by the same material only one `materials` struct is defined. The constitutive behavior considered in the first analysis case is the Rotated Engineering strain, then the field `hyperElasModel` is set to:
```
materials.hyperElasModel  = '1DrotEngStrain' ;
```
 and in the field `hyperElasParams` a vector with the parameters of the Engineering Strain model is set
```
materials.hyperElasParams = [ E nu ] ;
```
 which in the case of this model are the Young modulus and the Poisson ratio.
 The field `density` is not set, then the default $\rho = 0$ value is considered by ONSAS.

#### elements

Two different types of elements are required to create the model: `node` and `truss`, thus, the `elements` struct will have two entries. The type of the first entry is
```
elements(1).elemType = 'node' ;
```
 and the second entry is
```
elements(2).elemType = 'truss';
```
 for the geometries, the node has no geometry to assign, and the truss elements will be set as a square-cross section, then the elemTypeGeometry field is:
```
elements(2).elemCrossSecParams = { } ;
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = sqrt(A*4/pi) ;
elements(2).elemTypeParams = 1 ;
```

#### boundaryConds

 The elements are submitted to two different BoundaryConditions, then the struct `boundaryConds` will have length two.
 The nodes $1$ and $3$ are fixed, without loads applied (this is the first BC), and node $2$ has a constraint in displacement and an applied load (second BC).
 For the displacements, the first BC corresponds to a xyz-fixed displacement,

```
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
```
and the second BC corresponds to a zero displacement only in the $y$ direction.
```
boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
```
 Regarding the loads, the second BC is set so that the final load factor is $3 \cdot 10^8$ at 1 second. The default zero density is used, then no inertial effects are considered.

```
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 3.0e8*t     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ] ;
```

#### initial Conditions
 homogeneous initial conditions are considered, then an empty struct is set:
```
initialConds                = struct() ;
```

### mesh parameters

The coordinates of the nodes of the mesh are given by the matrix:
```
mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;
```
 where the columns 1,2 and 3 correspond to $x$, $y$ and $z$ coordinates, respectively, and the row $i$-th corresponds to the coordinates of node $i$.

The connectivity is introduced using the _conecCell_ cell. Each entry of the cell (indexed using {}) contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
```
mesh.conecCell = { } ;
```
 Then the entry of node $1$ is introduced:
```
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
```
 the first MEBI parameter (Material) is set as _zero_ (since nodes dont have material). The second parameter corresponds to the Element, and a _1_ is set since `node` is the first entry of the  `elements.elemType` cell. For the BC index, we consider that node $1$ is fixed, then the first index of the `boundaryConds` struct is used. Finally, no specific initial conditions are set for the node (0) and at the end of the vector the number of the node is included (1).
 A similar approach is used for node $3$,
```
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  3   ] ;
```
 and for node $2$ only the boundary condition is changed.
```
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
```
 Regarding the truss elements, the first material is considered, the second type of element, and no boundary conditions are applied.
```
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0 0  2 3 ] ;
```

### analysisSettings
 The method used in the analysis is the Newton-Raphson, then the field `methodName` must be introduced as:
```
analysisSettings.methodName    = 'newtonRaphson' ;
```
 and the following parameters correspond to the iterative numerical analysis settings
```
analysisSettings.deltaT        =   0.1  ;
analysisSettings.finalTime      =   1    ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
```

### otherParams
```
otherParams.problemName = 'staticVonMisesTruss_NR_RotEng';
otherParams.plotsFormat = 'vtk' ;
otherParams.controlDofs = [2 5 ];
```

### Analysis case 1: NR with Rotated Eng Strain
 In the first case ONSAS is run and the solution at the dof of interest is stored.
```
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNREngRot =  -matUs(11,:) ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;
```

 and the analytical value of the load factors is computed, as well as its difference with the numerical solution
```
analyticLoadFactorsNREngRot = @(w) -2 * E*A* ...
     ( (  (z2+(-w)).^2 + x2^2 - L^2 ) ./ (L * ( L + sqrt((z2+(-w)).^2 + x2^2) )) ) ...
     .*  (z2+(-w))                    ./ ( sqrt((z2+(-w)).^2 + x2^2) )  ;
difLoadEngRot = analyticLoadFactorsNREngRot( controlDispsNREngRot)' - loadFactorsNREngRot ;
```

### Analysis case 2: NR with Green Strain
 In order to perform a SVK case, the material is changed and the problemName is also updated
```
otherParams.problemName = 'staticVonMisesTruss_NR_Green';
materials.hyperElasModel  = 'SVK' ;
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasParams = [ lambda mu ] ;
```
 the load history is also changed
```
boundaryConds(2).loadsTimeFact = @(t) 1.5e8*t ;
```
 and the analysis is run
```
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRGreen =  -matUs(11,:) ;
loadFactorsNRGreen  =  loadFactorsMat(:,2) ;
% %md the analytic solution is computed
analyticLoadFactorsGreen = @(w) - 2 * E*A * ( ( z2 + (-w) ) .* ( 2*z2*(-w) + w.^2 ) ) ./ ( 2.0 * L^3 )  ;
difLoadGreen = analyticLoadFactorsGreen( controlDispsNRGreen )' - loadFactorsNRGreen ;
```

### Analysis case 3: NR-AL with Green Strain
```

elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ sqrt(A) sqrt(A)] ;

```
 In this case, the numerical method is changed for newtonRaphson arc length.
```
otherParams.problemName       = 'staticVonMisesTruss_NRAL_Green' ;
analysisSettings.methodName   = 'arcLength'                      ;
analysisSettings.finalTime     = 4.5                             ;
analysisSettings.incremArcLen = 0.1                             ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.posVariableLoadBC = 2 ;
```

```
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRALGreen =  -matUs(11,:) ;
loadFactorsNRALGreen  =  loadFactorsMat(:,2) ;
analyticLoadFactorsNRALGreen = analyticLoadFactorsGreen(controlDispsNRALGreen);
difLoadGreenNRAL = analyticLoadFactorsNRALGreen' - loadFactorsNRALGreen ;
```
## Verification
 the numerical resolution is validated for both strain measures.
```

verifBoolean =  ( ( norm( difLoadEngRot    ) / norm( loadFactorsNREngRot  ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreen     ) / norm( loadFactorsNRGreen   ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreenNRAL ) / norm( loadFactorsNRALGreen ) ) <  1e-4 ) ;

```
### Plots
 and solutions are plotted.
```
lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure
plot( controlDispsNREngRot, analyticLoadFactorsNREngRot( controlDispsNREngRot) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNREngRot, loadFactorsNREngRot, 'k-o' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, analyticLoadFactorsGreen( controlDispsNRALGreen ), 'g-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRGreen, loadFactorsNRGreen, 'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, loadFactorsNRALGreen, 'c-^' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement w(t)');   laby = ylabel('\lambda(t)') ;
legend( 'analytic-RotEng', 'NR-RotEng','analytic-Green', 'NR-Green','NRAL-Green', 'location','SouthEast')
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/vonMisesTrussCheck.png','-dpng')
```

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/vonMisesTrussCheck.png" alt="plot check" width="500"/>
```

