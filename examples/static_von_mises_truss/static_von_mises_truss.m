%md# Static Von-Mises Truss example
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/staticVonMisesTruss/onsasExample_staticVonMisesTruss.m)
%md
%mdIn this example the Static Von Mises Truss problem and its resolution using ONSAS are described. The aim of this example is to validate the implementations of the Newton-Raphson and Newton-Raphson-Arc-Length methods by comparing the results provided with the analytic solutions.
%md
%mdThe structural model is formed by two truss elements with length $L$ as it is shown in the figure, with node $2$ submitted to a nodal load $P$ and restrained to move in the $x-z$ plane, and nodes $1$ and $3$ fixed.
%md
%md```@raw html
%md<img src="../../assets/vonMisesTruss.svg" alt="von-mises truss structure diagram" width="500"/>
%md```
%md
%md## Analytic solutions
%md
%mdThe solutions for the nonlinear cases are developed in section 2.3 of [(Bazzano and Pérez Zerpa, 2017)](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf#section.2.3). The expressions obtained for different strain measures are:
%md * Rotated-Engineering: $P = \frac{EA_o(z_2+w)\left(\sqrt{(w+z_2)^2+x_2^2}-l_o\right)}{l_o\sqrt{(w+z_2)^2+x_2^2}}$
%md * SVK: $P = \frac{EA_o (z_2+w)\left( 2 z_2 w + w^2 \right) }{ 2 l_o^3 }$
%md where $x_2$ and $z_2$ are the coordinates of node 2 and $w$ is the vertical displacement, measured positive as $z$.
%md
%md## Numerical solutions
%md
%mdBefore defining the structs, the workspace is cleared, the ONSAS directory is added to the path and scalar auxiliar parameters are defined.
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ; nu = 0 ;

% x and z coordinates of node 2
x2 = cos( ang1*pi/180 ) * L ;
z2 = sin( ang1*pi/180 ) * L ;
%md
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md#### materials
%md Since both bars are formed by the same material only one `materials` struct is defined. The constitutive behavior considered in the first analysis case is the Rotated Engineering strain, then the field `hyperElasModel` is set to:
materials.hyperElasModel  = '1DrotEngStrain' ;
%md and in the field `hyperElasParams` a vector with the parameters of the Engineering Strain model is set
materials.hyperElasParams = [ E nu ] ;
%md which in the case of this model are the Young modulus and the Poisson ratio.
%md The field `density` is not set, then the default $\rho = 0$ value is considered by ONSAS.
%md
%md#### elements
%md
%mdTwo different types of elements are required to create the model: `node` and `truss`, thus, the `elements` struct will have two entries. The type of the first entry is
elements(1).elemType = 'node' ;
%md and the second entry is
elements(2).elemType = 'truss';
%md for the geometries, the node has no geometry to assign, and the truss elements will be set as a circle cross-section, then the elemCrossSecParams field is:
elements(2).elemCrossSecParams = { 'circle' , sqrt(A*4/pi) } ;
%md
%md#### boundaryConds
%md
%md The elements are submitted to two different BoundaryConditions, then the struct `boundaryConds` will have length two.
%md The nodes $1$ and $3$ are fixed, without loads applied (this is the first BC), and node $2$ has a constraint in displacement and an applied load (second BC).
%md For the displacements, the first BC corresponds to a xyz-fixed displacement,
%md
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
%mdand the second BC corresponds to a zero displacement only in the $y$ direction.
boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
%md Regarding the loads, the second BC is set so that the final load factor is $3 \cdot 10^8$ at 1 second. The default zero density is used, then no inertial effects are considered.
%md
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 3.0e8*t     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ] ;
%md
%md#### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;
%md
%md### mesh parameters
%md
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;
%md where the columns 1,2 and 3 correspond to $x$, $y$ and $z$ coordinates, respectively, and the row $i$-th corresponds to the coordinates of node $i$.
%md
%mdThe connectivity is introduced using the _conecCell_ cell. Each entry of the cell (indexed using {}) contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = cell(5,1) ;
%md Then the entry of node $1$ is introduced:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the first MEBI parameter (Material) is set as _zero_ (since nodes dont have material). The second parameter corresponds to the Element, and a _1_ is set since `node` is the first entry of the  `elements.elemType` cell. For the BC index, we consider that node $1$ is fixed, then the first index of the `boundaryConds` struct is used. Finally, no specific initial conditions are set for the node (0) and at the end of the vector the number of the node is included (1).
%md A similar approach is used for node $3$,
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  3   ] ;
%md and for node $2$ only the boundary condition is changed.
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
%md Regarding the truss elements, the first material is considered, the second type of element, and no boundary conditions are applied.
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0 0  2 3 ] ;
%md
%md### analysisSettings
%md The method used in the analysis is the Newton-Raphson, then the field `methodName` must be introduced as:
analysisSettings.methodName    = 'newtonRaphson' ;
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   0.1  ;
analysisSettings.finalTime     =   1    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;
%md
%md### otherParams
otherParams.problemName = 'staticVonMisesTruss_NR_RotEng';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;
%md
%md### Analysis case 1: Newton-Raphson with Rotated Eng Strain
%md In the first case ONSAS is run and the solution at the dof of interest is stored.
[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNREngRot =  -matUs(11,:) ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;
%md
%md### Analysis case 2: Newton-Raphson with linear elastic behavior
%mdIn this case a linear elastic behavior is assumed. Then the hyperelasmodel es overwritten
materials.hyperElasModel = 'linearElastic' ;
otherParams.problemName  = 'staticVonMisesTruss_linearElastic';
analysisSettings.finalTime  =   1.5    ;
%md and the analysis is run again
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md the displacements are extracted
controlDispsNRLinearElastic =  -matUs(11,:) ;
loadFactorsNRLinearElastic  =  loadFactorsMat(:,2) ;
%md
%md and the analytic values of the load factor are computed, as well as its difference with the numerical solution
%md
analyticLoadFactorsNREngRot = @(w) -2 * E*A* ...
     ( (  (z2+(-w)).^2 + x2^2 - L^2 ) ./ (L * ( L + sqrt((z2+(-w)).^2 + x2^2) )) ) ...
     .*  (z2+(-w))                    ./ ( sqrt((z2+(-w)).^2 + x2^2) )  ;
difLoadEngRot = analyticLoadFactorsNREngRot( controlDispsNREngRot)' - loadFactorsNREngRot ;
%md
%md### Analysis case 3: NR with Green Strain
%md In order to perform a SVK case analysis, the material is changed and the problemName is also updated
otherParams.problemName  = 'staticVonMisesTruss_NR_Green';
materials.hyperElasModel = 'SVK' ;
analysisSettings.finalTime =   1.0    ;
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasParams = [ lambda mu ] ;
%md the load history is also changed
boundaryConds(2).loadsTimeFact = @(t) 1.5e8*t ;
%boundaryConds(2).userLoadsFilename = 'myVMLoadFunc' ;
%md and the analysis is run
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md and the displacements are extracted
controlDispsNRGreen =  -matUs(11,:) ;
loadFactorsNRGreen  =  loadFactorsMat(:,2) ;
%md the analytic solution is computed
analyticLoadFactorsGreen = @(w) - 2 * E*A * ( ( z2 + (-w) ) .* ( 2*z2*(-w) + w.^2 ) ) ./ ( 2.0 * L^3 )  ;
difLoadGreen = analyticLoadFactorsGreen( controlDispsNRGreen )' - loadFactorsNRGreen ;
%md
%md### Analysis case 4: NR-AL with Green Strain
%md
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ sqrt(A) sqrt(A)] ;
%mdThe same loading conidition as before is used, but given by a user load function. The argument set in this case is:
%md
%md In this case, the numerical method is changed for newtonRaphson arc length.
otherParams.problemName       = 'staticVonMisesTruss_NRAL_Green' ;
analysisSettings.methodName   = 'arcLength'                      ;
analysisSettings.finalTime    = 1                               ;
analysisSettings.incremArcLen = 0.15                             ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.posVariableLoadBC = 2 ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRALGreen =  -matUs(11,:) ;
loadFactorsNRALGreen  =  loadFactorsMat(:,2) ;
analyticLoadFactorsNRALGreen = analyticLoadFactorsGreen(controlDispsNRALGreen);
difLoadGreenNRAL = analyticLoadFactorsNRALGreen' - loadFactorsNRALGreen ;
%md## Verification
%md the numerical resolution is validated for both strain measures.
%md
verifBoolean =  ( ( norm( difLoadEngRot    ) / norm( loadFactorsNREngRot  ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreen     ) / norm( loadFactorsNRGreen   ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreenNRAL ) / norm( loadFactorsNRALGreen ) ) <  1e-4 ) ;
%md
%md### Plots
%md and solutions are plotted.
lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure
plot( controlDispsNREngRot, analyticLoadFactorsNREngRot( controlDispsNREngRot) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNREngRot, loadFactorsNREngRot, 'k-o' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, analyticLoadFactorsGreen( controlDispsNRALGreen ), 'g-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRGreen, loadFactorsNRGreen, 'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, loadFactorsNRALGreen, 'c-^' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRLinearElastic, loadFactorsNRLinearElastic, 'm-+' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement w(t)');   laby = ylabel('\lambda(t)') ;
legend( 'analytic-RotEng', 'NR-RotEng','analytic-Green', 'NR-Green','NRAL-Green','LinearElastic', 'location','northoutside')
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/vonMisesTrussCheck.png','-dpng')
%print('../../docs/src/assets/vonMisesTrussCheck.png','-dpng')
%md
%md```@raw html
%md<img src="../../assets/vonMisesTrussCheck.png" alt="plot check" width="500"/>
%md```
%md
