%md# Uniaxial Extension Solid example using Cpp solver

clear all, close all
addpath( genpath( [ pwd '/../../src'] ) ) ;

% scalar parameters
E = 1 ; nu = 0.3 ; p = 3 ; Lx = 2 ; Ly = 1 ; Lz = 1 ;

lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasModel = 'SVK' ;
materials.hyperElasParams = [ lambda mu ] ;

elements(1).elemType = 'triangle';
elements(2).elemType = 'tetrahedron' ;
%md
boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) t ;
boundaryConds(1).loadsBaseVals = [p 0 0 0 0 0 ] ;
%
boundaryConds(2).imposDispDofs = [1] ;
boundaryConds(2).imposDispVals = [0] ;
%
boundaryConds(3).imposDispDofs = [3] ;
boundaryConds(3).imposDispVals = [0] ;
%
boundaryConds(4).imposDispDofs = [5] ;
boundaryConds(4).imposDispVals = [0] ;
%md
initialConds = struct();

mesh.nodesCoords = [ 0    0    0 ; ...
                     0    0   Lz ; ...
                     0   Ly   Lz ; ...
                     0   Ly    0 ; ...
                     Lx   0    0 ; ...
                     Lx   0   Lz ; ...
                     Lx  Ly   Lz ; ...
                     Lx  Ly    0 ] ;

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

%md### Analysis parameters
%md
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30     ;
analysisSettings.stopTolDeltau = 1.0e-8 ;
analysisSettings.stopTolForces = 1.0e-8 ;
analysisSettings.finalTime      = 1      ;
analysisSettings.deltaT        = .125   ;
analysisSettings.solverLang    = 'C++'  ;

%md### Output parameters
otherParams.plotsFormat = 'vtk' ;
otherParams.problemName = 'uniaxialExtension_HandMadeMesh' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md### Analytic solution computation
analyticFunc = @(w) 1/p *E * 0.5 * ( ( 1 + w/Lx ).^3 - ( 1 + w/Lx) ) ;
%md
analyticCheckTolerance = 1e-6 ;
analyticFunc           = @(w) 1/p * E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) ) ;
controlDisps = matUs(6*6+1,:) ;
analyticVals = analyticFunc( controlDisps ) ;
controlDispsValsCase1         = controlDisps  ;
loadFactorAnalyticalValsCase1 = analyticVals  ;
loadFactorNumericalValsCase1  = loadFactorsMat ;

%md
%md## Plot
%mdThe numerical and analytic solutions are plotted.
lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure, hold on, grid on
plot( controlDispsValsCase1, loadFactorAnalyticalValsCase1, 'r-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsValsCase1, loadFactorNumericalValsCase1,  'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('\lambda(t)') ;
legend( 'Analytic', 'Numeric-1', 'location', 'North' )
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/verifUniaxial.png','-dpng')
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/verifUniaxial.png" alt="plot check" width="500"/>
%md```
