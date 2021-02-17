% ======================================================================
% Von Mises truss example

close all, clear all        ;

dirOnsas = [ pwd '/../../src' ] ; % set ONSAS.m directory
addpath( dirOnsas )             ; % add ONSAS directory to path

% ===========================
% first case: linear analysis

problemName = 'staticVonMisesTrussLin' ; %#ok

% ----------------------------------------------------------------------
% scalar auxiliar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ; nu = 0 ;  rho = 0 ; 

% ----------------------------------------------------------------------
% Material Elements Boundary-condtions Initial-conditions parameters

% Materials properties
materials                  = struct()   ;  
materials.hyperElasModel   = {'SVK'}    ;
materials.hyperElasParams  = {[ E nu ]} ;

% Elements properties
elements                   = struct()   ;
elements.elemType          = { 'node', 'truss' } ;
elements.elemTypeParams    = {   []  ,    []   } ;
elements.elemTypeGeometry  = {   []  , [ 2 sqrt(A) sqrt(A) ] } ;

% Boundary Conditions parameters
boundaryConds               = struct()                  ;
boundaryConds.loadCoordSys = {    []    , 'global'   } ;
boundaryConds.loadTimeFact = {    []    , @(t) t     } ;
boundaryConds.loadBaseVals = {    []    , [ 0 0  0 0 -1 0 ] } ;
boundaryConds.impoDispDofs = { [ 1 3 5] ,   3        } ;
boundaryConds.impoDispVals = { [ 0 0 0] ,   0        } ;

% Initial Conditions parameters
initialConds                = struct() ;

% ----------------------------------------------------------------------
% mesh parameters

% nodes coordinates matrix and connectivity cell
auxx = cos( ang1*pi/180 ) * L ;        auxz = sin( ang1*pi/180 ) * L ;

mesh = struct();

% node coordinates matrix
mesh.Nodes = [      0  0     0  ; ...
            auxx  0  auxz  ; ...
          2*auxx  0     0  ] ;

% connectivity cell: MEBI + nodes
mesh.Conec = { [ 0 1 1 0  1   ] ; ... % fixed node
               [ 0 1 2 0  2   ] ; ... % loaded node
               [ 0 1 1 0  3   ] ; ... % fixed node
               [ 1 2 0 0  1 2 ] ; ... % truss element
               [ 1 2 0 0  2 3 ] } ;   % truss element

solVerifData = struct() ;
solVerifData.controlDofs = [ 2 5 -1 ] ; % [ node nodaldof scalefactor ]
solVerifData.analyticSolFlag = 2    ;
solVerifData.analyticFunc    = @(w) 2 * E * A * sin(ang1*pi/180)^2 * w / L ;

generalParams = struct() ;
generalParams.plotParamsVector = [3] ;
generalParams.reportBoolean    = 1     ;
generalParams.problemName      = problemName ;

numericalMethod = struct();
numericalMethod.name = 'linearStep';

ONSAS( mesh, materials, elements, boundaryConds, initialConds, generalParams, solVerifData, numericalMethod )
% ----------------------------------------------------------------------

return

% ====================================
% second case: newton raphson analysis

% ----------------------------------------------------------------------
% redefine Connectivity cell
Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 2  2   ] ; ... % loaded node
          [ 0 1 0 0 1  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % truss element
          [ 1 2 0 1 0  2 3 ] } ;   % truss element

problemName = 'staticVonMisesTrussNR' ;

% material change
materialsParams = {[ rho 3 E nu ]} ;

% Cross-Sections
crossSecsParams = { [ 3 sqrt(A*4/pi) ] } ;

% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 2.0e8  ;    nLoadSteps       = 6      ;
stopTolIts       = 30     ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

stabilityAnalysisBoolean = 2 ;

l0           = sqrt(auxx^2 + auxz^2) ;
analyticFunc = @(w) -2 * E*A* ( (  (auxz+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxz+(-w)).^2 + auxx^2) )) ) ...
 .* (auxz+(-w)) ./ ( sqrt((auxz+(-w)).^2 + auxx^2) )  ; 

ONSAS

controlDispsNR = controlDisps ;
loadFactorsNR  = loadFactors ;
% ----------------------------------------------------------------------


% ===============================================
% third case: NRarc-length analysis with dxf mesh

% ----------------------------------------------------------------------
problemName = 'staticVonMisesTrussNRAL_DXF' ;
[ Nodes, Conec ] = meshFileReader( 'geometry_vonMises.dxf' ) ;

% arc length params
targetLoadFactrNRAL   = 4e8  ;
incremArcLen          = 0.2  ;
numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 

% analytical solution using engineering strain
analyticFunc    = @(w) ...
  -2 * E*A* ( (  (auxz+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxz+(-w)).^2 + auxx^2) )) ) ...
 .* (auxz+(-w)) ./ ( sqrt((auxz+(-w)).^2 + auxx^2) )  ; 

ONSAS

controlDispsNRAL = controlDisps ;
loadFactorsNRAL  = loadFactors  ;
analyticNRAL     = analyticVals ;
% ----------------------------------------------------------------------


% ===============================================
% methods comparison

% ----------------------------------------------------------------------
% --- plots ---
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;

figure
plot( controlDispsNRAL, analyticNRAL ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNRAL, loadFactorsNRAL,'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNR, loadFactorsNR,'k-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic','NRAL-DXF','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
