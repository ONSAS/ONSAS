% ======================================================================
% Von Mises truss example

clear all, close all
dirOnsas    = [ pwd '/..' ] ; addpath( dirOnsas );
problemName = 'staticVonMisesTrussLin' ;
% ----------------------------------------------------------------------
% scalar auxiliar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ; nu = 0 ;  rho = 0 ; 

% ----------------------------------------------------------------------
% MELCS parameters

% Materials
materialsParams = {[ rho 1 E nu ]} ;

% Elements
elementsParams  = { 1; 2} ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 -1 0] } ;

% Cross-Sections
crossSecsParams = { [ 2 sqrt(A) sqrt(A) ] } ;

% Springs
springsParams   = { [ inf  0  inf  0  inf   0 ] ; ...
                    [ 0    0  inf  0    0   0 ] } ;

% ----------------------------------------------------------------------
% nodes coordinates matrix and connectivity cell
auxx = cos( ang1*pi/180 ) * L ;        auxz = sin( ang1*pi/180 ) * L ;

% node coordinates matrix
Nodes = [      0  0     0  ; ...
            auxx  0  auxz  ; ...
          2*auxx  0     0  ] ;

% connectivity cell
Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 2  2   ] ; ... % loaded node
          [ 0 1 0 0 1  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % truss element
          [ 1 2 0 1 0  2 3 ] } ;   % truss element

controlDofs      = [ 2 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;

analyticSolFlag = 2    ;
analyticFunc    = @(w) 2 * E * A * sin(ang1*pi/180)^2 * w / L ;

ONSAS

% ----------------------------------------------------------------------
% connectivity cell
Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 2  2   ] ; ... % loaded node
          [ 0 1 0 0 1  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % truss element
          [ 1 2 0 1 0  2 3 ] } ;   % truss element

problemName = 'staticVonMisesTrussNR' ;

materialsParams = {[ rho 3 E nu ]} ;

% Cross-Sections
crossSecsParams = { [ 3 sqrt(A*4/pi) ] } ;

% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 1.5e8  ;    nLoadSteps       = 6      ;
stopTolIts       = 30     ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

stabilityAnalysisBoolean = 2 ;

analyticSolFlag = 2    ;
l0           = sqrt(auxx^2 + auxz^2) ;
analyticFunc = @(w) -2 * E*A* ( (  (auxz+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxz+(-w)).^2 + auxx^2) )) ) ...
 .* (auxz+(-w)) ./ ( sqrt((auxz+(-w)).^2 + auxx^2) )  ; 

ONSAS

return
% arc length params
targetLoadFactrNRAL   = 4e7  ;
incremArcLen          = 0.2  ;
numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 

materialsParams = {[ rho 3 E nu ]} ;

% analytical solution using engineering strain
analyticSolFlag = 2    ;
analyticFunc    = @(w) ...
  -2 * E*A* ( (  (auxz+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxz+(-w)).^2 + auxx^2) )) ) ...
 .* (auxz+(-w)) ./ ( sqrt((auxz+(-w)).^2 + auxx^2) )  ; 

plotParamsVector = [ 3 ];

%% ONSAS execution
ONSAS

return
% ======================================================================

%~ Conec = {[ 0 1 0 0 1  1   ] ; ... % fixed node
         %~ [ 0 1 1 0 2  2   ] ; ... % loaded node
         %~ [ 0 1 0 0 1  3   ] ; ... % fixed node
         %~ [ 1 2 0 1 0  1 2 ] ; ... % truss element
         %~ [ 1 2 0 1 0  2 3 ]   ... % truss element
         %~ } ;
       %~ Conec
         
[ Nodes, Conec ] = meshFileReader( 'vonMises.dxf' ) ;

         
controlDispsNRAL = controlDisps ;
loadFactorsNRAL  = loadFactors ;
analyticNRAL = analyticVals ;

nLoadSteps       = 6    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNR nLoadSteps ] ; 
%~ plotParamsVector = [ 3 ];
plotParamsVector = [ 0 ];
problemName = 'staticVonMisesTrussNR_DXF' ;

acdir = pwd ; cd(dirOnsas); ONSAS; cd(acdir) ;

controlDispsNR = controlDisps ;
loadFactorsNR  = loadFactors ;


lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;

figure
plot( controlDispsNRAL, analyticNRAL ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNRAL, loadFactorsNRAL,'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNR, loadFactorsNR,'k-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic','NRAL','NR-DXF','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print( [ 'vonmises' ] ,'-dpdflatex','-tight') ;
cd(acdir);

