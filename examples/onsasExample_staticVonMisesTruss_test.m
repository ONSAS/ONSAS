% ======================================================================
% Von Mises truss example
clear all, close all
dirOnsas    = [ pwd '/..' ] ;
problemName = 'staticVonMisesTrussNR' ;
% ----------------------------------------------------------------------
% scalar auxiliar parameters
E = 210e9 ;  A = 2.5e-4 ; ang1 = 65 ; L = 2 ; nu = 0 ;  rho = 0 ; 

% ----------------------------------------------------------------------
% MELCS parameters
% Materials
materialsParams = {[ rho 2 E nu ]} ;
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
auxx = cos( ang1*pi/180 ) * L ;        auxy = sin( ang1*pi/180 ) * L ;
% nodes matrix
Nodes = [      0  0     0  ; ...
            auxx  0  auxy  ; ...
          2*auxx  0     0  ] ;

% connectivity cell
Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 2  2   ] ; ... % loaded node
          [ 0 1 0 0 1  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % truss element
          [ 1 2 0 1 0  2 3 ] } ;   % truss element
         
% ----------------------------------------------------------------------
% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 1.5e7  ;    nLoadSteps       = 6      ;
stopTolIts       = 30     ;
numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 
stabilityAnalysisBoolean = 2 ;

% ----------------------------------------------------------------------
% analysis parameters
controlDofs      = [ 2 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ];

addpath( dirOnsas );
ONSAS;


return


targetLoadFactrNRAL = 4e7    ; % arc length
incremArcLen     = .2     ;

%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactrNR nLoadSteps ] ; 

%% Analysis parameters


% analytical solution using engineering strain
analyticSolFlag        = 2    ;
analyticCheckTolerance = 1e-4 ;
l0 = sqrt(auxx^2 + auxy^2) ;
analyticFunc = @(w) -2 * E*A* ( (  (auxy+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxy+(-w)).^2 + auxx^2) )) ) ...
 .* (auxy+(-w)) ./ ( sqrt((auxy+(-w)).^2 + auxx^2) )  ; 

%% Output parameters
printFlag = 0 ;
%~ plotParamsVector = [ 3 10];
plotParamsVector = [ 0 ];

sectPar = [12 .1 .1]

reportBoolean = 0 ;

%% ONSAS execution
ONSAS

% ======================================================================

Conec = {[ 0 1 0 0 1  1   ] ; ... % fixed node
         [ 0 1 1 0 2  2   ] ; ... % loaded node
         [ 0 1 0 0 1  3   ] ; ... % fixed node
         [ 1 2 0 1 0  1 2 ] ; ... % truss element
         [ 1 2 0 1 0  2 3 ]   ... % truss element
         } ;
         
         
controlDispsNRAL = controlDisps ;
loadFactorsNRAL  = loadFactors ;
analyticNRAL = analyticVals ;

nLoadSteps       = 6    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNR nLoadSteps ] ; 
%~ plotParamsVector = [ 3 ];
plotParamsVector = [ 0 ];
problemName = 'staticVonMisesTrussNR' ;

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
legend('analytic','NRAL','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print( [ 'vonmises' ] ,'-dpdflatex','-tight') ;
cd(acdir);

