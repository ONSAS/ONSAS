%% Von Mises truss example using Newton-Raphson Arc-Length Method
%
%%
clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'staticVonMisesTrussNRAL' ;

auxx = cos(65*pi/180) * 2 ;
auxy = sin(65*pi/180) * 2 ;
imperfPerc = .0 ;

Nodes = [      0  0     0  ; ...
            auxx*(1+imperfPerc)  0  auxy  ; ...
          2*auxx  0     0  ] ;

Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 2  2   ] ; ... % loaded node
          [ 0 1 0 0 1  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % truss element
          [ 1 2 0 1 0  2 3 ]   ... % truss element
         } ;

% ======================================================================
% --- MELCS parameters ---
materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

E = 210e9 ;  nu = 0 ;  rho = 0 ;
materialsParams = {[ rho 3 E nu ]} ;

elementsParams = { 1; 2} ;

loadsParams = { [ 1 1   0 0 0 0 -1 0] } ;

A = 2.5e-4 ;
crossSecsParams = {[ A ] } ;

springsParams = { [ inf  0  inf  0  inf   0 ] ; ...
                  [ 0    0  inf  0    0   0 ]   ...
               } ;


%% Analysis parameters
% [ node nodaldof scalefactor(positive or negative) ]
controlDofs = [ 2 5 -1 ] ;

% analysis parameters
stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-8 ;
stopTolForces    = 1.0e-8  ;

targetLoadFactrNR   = 2.5e7    ; % newton
targetLoadFactrNRAL = 4e7    ; % arc length

nLoadSteps       = 60    ;
incremArcLen     = .2     ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 
%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactrNR nLoadSteps ] ; 

stabilityAnalysisBoolean = 0 ;

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
% move to onsas directory and ONSAS execution

acdir = pwd ; cd(dirOnsas);
ONSAS
cd(acdir) ;
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

