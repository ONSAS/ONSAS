%% Von Mises truss example using Newton-Raphson Arc-Length Method
%
%%
clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'staticVonMisesTrussNRAL' ;

% uncomment to delete variables and close windows
% clear all, close all

%% Structural properties
E = 210e9 ;  nu = 0 ;  rho = 0 ;
materialParams     =   cell(1,1)    ;
materialsParams{1} = [ rho 1 E nu ] ;

% each row shows the properties of each section: A, Iy Iz and J
A = 2.5e-4 ;
crossSecsParams = [ A 2 2 4 ] ;

auxx = cos(65*pi/180) * 2 ;
auxy = sin(65*pi/180) * 2 ;
imperfPerc = .0 ;

Nodes = [      0  0     0  ; ...
            auxx*(1+imperfPerc)  0  auxy  ; ...
          2*auxx  0     0  ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  0  inf  0  inf 0 ; ...
                 2    0  0  inf  0    0 0 ; ...
                 3  inf  0  inf  0  inf 0   ...
               ];

Conec = [ 1 2 0 0 1 1 1 ;
          2 3 0 0 1 1 1 ] ;

%% Loading parameters

nodalVariableLoads   = [ 2  0  0  0  0 -1  0 ];

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

stabilityAnalysisBoolean = 1 ;

% analytical solution using engineering strain
analyticSolFlag        = 2    ;
analyticCheckTolerance = 1e-4 ;
l0 = sqrt(auxx^2 + auxy^2) ;
analyticFunc = @(w) -2 * E*A* ( (  (auxy+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxy+(-w)).^2 + auxx^2) )) ) ...
 .* (auxy+(-w)) ./ ( sqrt((auxy+(-w)).^2 + auxx^2) )  ; 

%% Output parameters
printFlag = 0 ;
plotParamsVector = [ 3 10];

sectPar = [12 .1 .1]

%% ONSAS execution
% move to onsas directory and ONSAS execution

acdir = pwd ; cd(dirOnsas);
ONSAS
cd(acdir) ;

controlDispsNRAL = controlDisps ;
loadFactorsNRAL  = loadFactors ;
analyticNRAL = analyticVals ;

nLoadSteps       = 6    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNR nLoadSteps ] ; 
%~ plotParamsVector = [ 3 ];
plotParamsVector = [ 0 ];
problemName = 'staticVonMisesTrussNR' ;

acdir = pwd ; cd(dirOnsas);
ONSAS
cd(acdir) ;

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

