%% Von Mises truss example using Newton-Raphson Arc-Length Method
%
%%

clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'staticVonMisesTruss' ;

% uncomment to delete variables and close windows
% clear all, close all

%% Structural properties

Es = 210e9 ;
nu = 0 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu] ;

% each row shows the properties of each section: A, Iy Iz and J
A = 2.5e-4 ;
secGeomProps = [ A 2 2 4 ] ;

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
controlDofInfo = [ 2 5 -1 ] ;
%~ controlDofInfo = [ 2 1 +1 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ;  dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ; 

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-10 ;
%~ stopTolForces    = 1.0e-6  ;
stopTolForces    = 1.0e-10  ;

targetLoadFactrNR   = 2e7    ; % newton
targetLoadFactrNRAL = 4e7    ; % arc length

nLoadSteps       = 60    ;
incremArcLen     = .1     ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNRAL nLoadSteps incremArcLen ] ; 
%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactrNR nLoadSteps ] ; 

stabilityAnalysisBoolean = 1 ;

% analytical solution using engineering strain
analyticSolFlag        = 2    ;
analyticCheckTolerance = 1e-4 ;
l0 = sqrt(auxx^2 + auxy^2) ;
analyticFunc = @(w) -2 * Es*A* ( (  (auxy+(-w)).^2 + auxx^2 - l0^2 ) ./ (l0 * ( l0 + sqrt((auxy+(-w)).^2 + auxx^2) )) ) ...
 .* (auxy+(-w)) ./ ( sqrt((auxy+(-w)).^2 + auxx^2) )  ; 

%% Output parameters
printflag = 0 ;
plotParamsVector = [ 1 ];

%% ONSAS execution
% move to onsas directory and ONSAS execution

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;

compar = figure ;
plot( controlDisps, loadFactors,'b-x' )
hold on, grid on

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactrNR nLoadSteps ] ; 

plotParamsVector = [ 0 ];

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;


plot( controlDisps, loadFactors,'r-s' )
