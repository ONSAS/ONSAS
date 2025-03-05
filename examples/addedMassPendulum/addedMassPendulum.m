% nonlinear added mass pendulum eample

close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );

% input scalar parameters
EA = 1e8  ;  nu = 0       ; 
A  = 0.1  ;  l0  = 3.0443 ;
m  = 10   ;  angle_init = 2 ; % degrees

% computed scalar parameters
E   = EA / A ;              % young modulus
d   = 2 * sqrt(A/pi);       % cross-section diameter
rho_structure = 2 * m / ( A * l0 ) ;  % solid density

% analytic solution
massratio = 1;
AMcoef =  1 + 1/massratio ;
rhoFluid = rho_structure/massratio; nuFluid = 1e-6;
g = 9.80     ;

T_analy = 2*pi * sqrt( AMcoef * 1/g * ( d^2/(8*l0) + 2*l0/(3) ) );

% materials
materials.modelName   = 'elastic-rotEngStr' ;
materials.modelParams = [ E nu ] ;
materials.density     = rho_structure ;

% md### elements
elements(1).elemType = 'node' ;

elements(2).elemType = 'frame';
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [ d ] ;

anonymus_null = @(beta,Re) 0 ;
elements(2).aeroCoefFunctions = { anonymus_null, anonymus_null, anonymus_null };
elements(2).massMatType  = 'consistent';

% md### boundaryConds
boundaryConds(1).imposDispDofs = [ 1 2 3 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;

% md### initial Conditions
initialConds = {} ;

% md### mesh parameters
% mdThe coordinates considering a mesh of two nodes is:
x_ini = sind(angle_init)*l0 ;
mesh.nodesCoords = [  0      0  l0                       ; ...
                      x_ini  0  l0*(1-cosd(angle_init))  ] ;

mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1  1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0  1 2 ] ;

% md### analysisSettings
analysisSettings.deltaT        = T_analy/100  ;
analysisSettings.finalTime     = T_analy*.5 ;
analysisSettings.methodName    = 'newmark';
analysisSettings.stopTolDeltau = 1e-8 ;
analysisSettings.stopTolForces = 1e-12 ;
analysisSettings.stopTolIts    = 30    ;

analysisSettings.booleanSelfWeight = true ;

analysisSettings.fluidProps = {rhoFluid; nuFluid; @(x,t) zeros(3,1) } ;
analysisSettings.addedMassBool = true  ;

otherParams = struct();
otherParams.problemName     = 'addedMassPedulum';
otherParams.plots_format       = 'vtk' ;
% md
% mdFirst the input structs are converted to structs with the model information
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% mdAfter that the structs are used to perform the numerical time analysis
[ matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;
% md
% md the report is generated
outputReport( modelProperties.outputDir, modelProperties.problemName )

times  = (0:size(matUs,2)-1) * analysisSettings.deltaT ;
theta_ana = angle_init * cos( 2*pi / T_analy .* times);

x_num = x_ini + matUs(7,:)' ;
theta_num = asin( x_num ./ l0) * 180/pi ;

figure
plot( times, theta_ana,'r-o')
hold on, grid on
plot( times, theta_num,'b-x')
xlabel('time (s)')
ylabel('angle (degrees)')

legend('analytic','numeric')
verifBoolean = abs( theta_ana(end) - theta_num(end) ) < ( 3e-3 * abs( theta_ana(end)) )
