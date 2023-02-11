%md# Cantilever Modal Analysis
%md
%mdBefore defining the structs, the workspace is cleaned and the ONSAS directory is added to the path
close all, clear all, addpath( genpath( [ pwd '/../../src'] ) );
%md
%mdThe material scalar parameters are set.
E = 200e9 ; nu = 0.3;  rho = 700;
%mdThe cross-section of the beam is rectangular. The widths and other geometry scalar parameters are computed.
l = 10 ; diam = .01 ; 
numElements = 16 ; % Number of elements

tf     = 0.3     ; % s
deltat = 0.1   ; % s

global exportFirstMatrices
exportFirstMatrices = true;

global modalAnalysisBoolean
modalAnalysisBoolean = true;

%md### materials
materials.hyperElasModel  = '1DrotEngStrain' ;%'linearElastic';%
materials.hyperElasParams = [ E nu ]         ;
materials.density         = rho              ;
%md
%md### elements
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%elements(2).elemCrossSecParams = { 'rectangle' , [ty tz] } ;
elements(2).elemCrossSecParams = { 'circle' , diam } ;
%md The consistent mass approach is considered for the dynamic analysis
elements(2).massMatType = 'consistent';
%md
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%md and the second corresponds to a time dependant external force
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) t ;
boundaryConds(2).loadsBaseVals = [ 0 0 1 0 0 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;

mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1 ] ;

for i=1:numElements
  mesh.conecCell{ i+2, 1 } = [ 1 2 0 0  i i+1 ] ;
end

%md### analysisSettings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   deltat  ;
analysisSettings.finalTime     =   tf   ;
analysisSettings.stopTolDeltau =   1e-10 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   20   ;
%md
%md## otherParams
otherParams.problemName = 'cantilever_modal_analysis';
%md ONSAS execution
[coRotMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md

%%
addpath('./output');
filename = './output/matrices.mat';
load(filename);

KTred = KT( neumdofs, neumdofs ) ;
Mred  = massMat( neumdofs, neumdofs ) ;

%%
s = linspace(0,l,numElements+1)' ;
z = linspace(0,l,numElements*10+1)' ;
[a, b] = eig( full(KTred), full(Mred) ) ;
eigenvalues = flipud( diag(b) ) ;

numer_modes = [ zeros(6,(numElements)*6); fliplr(a) ] ;

% Analytical solution 
betavect = [1.8751 4.6941 7.8547 10.9955 14.13717];
sigvect  = [0.73409 1.01846 0.9992 1.00003 1     ]; 
modes = [1 2 3 4 5 ] ;
modes_errors = zeros(size(modes));
close all
for mode = modes 
  figure(mode)
  num_disp_mode = numer_modes(5:6:end, mode ) ;
  num_disp_mode = num_disp_mode ./ max(abs( num_disp_mode )) * sign( num_disp_mode(end) ) ;
  plot( s, num_disp_mode, 'rx', 'linewidth', 1.2,'markersize', 12)
  hold on, grid on, title(['mode ' num2str(mode) ])
  %freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma)/l^4)/(2*pi);
  S1 = cosh(betavect(mode)*z/l) - cos(betavect(mode)*z/l) - sigvect(mode)*sinh(betavect(mode)*z/l) + sigvect(mode)*sin(betavect(mode)*z/l) ;
  S1 = S1 ./ max(abs(S1)) * sign( S1(end) ) ;
  spanplot = 10;
  S1_points=S1(1:spanplot:end) ;
  modes_errors( mode ) = norm( S1_points - num_disp_mode ) / norm( S1_points ) ;
  plot( z(1:spanplot:end), S1_points , 'bo' , 'linewidth', 1.0,'markersize',11)
  plot( z, S1, 'b' , 'linewidth', 1.0,'markersize',11)
  legend('ONSAS','analytical')
end

modes_errors

verifBoolean = ( modes_errors(1) < 1e-3 )  &&  ( modes_errors(2) < 1e-3 )
