%md# Static Von-Mises Truss example

close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 10 ; L = 2 ;
%E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ;
Kplas = E*0 ;
sigma_Y_0 = 250e6 ;

% x and z coordinates of node 2
x2 = cos( ang1*pi/180 ) * L ;
z2 = sin( ang1*pi/180 ) * L ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(A*4/pi) } ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs =   3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 3.0e8*t     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ] ;

initialConds                = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0 0  2 3 ] ;

analysisSettings.methodName    = 'arcLength' ;

analysisSettings.deltaT        =   0.1  ;
analysisSettings.finalTime     =   2    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.incremArcLen = [ 0.003*ones(1, 4) 0.002*ones(1,20) ]                           ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'static_plastic_von_mises_truss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


deltas = -matUs(6+5,:)' ;
eles = sqrt( x2^2 + (z2-deltas).^2 ) ;

valsLin = -2*E*A * (z2-deltas ) ./ eles .* (eles - L) ./ L ;

valsFlu = 2*sigma_Y_0 * A * ( z2 - deltas ) ./ eles  ;

figure
plot( deltas , valsLin, 'b-x' )
grid on, hold on
plot( deltas , valsFlu, 'r-o' )
plot( deltas , loadFactorsMat(:,2), 'g-s' )
