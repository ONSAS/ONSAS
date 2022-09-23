%md# Static Von-Mises Truss example

close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ;
Kplas = E*.1 ;
sigma_Y_0 = 25e6 ;

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

analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   2e-5  ;
analysisSettings.finalTime     =   1e-3    ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'static_plastic_von_mises_truss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = -matUs(6+5,:)' ;
eles = sqrt( x2^2 + (z2-deltas).^2 ) ;

valsLin = -2* ( E * (eles - L) ./ L ) * A .* (z2-deltas ) ./ eles  ;

Etan = E*Kplas / ( E + Kplas) ;
sigmas_hard = sigma_Y_0 + Etan * ( abs( (eles - L)) ./ L - sigma_Y_0/E ) ;

valsFlu = 2*( sigmas_hard * A ) .* ( ( z2 - deltas ) ./ eles )  ;
valsP = min( valsLin, valsFlu) ;


% softening
Kplas = -E*.05 ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

analysisSettings.methodName    = 'arcLength' ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.incremArcLen = [ 2e-4 4e-5*ones(1,100)];

[matUsB, loadFactorsMatB ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
deltasB = -matUsB(6+5,:)' ;

eles = sqrt( x2^2 + (z2-deltasB).^2 ) ;

valsLin = -2*E*A * (z2-deltasB ) ./ eles .* (eles - L) ./ L ;

Etan = E*Kplas / ( E + Kplas) ;
sigmas_hard = sigma_Y_0 + Etan * ( abs( (eles - L)) ./ L - sigma_Y_0/E ) ;

valsFlu = 2*( sigmas_hard * A ) .* ( ( z2 - deltasB ) ./ eles )  ;
valsPB = min( valsLin, valsFlu) ;

figure
plot( deltas , valsP, 'b-x' )
grid on, hold on
plot( deltas , loadFactorsMat(:,2), 'r-o' )
plot( deltasB , loadFactorsMatB(:,2), 'g-o' )
plot( deltasB , valsPB, 'k-*' )
legend('analytic-hard','numeric-hard', 'analytic-soft','numeric-soft')