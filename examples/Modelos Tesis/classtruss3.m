%md# Static Von-Mises Truss example
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ;
Kplas = E*.1 ;
sigma_Y_0 = 25e6 ;

% x and z coordinates of node 2
x2 = cos( ang1*pi/180 ) * L ;
z2 = sin( ang1*pi/180 ) * L ;

materials = struct();
materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements = struct();
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(A*4/pi) } ;

boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs =   3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 3.0e8*t     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ] ;

mesh = struct();
mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1  3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2  2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0  2 3 ] ;

initialConds                = struct() ;

analysisSettings = struct();
analysisSettings.deltaT        =   2e-5  ;
analysisSettings.finalTime     =   1e-3    ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

% an Eternal Golden Braid

otherParams.problemName       = 'staticVonMisesTruss_NRAL_Jirasek_Green' ;
analysisSettings.methodName   = 'arcLength'                      ;
analysisSettings.finalTime    = 1                               ;
analysisSettings.incremArcLen = [0.15]                     ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 11 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRAL_Jirasek_Green =  -matUs(11,:) ;
loadFactorsNRAL_Jirasek_Green  =  loadFactorsMat(:,2) ;
analyticLoadFactorsNRAL_Jirasek_Green = analyticLoadFactorsGreen(controlDispsNRAL_Jirasek_Green);
difLoadGreenNRAL_Jirasek = analyticLoadFactorsNRAL_Jirasek_Green' - loadFactorsNRAL_Jirasek_Green ;

verifBoolean =  ( ( norm( difLoadEngRot            ) / norm( loadFactorsNREngRot           ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreen             ) / norm( loadFactorsNRGreen            ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreenNRAL         ) / norm( loadFactorsNRALGreen          ) ) <  1e-4 ) ...
             && ( ( norm( difLoadGreenNRAL_Jirasek ) / norm( loadFactorsNRAL_Jirasek_Green ) ) <  1e-4 ) ;

lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure
plot( controlDispsNREngRot, analyticLoadFactorsNREngRot( controlDispsNREngRot) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNREngRot, loadFactorsNREngRot, 'k-o' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, analyticLoadFactorsGreen( controlDispsNRALGreen ), 'g-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRGreen, loadFactorsNRGreen, 'r-s' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRALGreen, loadFactorsNRALGreen, 'c-^' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRAL_Jirasek_Green, loadFactorsNRAL_Jirasek_Green, 'y-*' , 'linewidth', lw,'markersize',ms )
plot( controlDispsNRLinearElastic, loadFactorsNRLinearElastic, 'm-+' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement w(t)');   laby = ylabel('\lambda(t)') ;
legend( 'analytic-RotEng', 'NR-RotEng','analytic-Green', 'NR-Green','NRAL-Green', 'NRAL-Jirasek-Green','LinearElastic', 'location','northoutside')
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/vonMisesTrussCheck.png','-dpng')



