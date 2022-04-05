%md# Beam truss pendulum problem
close all, clear all
%mdProblem name:
otherParams.problemName = 'beamTrussPendulum' ;
addpath( genpath( [ pwd '/../../src'] ) );
%mdTuss element material and geometrical scalar properties
Et = 200e9 ; nut = 0.3 ;  rhot = 8e3 ;
dt = .02 ; At = pi*dt^2/4 ; lt = 1 ;
%
%mdFrame element geometrical properties :
Ef = Et/20 ; nuf = 0.3;  rhof = rhot;
df = dt*2; Ab = pi*df^2/4;  lf = lt; If = pi*df^4/64 ;
%
%mdScalar parameters of the mesh:
numElemF = 10; numElemT = 1;
%
%md### MEBI parameters
%md
%md### materials
%mdtruss material:
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ Et nut ] ;
materials(1).density = rhot ;
%mdframe material:
materials(2).hyperElasModel  = '1DrotEngStrain' ;
materials(2).hyperElasParams = [ Ef nuf ] ;
materials(2).density = rhof ;
%md
%md### elements
%md
%mdNodes:
elements(1).elemType = 'node'           ;
%md
%mdTruss:
elements(2).elemType = 'truss'                    ;
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ dt dt ]   ;
elements(2).massMatType = 'consistent'            ;
%md
%mdFrame:
elements(3).elemType = 'frame'                  ;
elements(3).elemCrossSecParams{1,1} = 'circle'  ;
elements(3).elemCrossSecParams{2,1} =  df       ;
elements(3).massMatType = 'consistent'       ;
%md
%md### boundaryConds
%md
%mdFixed node:
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdLoaded node:
massPendulum = 5; g = 9.8;
boundaryConds(2).loadsCoordSys = 'global'                       ;
boundaryConds(2).loadsTimeFact = @(t) 1.0                       ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -massPendulum*g 0 ]  ;
%md
%md### initial Conditions
%md
initialConds                = struct() ;
%md### mesh parameters
%md
%mdnodesCoords:
mesh.nodesCoords = [  ( 0:(lf/numElemF):lf)'        zeros(numElemF + 1, 1)  zeros(numElemF + 1,1)  ; ...
		                  ( lf+(1:(lt/numElemT):lt)' )  zeros(numElemT    , 1)  zeros(numElemT,1)      ] ;
%md
%mdConec matrix:
%mdnodes conectivity:
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1                   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElemF+1+numElemT ] ;
%
%mdelements conectivity:
auxConecElem  = [ %MEBI frame elements
                  [ (ones(numElemF,1)*2 )    (ones(numElemF,1)*3)      (zeros(numElemF,1))    (zeros(numElemF,1)) ...
                    %ElemNodes...
                    (1:(numElemF))'                         (2:numElemF+1)' ];
                    %MEBI truss elements
                  [ (ones(numElemT,1)*1 )    (ones(numElemT,1)*2)      (zeros(numElemT,1))     (zeros(numElemT,1)) ...
                    %ElemNodes...
                    (numElemF + 1 : numElemF + numElemT)'	 	(numElemF + 2 : numElemF + numElemT + 1)'  ]
                ] ;
%md
%mdFill conec cell with auxConecElem:
for i =  1:numElemF + numElemT
    mesh.conecCell{2 + i,1} = auxConecElem(i,:);
end
%md
%md### analysisSettings
%md
analysisSettings.deltaT        =   0.025 ;
analysisSettings.finalTime      =   1.0 ;
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   30   ;
%md
%md
%md### otheParams
%md
otherParams.plotsFormat = 'vtk' ;

%mdRun ONSAS
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%mdExtract numerical results:
%mdTime vector:
timeVec = linspace(0, analysisSettings.finalTime, ceil(analysisSettings.finalTime / analysisSettings.deltaT) + 1 )' ;
%md
%mdDisp Z at frame end node
controlDofF  = (numElemF + 1)*6 - 1    ;
dispZF       = matUs(controlDofF,:)' ;
%mdDisp Z at truss end node
controlDofT  = (numElemF + numElemT + 1)*6 - 1 ;
dispZT       = matUs(controlDofT,:)' ;
%md
%md### Plots
%mdPlot parameters:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
%md
%mdOutput diplacments vs time are plotted:
figure
plot( timeVec, dispZF ,'r-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time (s)');   laby = ylabel('uF_z (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% print('output/rightAngleCantilever_uA_y.png','-dpng')

figure
plot( timeVec, dispZT ,'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time (s)');   laby = ylabel('uT_z (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% print('output/rightAngleCantilever_uA_z.png','-dpng')
