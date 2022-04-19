%md# Dynamic linear cantilever beam example
%md---
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
E = 210e9 ;  nu = 0.3 ; rho = 2100 ;
% geometrical scalar parameters
l = 5 ; ty = .3 ;  tz = .1 ; 
% the number of elements of the mesh
numElements = 10 ;
% reference bending periods in y and z direction 
Iyy = ty * tz^3 / 12 ;
Tbyy = 2 * pi / sqrt (4.73 * E * Iyy / (rho * (tz * ty *l ) ) ) ;
Izz = ty^3 * tz / 12 ;
Tbzz = 2 * pi / sqrt (4.73 * E * Izz / (rho * (tz * ty *l ) ) ) ;
%md
%md### materials
%md
materials.hyperElasParams = [ E nu ] ;
materials.density = rho ;
%md
%md### elements
%md
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ty tz]     ;
elements(2).massMatType             = 'consistent';
%md
%md### boundaryConds
%md
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdand the second corresponds to an incremental nodal moment, where the target load produces a circular form of the deformed beam.
boundaryConds(2).loadsCoordSys = 'global'        ;
timeZeroMoment = Tbyy / 5 ;
boundaryConds(2).loadsTimeFact = @(t) E * Iyy * 2 * pi / l * t * (t <= timeZeroMoment) ;
%md
%md### initialConds
%md
initialConds                = struct() ;
%md
%md### mesh parameters
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1 ] ;
for i=1:numElements,
  mesh.conecCell{ i+2,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
%md
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   Tbyy/100 ;
analysisSettings.finalTime      =  Tbyy  ;
analysisSettings.stopTolDeltau =   1e-12 ;
analysisSettings.stopTolForces =   1e-12 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
%md
otherParams.plotsFormat = ' ' ;
%md
%md## Case 1: Linear beam element submitted to a bending moment in $y$ 
%md
%md## materials
%md
materials.hyperElasModel  = 'linearElastic' ;
%md
%md## boundaryConds
%md 
boundaryConds(2).loadsBaseVals = [ 0 0 0 1 0 0 ] ;
%md
%md Others params
%md
otherParams.problemName = 'linearCantileverMomentY';
%md
%md execute ONSAS
%md
[matUsYLinear, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Case 2: Linear beam element submitted to a bending moment in $z$ 
%md
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 0 1 ] ;
%md
%md Others params
%md
otherParams.problemName = 'linearCantileverMomentZ';
%md
%md execute ONSAS
%md
[matUsZLinear, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Case 3: Co-rotational beam element submitted to a bending moment in $z$ 
%md
materials.hyperElasModel  = '1DrotEngStrain' ;
%md
%md
%md Others params
%md
otherParams.problemName = 'coRotCantileverMomentZ';
%md
%md execute ONSAS
%md
[matUsZCorot, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Case 4: Co-rotational beam element submitted to a bending moment in $y$ 
%md
%md## boundaryConds
%md 
boundaryConds(2).loadsBaseVals = [ 0 0 0 1 0 0 ] ;
%md
%md execute ONSAS
%md
%md Others params
%md
otherParams.problemName = 'coRotCantileverMomentY';
%md
[matUsYCorot, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md## verification
%mdFor all cases the displacements solution of co-rotational and linear beam element must be the same.
% dofs to plot at the end node
dofZendNode = 6*(numElements +1) - 1 ;
dofYendNode = 6*(numElements +1) - 3 ;
% Solution case 1:
uZcase1 = matUsYLinear(dofZendNode, :);
% Solution case 2:
uYcase2 = matUsZLinear(dofYendNode, :);
% Solution case 3:
uYcase3 = matUsZCorot(dofYendNode, :) ;
% Solution case 4:
uZcase4 = matUsYCorot(dofZendNode, :) ;
% find the maximum displacement magnitude of node A
maxUz = max( max( abs(uZcase1) ), max( abs( uZcase4 ) ) ) ;
maxUy = max( max( abs(uYcase2) ), max( abs( uYcase3 ) ) ) ;
%md
%md set the tolerance and define the boolean: 
tolVerifDisp = 1e-2 ;
verifBooleanZ =  ( norm( uZcase1(end) - uZcase4(end) ) <  tolVerifDisp * maxUz ) ;
verifBooleanY =  ( norm( uYcase2(end) - uYcase3(end) ) <  tolVerifDisp * maxUy ) ;
%md all cases must be verifyed, so then:
verifBoolean    = verifBooleanZ && verifBooleanY ;
%md
%md### Plots
%mdPlot parameters:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
%md
figure, hold on, grid on
%md time vector
timeVec = linspace( 0,analysisSettings.finalTime, size(matUsYCorot,2) );

% plot linear cases
plot(timeVec, matUsYLinear(dofZendNode, :),'r-x' , 'linewidth', lw,'markersize',ms )
plot(timeVec, matUsZLinear(dofYendNode, :),'b-x' , 'linewidth', lw,'markersize',ms )
% plot co-rotational cases
plot(timeVec, matUsYCorot(dofZendNode, :),'k-' , 'linewidth', lw,'markersize',ms )
plot(timeVec, matUsZCorot(dofYendNode, :),'c-' , 'linewidth', lw,'markersize',ms )
%plot time zero moment
plot([timeZeroMoment timeZeroMoment], [ min(  matUsYCorot(dofZendNode, :) ) max( matUsZLinear(dofYendNode, :) ) ],...
     'g:','linewidth', lw+5, 'markersize', ms )
% legends
legend('linearMomentZ','linearMomentY', 'co-rotMomentZ',  'co-rotMomentY', 'time moment 0 N.m',...
       'location', 'southeast')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
axis([ 0 1.5*analysisSettings.finalTime min(  matUsYCorot(dofZendNode, :) ) max( matUsZLinear(dofYendNode, :) ) ])
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;