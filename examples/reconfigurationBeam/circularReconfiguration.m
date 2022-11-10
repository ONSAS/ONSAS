%% Add in assembler at 338
%  global FDrag
% FDrag(timeVar) = sum(Faero(1:6:end)) ;
% # Reconfiguration cantilever beam example 
%----------------------------
close all, clear all ;
%----------------------------
% add dynamic case boolean for non test executions: 
testBool = true; 
%----------------------------
% add path
addpath( genpath( [ pwd '/../../src'] ) ); 
addpath( genpath( [ pwd ] ) ); 
% General  problem parameters
%----------------------------
% we load the given parameters:
[l, d, Izz, E, nu, rhoS, rhoF, nuF, dragCoefFunction, NR, cycd_vec, uy_vec ] = loadParametersCirc();
%
numElements = 10 ;
%
% materials
%----------------------------
% Since the example contains only one material and co-rotational strain element so then `materials` struct is:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]         ;
materials.density         = rhoS             ;
%
% elements
%----------------------------
% Two different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
% for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a circular section with $d$ diameter
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [ d ] ;% number of Gauass integration points and elemTypeAero field:
numGaussPoints = 4 ;
computeAeroTangentMatrix = true ;
elements(2).elemTypeAero   = [0 d 0 numGaussPoints computeAeroTangentMatrix ] ;
% The drag function name is:
elements(2).aeroCoefs = {dragCoefFunction; []; [] } ;
%
% boundaryConds
%----------------------------
% The elements are submitted to only one different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
% initial Conditions
%----------------------------
% any non homogeneous initial conditions is considered, then an empty struct is set:
initialConds = struct() ;
%
% analysisSettings Static
%----------------------------
analysisSettings.fluidProps = {rhoF; nuF; 'windVelCircStatic'} ;
% The geometrical non-linear effects are not considered in this case to compute the aerodynamic force. As consequence the wind load forces are computed on the reference configuration, and remains constant during the beam deformation. The field  _geometricNonLinearAero_ into  `analysisSettings` struct is then set to:
analysisSettings.geometricNonLinearAero = true;
% since this problem is static, then a N-R method is employed. The convergence of the method is accomplish with ten equal load steps (equally increasing the cycd order). The time variable for static cases is a load factor parameter that must be configured into the `windVelCircStatic.m` function. The load step represents the index in the fluid velocity for a given value of cycd:
analysisSettings.deltaT        =   1             ; % needs to be 1
analysisSettings.finalTime     =   NR            ;
analysisSettings.methodName    = 'newtonRaphson' ;
% Next the maximum number of iterations per load(time) step, the residual force and the displacements tolerances are set to: 
analysisSettings.stopTolDeltau =   0             ;
analysisSettings.stopTolForces =   1e-8          ;
analysisSettings.stopTolIts    =   50            ;
%
% otherParams
%----------------------------
otherParams.problemName = 'staticReconfigurationCircle';
otherParams.plots_format = 'vtk' ;
%
%
% meshParams
%----------------------------
%The coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))' * l / numElements  zeros(numElements+1,2) ];
%The connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
% then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1] ;
% Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md 
%md## Run ONSAS
%md---------------------
%md Declare a global variable to store drag 
%
global globalFDrag
global globalNIter
globalFDrag = zeros(analysisSettings.finalTime, 1) ;
globalNIter = zeros(analysisSettings.finalTime + 1, 1) ;
%
% Run ONSAS 
%
[matUs] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%### Numeric solution
%
% The numerical solution is extracted. First the reference coordinaes
xref = mesh.nodesCoords(:,1) ;
yref = mesh.nodesCoords(:,2) ;
zref = mesh.nodesCoords(:,3) ;
% 
%## Validation results
%---------------------
% Compute the values of R and CyCd:
numLoadSteps = size(matUs, 2) ;
timeVec = linspace(0,analysisSettings.finalTime, numLoadSteps) ;
% initialize vectors
Cy = zeros(numLoadSteps-1, 1)                 ;
R  = zeros(numLoadSteps-1, 1)                 ;
C_d = feval( elements(2).aeroCoefs{1}, 0 , 0) ;
% fill them
for windVelStep = 1:numLoadSteps - 1
    % Compute dimensionless magnitudes 
    windVel         = feval( analysisSettings.fluidProps{3,:}, 0, timeVec(windVelStep + 1 ) ) ;
    normWindVel     = norm( windVel )                                                         ;
    dirWindVel      = windVel / normWindVel                                                   ;
    Cy(windVelStep) =  1/2 * rhoF * normWindVel^2 * (l)^3 *d / (E*Izz)                        ;

    % numeric drag 
    FDragi = globalFDrag(windVelStep)                 ;
    FDRef  = 1/2 * rhoF * normWindVel^2 * C_d * d * l ;
    R(windVelStep) =  abs(FDragi)/(FDRef )            ;

end
%
%### Gosselin et.Al 2010 solution
%
% resudrag (cycd, R) and def wich contains de deformed configuration for 10^i cycyd values: 
load( 'Gosselin2010_data.mat', 'def', 'resudrag')
%md
%md### Validation plots
%md
%md The plot parameters are:
lw = 4 ; ms = 5 ;
axislw = 1  ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ;
Gline = 'k-'; ONSASline = 'bo'  ;     
folderPathFigs = './output/figs/' ;
mkdir(folderPathFigs) ;
%md The modified Cauchy number vs R is plotted:  
fig1 = figure(1) ;
hold on
loglog(C_d*Cy       , R             , ONSASline , 'linewidth', lw, 'markersize', ms ) ;
loglog(resudrag(:,1), resudrag(:,2) , Gline     , 'linewidth', lw, 'markersize', ms ) ;
% add legend
legend('ONSAS', 'Gosselin2010')
% set labels legend
labx=xlabel(' Cy* ');    laby=ylabel('R');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
% set fonts
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
grid on
% save figure
namefig1 = strcat(folderPathFigs, 'CyR.png') ;
% print(namefig1,'-dpng')
%
%md The number of iterations vs Cauchy number is then plotted:  
%
fig2 = figure(2) ;
hold on
semilogx(C_d*Cy, globalNIter(3:end), ONSASline, 'linewidth', lw, 'markersize', ms );
% add legend
legend('iters')
labx=xlabel(' Cy* ');    laby=ylabel('Number of iteration');
% set fonts
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
grid on
% save figure
namefig2 = strcat(folderPathFigs, 'CyNiter.png') ;
% print(namefig2,'-dpng')
%
%md The modified Cauchy number vs R is plotted:  
%
fig3 = figure(3) ;
hold on
plot(xref, yref  , 'k--' , 'linewidth', lw, 'markersize', ms   );
for nr = 1:NR 
  % Numerical deformed coordinates solution
  xdef = xref + matUs(1:6:end,nr+1) ;
  ydef = yref + matUs(3:6:end,nr+1) ;
  zdef = zref + matUs(5:6:end,nr+1) ;
  thetaXdef = matUs(2:6:end,nr+1)   ;
  thetaYdef = matUs(4:6:end,nr+1)   ;
  thetaZdef = matUs(6:6:end,nr+1)   ;
  % Gosselin deformed coordinates solution
  xdefG = def(1,:,nr)               ;
  ydefG = -def(2,:,nr)              ;
  % Plot 
  plot(xdef ,  ydef,  ONSASline, 'linewidth', lw, 'markersize', ms );
  plot(xdefG, ydefG,  Gline    , 'linewidth', lw, 'markersize', ms );
end
% add legend
legend('Gosselin2010', 'ONSAS')
labx=xlabel('x [m]');    laby=ylabel('y [m]');
% set fonts
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
grid on
axis equal
% save fig
namefig3 = strcat(folderPathFigs, 'xy.png') ;
% print(namefig3,'-dpng')
%
%### Verification boolean
%
% The verification boolean is computed as for the deformed configurations and the cycd curve
% deformed coordinates dif norm
vecDifDeform =  [ norm( ydef - ydefG(1:numElements*10:end)') ;...
                  norm( xdef - xdefG(1:numElements*10:end)') ] ;

% verification boolean deformed 
verifBooleanDef =  vecDifDeform <=  2e-2 * l ;
% cycd vs R verification boolean is: 
verifBooleanR = abs(R(end) - resudrag(end,2) ) <  5e-3 ;
% The example verifboolean is:
verifBoolean = verifBooleanR && all(verifBooleanDef)
if ~testBool
  %
  % Dynamic Case
  %----------------------------
  %
  % since this case is highly dynamic the tangent matrix of the aerodynamic force 
  % vector is not necessary 
  %
  % elements
  %----------------------------
  elements(2).massMatType = 'consistent'                                        ;
  computeAeroTangentMatrix = false                                              ;
  elements(2).elemTypeAero   = [0 d 0 numGaussPoints computeAeroTangentMatrix ] ;
  %
  % analysisSettings Dynamic
  %----------------------------
  analysisSettings.fluidProps = {rhoF; nuF; 'windVelCircDynamic'}             ;
  analysisSettings.finalTime     =   4                                        ;
  analysisSettings.deltaT        =   0.0025                                   ;
                                       
  numTimeSteps  = round(analysisSettings.finalTime / analysisSettings.deltaT) ;
  analysisSettings.methodName    = 'newmark'                                  ;
  analysisSettings.stopTolDeltau =   1e-11                                    ;
  analysisSettings.stopTolForces =   1e-6                                     ;
  analysisSettings.stopTolIts    =   15                                       ;
  %
  % otherParams
  %----------------------------
  otherParams.problemName = 'dynamicReconfigurationCircle';
  %
  % Run ONSAS 
  %
  [matUs] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

end