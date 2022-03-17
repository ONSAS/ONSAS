%md# Aerodynamic linear static cantilever beam example
close all, clear all ;
% add path
accDir = pwd ;
addpath( genpath( [ accDir '/../../src'] ) );
% measure execution time
tic
% material scalar parameters
E = 3e7 ;  nu = 0.3 ; rho = 700 ; G = E / (2 * (1+nu)) ;
% geometrical scalar parameters
l = 2 ; d = l/10 ; J = pi * d ^ 4 / 64 ; Iyy = J / 2 ; Izz = J / 2 ;  
% the number of elements of the mesh
numElements = 5 ;
%md##Numerical solution
%md### MEBI parameters
%md
%md### materials
%md Since the example contains only one rod the fields of the `materials` struct will have only one entry. In order to validate against and analytical solution solved with DiffEq.jl file the Euler-Bernoulli linear element is employed:
materials.hyperElasModel  = 'linearElastic' ;
%mdthis moddel is valid only under small displacements regime.
%mdalthough, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the beam elements will be set as a circlular-cross section with $d$ diamater, then the first entery of the elemCrossSecParams cell is filled with 'circular' and then the second with a the diameter length 'd':
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = d          ;
%md now the element aerodynmamic properties are defined
%md The number of Guass integration points selected is
numGaussPoints  = 4 ;
%the fluid density of the fluid is conisdered as air at 20 degreess and atmospheric pressure
formulCase = 4 ;
%mdthe element typeAero in the first 3 entries contain the direction and the norm of the chord vector. Since a cylinder is simulated for this case and the orientation is alog y local axis, then:
elements(2).elemTypeAero   = [0 d 0 numGaussPoints formulCase ] ;
%mdthe drag function dependeing on the angle of incidence is
elements(2).userDragCoef   = 'dragCoefNonLinear'   ;
%md no lift and moment is condiered so
% elements(2).userLiftCoef   = 'liftCoefNonLinear'   ;
% elements(2).userMomentCoef = 'momentCoefNonLinear' ;
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first and unique BC corresponds to a welded condition for a cantilever beam (all 6 dofs set to zero)
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%md the following beam elements only differ in the second and first enerty (MEBI), the first material is considered, then the second enerty refers to the beam aerodynamic element and finally no initial and boundary condition is used, as consequence
for i=1:numElements
  conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md A first static small displacements case is executed 
%md-------------------------------------
%md## Static small displacements(SD) case
%md -------------------------------------
%md### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   .1            ;
analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
analysisSettings.stopTolDeltau =   0             ;
analysisSettings.stopTolForces =   1e-9          ;
analysisSettings.stopTolIts    =   50            ;
%md the name of the wind velocity function is: 
analysisSettings.userWindVel   = 'windVelNonLinearStaticSD';
%md geometrical nonlinearity in the wind force is taken into account in this example, as consequence the drag force is computed in the rigid configuration:
%mdloadsSettings
analysisSettings.booleanSelfWeight      = false ;
analysisSettings.geometricNonLinearAero = true  ;
%md
%md## otherParams
otherParams.problemName = strcat( 'onsasExample_nonLinearCantilever_staitcSD_formulation=', num2str(formulCase) ) ;
otherParams.controlDofs = [ numElements+1  4 ] ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUsSD, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Assamble static Julia solution
%mdread julia solution
solCell     = dlmread ('solJDiffEq.csv', ',', 1, 0) ;
xJulia      = solCell(:,1) ;
thetaZJulia = solCell(:,4) ;
uYJulia     = solCell(:,5) ;
%md fill julia solution vector
dSolJulia          = zeros( 6*size( xJulia, 1 ), 1 ) ;
dSolJulia(3:6:end) = -1 * uYJulia'                   ;
dSolJulia(6:6:end) = -1 * thetaZJulia'               ;

%md create and fill solution with the same size of the vector lacating their cooridnates
dSol = zeros( (numElements + 1) * 6 ,1 );
numElemJulia = size( dSolJulia(1:6:end) ) - 1;

for elem = 1:size(conecElemMatrix,1)
  % Localizate dofs and element coordinates
  %Cords element
  x = mesh.nodesCoords(conecElemMatrix(elem,5:6)',1) ;
  y = mesh.nodesCoords(conecElemMatrix(elem,5:6)',2) ;	
  z = mesh.nodesCoords(conecElemMatrix(elem,5:6)',3) ;

  elemCoords = [x(1);y(1);z(1);x(2);y(2);z(2)]; % def segun la funcion fint_e

  %Dofs element
  dofElemVec        = 6*conecElemMatrix( elem,5:6 )-5                             ;
  dofElem           = [ dofElemVec(1) : dofElemVec(2) + 5 ]'                      ; % ocurre sii son consecutivos
  xelemCoords       = elemCoords(1:3:end)                                         ;
  solJuliaIndex     = round (xelemCoords * numElemJulia / l)                      ;
  dofNode1ElemJulia = ( solJuliaIndex(1)*6 + 1 : 1:( solJuliaIndex(1)+1 )*6)'     ;
  dofNode2ElemJulia = ( solJuliaIndex(2)*6 + 1 : 1 :( solJuliaIndex(2) + 1 )*6 )' ;
  dofElemJulia      = [ dofNode1ElemJulia; dofNode2ElemJulia ]                    ;
  % Read displacments from julia solution
  UeSol             = dSolJulia( dofElemJulia ) ;
  dSol(dofElem)     = UeSol ;
end

%md
%md## Evaluate analytical solutions
rhoAire = 1.2;
%evaluate drag/lift and moment coefficents
betaRel = acos(dot(elements(2).elemTypeAero(1:3) , [0 0 1] )) ;


if isfield(elements(2), 'userDragCoef')
  c_d = feval(elements(2).userDragCoef, betaRel) ;
else
  c_d = 0;
end
if isfield(elements(2), 'userLiftCoef')
  c_l = feval(elements(2).userLiftCoef, betaRel) ;
else
  c_l = 0;
end
if isfield(elements(2), 'userMomentCoef')
  c_m = feval(elements(2).userMomentCoef, betaRel) ;
else
  c_m = 0;
end

%mdget wind velocity
windVel = feval(analysisSettings.userWindVel, mesh.nodesCoords(1,1), analysisSettings.finalTime) ;
%mdcaracteristicDimension
dimCaracteristic = norm(elements(2).elemTypeAero(1:3)) ;

%reference coordinates
xref = mesh.nodesCoords(:,1) ;
yref = mesh.nodesCoords(:,2) ;
zref = mesh.nodesCoords(:,3) ;

%Analytic x vector
sizeAnalyticX = 100 ;
xanal = linspace( 0, l , sizeAnalyticX )' ;


%linear disp julia Diff eq solution
xdefJulia = xJulia +  dSolJulia(1:6:end) ;
ydefJulia = dSolJulia(3:6:end)           ;
zdefJulia = dSolJulia(5:6:end)           ;
%angular disp julia Diff eq solution
thetaXdefJulia = dSolJulia(2:6:end) ;
thetaYdefJulia = dSolJulia(4:6:end) ;
thetaZdefJulia = dSolJulia(6:6:end) ;

% Load numerical solution
%linear disp
xdefNum = mesh.nodesCoords(:,1) + matUsSD(1:6:end, end) ;
ydefNum = mesh.nodesCoords(:,2) + matUsSD(3:6:end, end) ;
zdefNum = mesh.nodesCoords(:,2) + matUsSD(5:6:end, end) ;
%angular disp
thetaXdefNum = matUsSD(2:6:end, end) ;
thetaYdefNum = matUsSD(4:6:end, end) ;
thetaZdefNum = matUsSD(6:6:end, end) ;

% % Plot parameters:
lw = 5 ; ms = 5 ;
% labels parameters:
labelTitle = [' Validating solution with ' num2str(numElements) ' elements' ] ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ;    
% folder to save plots
folderSDpath = './output/SD/' ;
mkdir(folderSDpath) ;

% Plot linear displacements
fig1 = figure(1) ;
hold on  
grid on
plot(xdefNum    , ydefNum,        'bo' , 'linewidth', lw, 'markersize', ms+5   );
plot(xdefJulia  , ydefJulia,      'b-' , 'linewidth', lw, 'markersize', ms     );
legend('y numeric', 'y semi-analytic', 'location', 'northEast')
labx=xlabel(' x[m] ');    laby=ylabel('y ');
% title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig1 = strcat(folderSDpath, 'linDispSD.png') ;
print(fig1, namefig1,'-dpng') ;
% close(1)
% Plot angular displacements
fig2 = figure(2) ;
hold on  
grid on
plot(xdefNum,     rad2deg(thetaZdefNum),      'bo' , 'linewidth', lw,'markersize', ms+5   );
plot(xdefJulia,   rad2deg(thetaZdefJulia),    'b-' , 'linewidth', lw, 'markersize', ms    );
legend('\theta_z numeric SD', '\theta_z semi-analyitc SD', 'location','northEast')
labx=xlabel(' x[m] ');    laby=ylabel('Angle');
% title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig2 = strcat(folderSDpath, 'angDispSD.png') ;
print(fig2, namefig2,'-dpng')
% close(2)
% %Plot 3D deformed
fig3 = figure(3) ;
hold on
plot(xref,      yref, 'k-' , 'linewidth', lw+300,'markersize', ms+200);
plot(xdefJulia, ydefJulia, 'b-' , 'linewidth', lw,'markersize', ms);
plot(xdefNum,   ydefNum, 'bo' , 'linewidth', lw,'markersize', ms+5);
legend('Reference config', 'DiffEq.jl defomred config SD',  'Numerical defomred config SD', 'location','northEast')
labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
grid on
namefig3 = strcat(folderSDpath, 'defSD.png') ;
print( fig3, namefig3, '-dpng' ) ;
%md-------------------------------------
%md## Static 2D large displacements case
%md -------------------------------------
numElementsLD2DStatic = 5 ;
%md## mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElementsLD2DStatic))'*l/numElementsLD2DStatic  zeros(numElementsLD2DStatic+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
% mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElementsLD2DStatic+1   ] ;
%md the following case only differs in the boundary condition and the node number
for i=1:numElementsLD2DStatic
  conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
materials.hyperElasModel  = '1DrotEngStrain' ;
%
% elements
%md element elemTypeAero if it is conssistent is equal to 1 and lumped 0:
elements(2).elemTypeParams = 1 ;
%md### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   .1            ;
analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
analysisSettings.stopTolDeltau =   1e-6          ;
analysisSettings.stopTolForces =   1e-6          ;
analysisSettings.stopTolIts    =   35            ;
%md the name of the wind velocity function is: 
analysisSettings.userWindVel   = 'windVelNonLinearStaticLD';
%md geometrical nonlinearity in the wind force is not taken into account in this example:
%mdloadsSettings
analysisSettings.booleanSelfWeight      = false ;
analysisSettings.geometricNonLinearAero = true  ;
%md
%md## otherParams
otherParams.problemName = strcat( 'onsasExample_nonLinearCantilever_staitcLD_formulation=', num2str(formulCase) ) ;
otherParams.controlDofs = [ numElementsLD2DStatic + 1  4 ] ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUsLD2DStatic, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%mdget wind velocity
windVel = feval(analysisSettings.userWindVel, mesh.nodesCoords(1,1), analysisSettings.finalTime) ;
% Load numerical solution
%linear disp
xdefNum = mesh.nodesCoords(:,1) + matUsLD2DStatic(1:6:end, end) ;
ydefNum = mesh.nodesCoords(:,2) + matUsLD2DStatic(3:6:end, end) ;
zdefNum = mesh.nodesCoords(:,2) + matUsLD2DStatic(5:6:end, end) ;
%angular disp
thetaXdefNum = matUsLD2DStatic(2:6:end, end) ;
thetaYdefNum = matUsLD2DStatic(4:6:end, end) ;
thetaZdefNum = matUsLD2DStatic(6:6:end, end) ;

% Folder path
folderLD2Dpath = strcat('./output/', 'LD/', '2D/') ;
mkdir(folderLD2Dpath) ;
folderLD2DStaticPath = strcat(folderLD2Dpath, 'static/') ;
mkdir(folderLD2DStaticPath) ;
% Plot linear displacements
fig4 = figure(4) ;
hold on  
grid on
plot(xdefNum,   ydefNum,                   'b-o' , 'linewidth', lw,'markersize', ms+5 );
legend('y numeric LD', 'location', 'north')
labx=xlabel(' x [m] ');    laby=ylabel('y [m] ');
% title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig4 = strcat(folderLD2DStaticPath, 'linDispLD.png') ;
print( fig4, namefig4, '-dpng' ) ;
% Plot angular displacements
fig5 = figure(5) ;
hold on  
grid on
plot(xdefNum,     thetaZdefNum,               'b-o' , 'linewidth', lw,'markersize', ms+5 );
legend('\theta z numeric LD', 'location','north')
labx=xlabel(' x [m] ');    laby=ylabel('\theta z [rad]');
% title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig5 = strcat(folderLD2Dpath, 'static/', 'angDispLD.png');
print( fig5, namefig5, '-dpng' ) ;  
%Plot 3D deformed
fig6 = figure(6) ;
hold on
grid on
plot(xref,      yref, 'k-' , 'linewidth', lw+300,'markersize', ms+200);
plot(xdefNum,   ydefNum, 'b-o' , 'linewidth', lw,'markersize', ms+5);
legend('Reference config', 'Numerical defomred config LD', 'location','northwest')
labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
namefig6 = strcat(folderLD2Dpath, 'static/', 'defLD.png') ;
print( fig6, namefig6, '-dpng' ) ;  
%md-------------------------------------
%md Dynamic 2D Case 
%md-------------------------------------
numElements2DLDDynamic = 10 ;
%md## mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements2DLDDynamic))'*l/numElements2DLDDynamic  zeros(numElements2DLDDynamic+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
% mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements2DLDDynamic+1   ] ;
%md the following case only differs in the boundary condition and the node number
for i=1:numElements2DLDDynamic
  conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end

% MEBI parameters
%
% analysisSettings
analysisSettings.finalTime   =   5                               ;
analysisSettings.deltaT      =   analysisSettings.finalTime /400 ;
analysisSettings.methodName  = 'newmark'                         ;
analysisSettings.alphaNM     =   0.25                            ;
analysisSettings.deltaNM     =   0.5                             ;
analysisSettings.userWindVel = 'windVelNonLinearDynamic'         ;
analysisSettings.stopTolIts  =   10                              ;

  % create folder 
folderLD2DDynamicPath = strcat(folderLD2Dpath, 'dynamic/') ;
mkdir(folderLD2DDynamicPath) ;

%md 
%
% otherParams
%
otherParams.problemName      = strcat( 'onsasExample_nonLinearCantilever_dynamic', '_formulation=', num2str(formulCase) ) ;
otherParams.plotsFormat      = 'vtk' ;
%
% Run ONSAS
%
[ matUsDynLD2D, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%
% Extract numerical solution
%linear disp
xdefNum = mesh.nodesCoords(:,1) + matUsDynLD2D(1:6:end, :) ;
ydefNum = mesh.nodesCoords(:,2) + matUsDynLD2D(3:6:end, :) ;
zdefNum = mesh.nodesCoords(:,2) + matUsDynLD2D(5:6:end, :) ;
%angular disp
thetaXdefNum = matUsDynLD2D(2:6:end, :) ;
thetaYdefNum = matUsDynLD2D(4:6:end, :) ;
thetaZdefNum = matUsDynLD2D(6:6:end, :) ;
%
%Plot deformed configurations at different times
timVec = linspace(0, analysisSettings.finalTime, size( matUsDynLD2D, 2 ) ) ;
numTimesToPlot = 4 ;
timeIndexToPlot = floor( linspace( 1, size( matUsDynLD2D, 2 ), numTimesToPlot ) ) ; 
%
% Plot time evolution 
% select node and dof 
node = numElements2DLDDynamic +1 ;
dofNodeY = 3 ;
fig7 = figure(7) ;
hold on
grid on
plot ( timVec, matUsDynLD2D( (node-1)*6 + dofNodeY, :), 'linewidth', lw, 'markersize' , ms ) 
%
% Add legend to yA displacements plot
plot ( timVec, ones(size(timVec, 2),1) * matUsLD2DStatic( (numElementsLD2DStatic)*6 + dofNodeY, end), 'g:', 'linewidth', lw, 'markersize' , ms ) 
legend('Dynamic solution', 'Static solution', 'location', 'South' ) ;
labx = xlabel(' time (s)');    laby=ylabel('Displacements ');
set(legend, 'linewidth' , axislw    , 'fontsize', legendFontSize )  ;
set(gca   , 'linewidth' , axislw    , 'fontsize', curveFontSize )   ;
set(labx  , 'FontSize'  , axisFontSize) ; 
set(laby  , 'FontSize'  , axisFontSize) ;
namefig5 = strcat(folderLD2DDynamicPath, 'uyA.png') ;
print( fig5, namefig5, '-dpng' ) ;
cd(accDir)
% close(2)
% elapsed time in seconds
elapTime = toc ;
fprintf( strcat('\n The elapsed execution time was ', num2str(elapTime/60), '  minutes', '\n' ) );