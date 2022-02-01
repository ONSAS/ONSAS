%md# Aerodynamic linear static cantilever beam example
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
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
%md Since the example contains only aeroFoone rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = 'linearElastic' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [3 d] ;
%md The drag and lift section function names are:
numGaussPoints  = 3 ;
%mdTest different formulations
for formulCase = [2]
  elements(2).elemTypeAero   = [0 d 0 numGaussPoints formulCase ] ;
  elements(2).userDragCoef   = 'dragCoefNonLinear'   ;
  % elements(2).userLiftCoef   = 'liftCoefNonLinear'   ;
  % elements(2).userMomentCoef = 'momentCoefNonLinear' ;
  %md
  %md### boundaryConds
  %md
  %md The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
  boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
  boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
  %md
  %md### initial Conditions
  %md homogeneous initial conditions are considered, then an empty struct is set:
  initialConds                = struct() ;
  %md
  %md### mesh parameters
  %mdThe coordinates of the nodes of the mesh are given by the matrix:
  mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
  %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
  mesh.conecCell = { } ;
  %md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  % mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
  %md the following case only differs in the boundary condition and the node number
  for i=1:numElements
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  end
  %md-------------------------------------
  %md## Static small displacements case
  %md -------------------------------------
  %md### analysisSettings
  analysisSettings.methodName    = 'newtonRaphson' ;
  analysisSettings.deltaT        =   .1            ;
  analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
  analysisSettings.stopTolDeltau =   1e-8          ;
  analysisSettings.stopTolForces =   1e-8          ;
  analysisSettings.stopTolIts    =   30            ;
  %md the name of the wind velocity function is: 
  analysisSettings.userWindVel   = 'windVelNonLinearStaticSD';
  %md geometrical nonlinearity in the wind force is not taken into account in this example:
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
  solCell     = dlmread ('output/solJDiffEq.csv', ',', 1, 0) ;
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

  % linear disp
  ydefAnalyticLin = @(windVel) 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_d * dimCaracteristic / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4) ;
  zdefAnalyticLin = @(windVel) 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_l * dimCaracteristic  / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4) ;
  % angular disp
  thetaXAnalyticLin = @(windVel)  1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_m * dimCaracteristic   / (2 * (Izz + Iyy) * G) * ( l^2  - ( xanal - l).^2 )  ;
  thetaYAnalyticLin = @(windVel) -1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_l * dimCaracteristic   / (6*E*Iyy) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;
  thetaZAnalyticLin = @(windVel)  1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_d * dimCaracteristic   / (6*E*Izz) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;

  %linear disp
  xdefJulia = xJulia +  dSolJulia(1:6:end) ;
  ydefJulia = dSolJulia(3:6:end)           ;
  zdefJulia = dSolJulia(5:6:end)           ;
  %angular disp
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
  % Plot linear displacements
  figure(1)
  hold on  
  grid on
  plot(xdefNum    , ydefNum,        'bo' , 'linewidth', lw, 'markersize', ms+5   );
  plot(xanal      , ydefAnalyticLin(windVel),'b:' , 'linewidth', lw, 'markersize', ms     );
  plot(xdefJulia  , ydefJulia,      'b-' , 'linewidth', lw, 'markersize', ms     );
  legend('y_n', 'y_{lin}', 'y_{sa}', 'location', 'north')
  labx=xlabel(' x ');    laby=ylabel('Displacements ');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/linDispSD.png') ) ;
  print( strcat(folderSDpath, 'linDispSD.png') ) ;
  close(1)
  % Plot angular displacements
  figure(2)
  hold on  
  grid on
  plot(xdefNum,     rad2deg(thetaZdefNum),      'bo' , 'linewidth', lw,'markersize', ms+5   );
  plot(xanal,       rad2deg( thetaZAnalyticLin(windVel) ), 'b:' , 'linewidth', lw, 'markersize', ms    );
  plot(xdefJulia,   rad2deg(thetaZdefJulia),    'b-' , 'linewidth', lw, 'markersize', ms    );
  legend('\theta z_n SD', '\theta z_{lin} SD', '\theta z_{sa} SD', 'location','eastoutside')
  labx=xlabel(' x ');    laby=ylabel('Angle');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/angDispSD.png') ) ;
  print( strcat(folderSDpath, 'angDispSD.png') ) ;
  close(2)
  % %Plot 3D deformed
  figure(3)
  hold on
  plot(xref,      yref, 'k-' , 'linewidth', lw+300,'markersize', ms+200);
  plot(xanal,     ydefAnalyticLin(windVel), 'b:' , 'linewidth', lw,'markersize', ms);
  plot(xdefJulia, ydefJulia, 'g-' , 'linewidth', lw,'markersize', ms);
  plot(xdefNum,   ydefNum, 'ro' , 'linewidth', lw,'markersize', ms+5);
  legend('Reference config','Linear defomred config SD', 'DiffEq.jl defomred config SD',  'Numerical defomred config SD', 'location','northEast')
  labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
  grid on
  print( strcat('./output/', otherParams.problemName, '/defSD.png') ) ;
  print( strcat(folderSDpath, 'defSD.png') ) ;
  close(3)
  %md-------------------------------------
  %md## Static 2D large displacements case
  %md -------------------------------------
  numElements = 5 ;
  %md## mesh parameters
  %mdThe coordinates of the nodes of the mesh are given by the matrix:
  mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
  %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
  mesh.conecCell = { } ;
  %md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  % mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
  %md the following case only differs in the boundary condition and the node number
  for i=1:numElements
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  end
  materials.hyperElasModel  = '1DrotEngStrain' ;
  %md### analysisSettings
  analysisSettings.methodName    = 'newtonRaphson' ;
  analysisSettings.deltaT        =   .1            ;
  analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
  analysisSettings.stopTolDeltau =   1e-7          ;
  analysisSettings.stopTolForces =   1e-7          ;
  analysisSettings.stopTolIts    =   30            ;
  %md the name of the wind velocity function is: 
  analysisSettings.userWindVel   = 'windVelNonLinearStaticLD';
  %md geometrical nonlinearity in the wind force is not taken into account in this example:
  %mdloadsSettings
  analysisSettings.booleanSelfWeight      = false ;
  analysisSettings.geometricNonLinearAero = true  ;
  %md
  %md## otherParams
  otherParams.problemName = strcat( 'onsasExample_nonLinearCantilever_staitcLD_formulation=', num2str(formulCase) ) ;
  otherParams.controlDofs = [ numElements+1  4 ] ;
  otherParams.plotsFormat = 'vtk' ;
  %md
  [matUsLD, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
  %md
  %mdget wind velocity
  windVel = feval(analysisSettings.userWindVel, mesh.nodesCoords(1,1), analysisSettings.finalTime) ;
  % Load numerical solution
  %linear disp
  xdefNum = mesh.nodesCoords(:,1) + matUsLD(1:6:end, end) ;
  ydefNum = mesh.nodesCoords(:,2) + matUsLD(3:6:end, end) ;
  zdefNum = mesh.nodesCoords(:,2) + matUsLD(5:6:end, end) ;
  %angular disp
  thetaXdefNum = matUsLD(2:6:end, end) ;
  thetaYdefNum = matUsLD(4:6:end, end) ;
  thetaZdefNum = matUsLD(6:6:end, end) ;

  % Folder path
  folderLD2Dpath = strcat('./output/', 'LD/', '2D/') ;
  % Plot linear displacements
  figure(1)
  hold on  
  grid on
  plot(xdefNum,   ydefNum,                   'ro' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,     ydefAnalyticLin(windVel),  'r:' , 'linewidth', lw,'markersize', ms   );
  legend('y_n LD', 'y_{lin} LD', 'location', 'north')
  labx=xlabel(' x ');    laby=ylabel('Displacements ');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/linDispLD.png') ) ;
  print( strcat(folderLD2Dpath, 'static/', 'linDispLD.png') ) ;
  close(1)
  % Plot angular displacements
  figure(2)
  hold on  
  grid on
  plot(xdefNum,     rad2deg(thetaZdefNum),               'bo' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,       rad2deg(thetaZAnalyticLin(windVel)), 'b:' , 'linewidth', lw, 'markersize', ms  );
  legend('\theta z_n LD', '\theta z_{lin} LD', 'location','eastoutside')
  labx=xlabel(' x ');    laby=ylabel('Angle');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/angDispLD.png') ) ;
  print( strcat(folderLD2Dpath, 'static/', 'angDispLD.png') ) ;
  close(2)
  %Plot 3D deformed
  figure(3)
  hold on
  grid on
  plot(xref,      yref, 'k-' , 'linewidth', lw+300,'markersize', ms+200);
  plot(xanal,     ydefAnalyticLin(windVel), 'b:' , 'linewidth', lw,'markersize', ms);
  plot(xdefNum,   ydefNum, 'ro' , 'linewidth', lw,'markersize', ms+5);
  legend('Reference config','Linear defomred config LD',  'Numerical defomred config LD', 'location','northEast')
  labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/defLD.png') ) ;
  print( strcat(folderLD2Dpath, 'static/', 'defLD.png') ) ;
  close(3)
  %md-------------------------------------
  %md Dynamic 2D Case 
  %md-------------------------------------
  numElements = 10 ;
  %md## mesh parameters
  %mdThe coordinates of the nodes of the mesh are given by the matrix:
  mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
  %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
  mesh.conecCell = { } ;
  %md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  % mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
  %md the following case only differs in the boundary condition and the node number
  for i=1:numElements
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  end
  % MEBI parameters
  % analysisSettings
  %
  analysisSettings.deltaT      =   0.01    ;
  analysisSettings.finalTime   =   5       ;
  analysisSettings.methodName  = 'newmark' ;
  analysisSettings.alphaNM     =   0.25    ;
  analysisSettings.deltaNM     =   0.5     ;
  analysisSettings.userWindVel = 'windVelNonLinearDynamic';
  %md 
  %mdRun with different nodal disp dampings
  dampingVec = [ 0 1 ] ;
  for dampingCase = dampingVec
    %
    % otherParams
    %
    otherParams.nodalDispDamping = dampingCase ;
    otherParams.problemName      = strcat( 'onsasExample_nonLinearCantilever_dynamc_c=', num2str(dampingCase), '_formulation=', num2str(formulCase) ) ;
    otherParams.plotsFormat      = 'vtk' ;
    %
    % Run ONSAS
    %
    [ matUsDyn, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 

    % Extract numerical solution
    %linear disp
    xdefNum = mesh.nodesCoords(:,1) + matUsDyn(1:6:end, :) ;
    ydefNum = mesh.nodesCoords(:,2) + matUsDyn(3:6:end, :) ;
    zdefNum = mesh.nodesCoords(:,2) + matUsDyn(5:6:end, :) ;
    %angular disp
    thetaXdefNum = matUsDyn(2:6:end, :) ;
    thetaYdefNum = matUsDyn(4:6:end, :) ;
    thetaZdefNum = matUsDyn(6:6:end, :) ;
    
    %Plot deformed configurations at different times
    timVec = linspace(0, analysisSettings.finalTime, size( matUsDyn, 2 ) ) ;
    numTimesToPlot = 4 ;
    timeIndexToPlot = floor( linspace( 1, size( matUsDyn, 2 ), numTimesToPlot ) ) ; 

    % Plot 2D deformed
    legendsText = [] ;
    figure(1)
    hold on
    % plot3(xref, yref, zref,'k-' , 'linewidth', lw+300,'markersize', ms+200);
    plot(xref, yref, 'k-' , 'linewidth', lw+300,'markersize', ms+200);
    for timePlot = timeIndexToPlot
      plot(xdefNum(:, timePlot), ydefNum(:, timePlot) , 'linewidth', lw,'markersize', ms+5) ;
      % plot3(xdefNum(:, timePlot), ydefNum(:, timePlot), zdefNum(:, timePlot) , 'linewidth', lw,'markersize', ms+5) ;
      legendsText  = [ legendsText; strcat( 'deformed configuration t = ', num2str( timVec(timePlot) ), ' s' ) ] ;
    end
    % view([0 0 1])
    legend('Reference configuration',legendsText(1,:), legendsText(2,:), legendsText(3,:), legendsText(4,:), 'location','north')
    labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
    set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
    set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
    set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
    grid on
    print( strcat('./output/', otherParams.problemName, '/def.png') ) ;
    print( strcat(folderLD2Dpath, 'static/', 'def.png') ) ;
    close(1)
    % Plot time evolution 
    % select node and dof 
    node = numElements +1 ;
    dofNodeY = 3 ;
    figure(2)
    hold on
    grid on
    plot ( timVec, matUsDyn( (node-1)*6 + dofNodeY, :), 'linewidth', lw, 'markersize' , ms ) 
  
  end
    % Add legend to yA displacements plot
    figure(2)
    % legend(strcat(' uyA_with c =', num2str( dampingVec(1) ) ), strcat(' uyA c =', num2str( dampingVec(2) ), strcat(' uyA c =', num2str( dampingVec(3) ) ), 'location', 'South' ) ;
    labx = xlabel(' time (s)');    laby=ylabel('Displacements ');
    set(legend, 'linewidth' , axislw    , 'fontsize', legendFontSize )  ;
    set(gca   , 'linewidth' , axislw    , 'fontsize', curveFontSize )   ;
    set(labx  , 'FontSize'  , axisFontSize) ; 
    set(laby  , 'FontSize'  , axisFontSize) ;
    print( strcat(folderLD2Dpath, 'static/', 'uyA.png') ) ;
  
  %md-------------------------------------
  %md## Static 3D large displacements case
  %md -------------------------------------
  numElements = 5 ;
  %md## mesh parameters
  %mdThe coordinates of the nodes of the mesh are given by the matrix:
  mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
  %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
  mesh.conecCell = { } ;
  %md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  % mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
  %md the following case only differs in the boundary condition and the node number
  for i=1:numElements
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  end
  materials.hyperElasModel  = '1DrotEngStrain' ;
  %md### add lif
  elements(2).userLiftCoef   = 'liftCoefNonLinear'   ;
  elements(2).userMomentCoef = 'momentCoefNonLinear' ;
  
  %md### analysisSettings
  analysisSettings.methodName    = 'newtonRaphson' ;
  analysisSettings.deltaT        =   .1            ;
  analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
  analysisSettings.stopTolDeltau =   1e-7          ;
  analysisSettings.stopTolForces =   1e-7          ;
  analysisSettings.stopTolIts    =   40            ;
  %md the name of the wind velocity function is: 
  analysisSettings.userWindVel   = 'windVelNonLinearStaticLD';
  %md geometrical nonlinearity in the wind force is not taken into account in this example:
  %mdloadsSettings
  analysisSettings.booleanSelfWeight      = false ;
  analysisSettings.geometricNonLinearAero = true  ;
  %md
  %md## otherParams
  otherParams.problemName = strcat( 'onsasExample_nonLinearCantilever_staitcLD_formulation=', num2str(formulCase),' 3D' ) ;
  otherParams.controlDofs = [ numElements+1  4 ] ;
  otherParams.plotsFormat = 'vtk' ;
  %md
  %md In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
  [matUsLD3D, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
  %md
  %Plot solution
  folderLD3Dpath = strcat('./output/', 'LD/', '3D/') ;
  % Build analytical solution 3D
    if isfield(elements(2), 'userDragCoef')
    c_d = feval(elements(2).userDragCoef, 0) ;
  else
    c_d = 0;
  end
  if isfield(elements(2), 'userLiftCoef')
    c_l = feval(elements(2).userLiftCoef, 0) ;
  else
    c_l = 0;
  end
  if isfield(elements(2), 'userMomentCoef')
    c_m = feval(elements(2).userMomentCoef, 0) ;
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

  % linear disp
  ydefAnalyticLin = @(windVel) 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_d * dimCaracteristic / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4) ;
  zdefAnalyticLin = @(windVel) 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_l * dimCaracteristic  / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4) ;
  % angular disp
  thetaXAnalyticLin = @(windVel)  1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_m * dimCaracteristic   / (2 * (Izz + Iyy) * G) * ( l^2  - ( xanal - l).^2 )  ;
  thetaYAnalyticLin = @(windVel) -1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_l * dimCaracteristic   / (6*E*Iyy) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;
  thetaZAnalyticLin = @(windVel)  1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) * c_d * dimCaracteristic   / (6*E*Izz) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;
  
  % Extract numerical solution
  xdefNum = mesh.nodesCoords(:,1) + matUsLD3D(1:6:end, end) ;
  ydefNum = mesh.nodesCoords(:,2) + matUsLD3D(3:6:end, end) ;
  zdefNum = mesh.nodesCoords(:,2) + matUsLD3D(5:6:end, end) ;
  %angular disp
  thetaXdefNum = matUsLD3D(2:6:end, end) ;
  thetaYdefNum = matUsLD3D(4:6:end, end) ;
  thetaZdefNum = matUsLD3D(6:6:end, end) ;
  
  % Plot linear displacements
  figure(1)
  hold on  
  grid on
  plot(xdefNum,   ydefNum,                  'ro' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,     ydefAnalyticLin(windVel), 'r:' , 'linewidth', lw,'markersize', ms   );
  plot(xdefNum,   zdefNum,                  'bo' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,     zdefAnalyticLin(windVel), 'b:' , 'linewidth', lw,'markersize', ms   );
  legend('uy_n', 'uy_{lin}', 'uy_n', 'uz_{lin}', 'location', 'north')
  labx=xlabel(' x ');    laby=ylabel('Displacements ');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/linDispLD3D.png') ) ;
  print( strcat(folderLD3Dpath, 'static/', 'linDispLD3D.png') ) ;
  close(1)
  % Plot angular displacements
  figure(2)
  hold on  
  grid on
  plot(xdefNum, rad2deg(thetaZdefNum),               'ro' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,   rad2deg(thetaZAnalyticLin(windVel)), 'r:' , 'linewidth', lw, 'markersize', ms  );
  plot(xdefNum, rad2deg(thetaYdefNum),               'bo' , 'linewidth', lw,'markersize', ms+5 );
  plot(xanal,   rad2deg(thetaYAnalyticLin(windVel)), 'b:' , 'linewidth', lw, 'markersize', ms  );
  legend('\theta z_n', '\theta z_{lin}', '\theta y_n', '\theta y_{lin}', 'location','eastoutside')
  labx=xlabel(' x ');    laby=ylabel('Angle');
  % title (labelTitle)
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
  print( strcat('./output/', otherParams.problemName, '/angDispLD3D.png') ) ;
  print( strcat(folderLD3Dpath, 'static/', 'angDispLD3D.png') ) ;
  close(2)
  %Plot 3D deformed
  figure(3)
  hold on
  grid on
  plot3(xref,    yref,                     zref,                     'k-' , 'linewidth', lw+300,'markersize', ms+200);
  plot3(xanal,   ydefAnalyticLin(windVel), zdefAnalyticLin(windVel), 'b:' , 'linewidth', lw,'markersize', ms        );
  plot3(xdefNum, ydefNum,                  zdefNum,                  'r' , 'linewidth', lw,'markersize', ms+5      );
  legend('Reference config','Linear defomred config LD',  'Numerical defomred config LD', 'location','northEast')
  labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
  set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
  set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
  set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
  view([-1.5 -3 1])
  print( strcat('./output/', otherParams.problemName, '/defLD3D.png') ) ;
  print( strcat(folderLD3Dpath, '/defLD3D.png') ) ;
  close(3)
  %md-------------------------------------
  %md Dynamic 3D Case 
  %md-------------------------------------
  numElements = 10 ;
  %md## mesh parameters
  %mdThe coordinates of the nodes of the mesh are given by the matrix:
  mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
  %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
  mesh.conecCell = { } ;
  %md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  % mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
  %md the following case only differs in the boundary condition and the node number
  for i=1:numElements
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  end
  % MEBI parameters
  % analysisSettings
  %
  analysisSettings.finalTime   =   5        ;
  analysisSettings.deltaT      =   analysisSettings.finalTime / 100 ;
  analysisSettings.methodName  = 'newmark'  ;
  analysisSettings.alphaNM     =   0.25     ;
  analysisSettings.deltaNM     =   0.5      ;
  analysisSettings.methodName  =  'alphaHHT';
  analysisSettings.alphaHHT    =  -.05      ;
  analysisSettings.userWindVel = 'windVelNonLinearDynamic';
  %md 
  %mdRun with different nodal disp dampings
  dampingVec = [ 0 ] ;
  for dampingCase = dampingVec
    %
    % otherParams
    %
    otherParams.nodalDispDamping = dampingCase ;
    otherParams.problemName      = strcat( 'onsasExample_nonLinearCantilever_dynamc_c=', num2str(dampingCase), '_formulation=', num2str(formulCase), ' 3D') ;
    otherParams.plotsFormat      = 'vtk' ;
    %
    % Run ONSAS
    %
    [ matUsDyn3D, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 

    % Extract numerical solution
    %linear disp
    xdefNum = mesh.nodesCoords( :,1 ) + matUsDyn3D( 1:6:end, : ) ;
    ydefNum = mesh.nodesCoords( :,2 ) + matUsDyn3D( 3:6:end, : ) ;
    zdefNum = mesh.nodesCoords( :,2 ) + matUsDyn3D( 5:6:end, : ) ;
    %angular disp
    thetaXdefNum = matUsDyn3D( 2:6:end, : ) ;
    thetaYdefNum = matUsDyn3D( 4:6:end, : ) ;
    thetaZdefNum = matUsDyn3D( 6:6:end, : ) ;
    
    %Plot deformed configurations at different times
    timVec = linspace(0, analysisSettings.finalTime, size( matUsDyn3D, 2 ) ) ;
    numTimesToPlot = 4 ;
    timeIndexToPlot = floor( linspace( 1, size( matUsDyn3D, 2 ), numTimesToPlot ) ) ; 
    % Plot 3D deformed
    legendsText = [] ;
    figure(1)
    hold on
    plot3(xref, yref, zref,'k-' , 'linewidth', lw+300,'markersize', ms+200);
    for timePlot = timeIndexToPlot
      plot3(xdefNum(:, timePlot), ydefNum(:, timePlot), zdefNum(:, timePlot) , 'linewidth', lw,'markersize', ms+5) ;
      legendsText  = [ legendsText; strcat( 'deformed configuration t = ', num2str( timVec(timePlot) ), ' s' ) ] ;
    end
    % view([0 0 1])
    legend('Reference configuration',legendsText(1,:), legendsText(2,:), legendsText(3,:), legendsText(4,:), 'location','north')
    labx=xlabel('x ');    laby=ylabel('y'); labz=zlabel('z');
     view([-1.5 -3 1])
    set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
    set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
    set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
    grid on
    print( strcat('./output/', otherParams.problemName, '/def3D.png') ) ;
    print( strcat(folderLD3Dpath, '/dynamic', '/def3D.png') ) ;
    close(1)
    % Plot time evolution 
    % select node and dof 
    node = numElements +1 ;
    dofNodeY = 3 ;
    dofNodeZ = 5 ;
    figure(2)
    hold on
    plot ( timVec, matUsDyn3D( (node-1)*6 + dofNodeY, :), 'linewidth', lw, 'markersize' , ms ) 
    figure(3)
    hold on 
    plot ( timVec, matUsDyn3D( (node-1)*6 + dofNodeZ, :), 'linewidth', lw, 'markersize' , ms ) 
  end
    % Add legend to yA displacements plot
    figure(2)
    grid on
    plot ( timVec, ones(size(timVec,2),1)* matUsLD3D( (6-1)*6 + dofNodeY, end), 'r:', 'linewidth', lw, 'markersize' , ms ) 
    legend(strcat(' uyA c =', num2str( dampingVec(1) ) ), strcat(' uyA c =', num2str( dampingVec(2) ) ), strcat(' uyA c =', num2str( dampingVec(3) ) ), ' uyA static' , 'location', 'South' ) ;
    labx = xlabel(' time (s)');    laby=ylabel('Displacements ');
    set(legend, 'linewidth' , axislw    , 'fontsize', legendFontSize )  ;
    set(gca   , 'linewidth' , axislw    , 'fontsize', curveFontSize )   ;
    set(labx  , 'FontSize'  , axisFontSize) ; 
    set(laby  , 'FontSize'  , axisFontSize) ;
    print( strcat(folderLD3Dpath, 'dynamic/', 'uyA.png') ) ;
    % Add legend to yA displacements plot
    figure(3)
    grid on
    plot ( timVec, ones(size(timVec,2),1)* matUsLD3D( (6-1)*6 + dofNodeZ, end), 'r:', 'linewidth', lw, 'markersize' , ms ) 
    legend(strcat(' uzA c =', num2str( dampingVec(1) ) ), strcat(' uzA c =', num2str( dampingVec(2) ) ), strcat(' uzA c =', num2str( dampingVec(3) ) ), ' uzA static' , 'location', 'South' ) ;
    labx = xlabel(' time (s)');    laby=ylabel('Displacements ');
    set(legend, 'linewidth' , axislw    , 'fontsize', legendFontSize )  ;
    set(gca   , 'linewidth' , axislw    , 'fontsize', curveFontSize )   ;
    set(labx  , 'FontSize'  , axisFontSize) ; 
    set(laby  , 'FontSize'  , axisFontSize) ;
    print( strcat(folderLD3Dpath, 'dynamic/', 'uzA.png') ) ;

end