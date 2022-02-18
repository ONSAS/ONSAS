%md# Wind turbine example
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );
%
% General  problem parameters
%----------------------------
% material scalar parameters
E = 210e9 ;  nu = 0.3 ; rho = 7850 ; G = E / (2 * (1+nu)) ;
% geometrical scalar parameters
% l = 5 ; a = 0.2 ; J = 1/3 * 0.40147 * a^4 ; Iyy = a ^ 4 / 12  ; Izz = Iyy ;  
l = 5 ; d = 0.3;  
% the number of elements of the mesh for static case
numElementsBlade = 1;

% second carindal:
rhoA = 1.225 ;
c_l = feval('liftCoefS809', 0) ;
vwind = feval('windVelDynamic', 0,0) ;
fl = 1 / 2 * rhoA * norm(vwind) ^ 2 * d ;
axialMoment = 3 * fl * l * l / 2 ;

mass = rho * l * pi * d ^2 /4 ; 
Jrho =  3 * mass  * l ^ 2 ; 
angularAcel = axialMoment / Jrho ;
timeT = 1 ;
angleTimeT =  angularAcel * timeT ^ 2 / 2 ;

%
% materials
%----------------------------
% Since the example contains only aeroFoone rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%
% elements
%----------------------------
%Two different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
% for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [2 d d] ;
elements(3).elemType = 'frame' ;
elements(3).elemTypeGeometry = [2 d d] ;
% Test different formulations
for formulCase = [2]
  % boundaryConds
  %----------------------------
  % The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
  boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
  boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
  %
  % initial Conditions
  %----------------------------
  % homogeneous initial conditions are considered, then an empty struct is set:
  initialConds = struct() ;
  %
  % mesh parameters
  %----------------------------
  %The coordinates of the nodes of the mesh are given by the matrix:
  localAxialBladeCords = ( 0:( numElementsBlade ) )'*l/ numElementsBlade ;
  nodesLocalBladeZ = [ zeros( numElementsBlade +1,2) -localAxialBladeCords ] ;
  nodesLocalBlade120 = [ zeros( numElementsBlade +1,1), +sin( deg2rad(0) )*localAxialBladeCords,  +cos( deg2rad(0) )*(localAxialBladeCords) ] ;
  nodesLocalBlade120(1,:) = [] ;
  nodesLocalBlade240 = [ zeros( numElementsBlade +1,1) -sin( deg2rad(60) )*localAxialBladeCords +cos( deg2rad(60) )*(localAxialBladeCords) ] ;
  nodesLocalBlade240(1,:) = [] ;
  %The final nodes coordinates matrix is:
  % the mesh conecitvity cell is
  % mesh.nodesCoords = [   nodesLocalBladeZ;  nodesLocalBlade240; nodesLocalBlade120; ] ;
  mesh.nodesCoords = [   nodesLocalBladeZ;  nodesLocalBlade120; ] ;
  mesh.conecCell = { } ;
  mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
  %The conec cell is assamble as:
  for i=1:numElementsBlade
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
    if i == 1 
      conecElemMatrix( i + numElementsBlade  , : ) = [ 1 3 0 0  1   numElementsBlade + 1 + i ] ;
      % conecElemMatrix( i + 2*numElementsBlade, : ) = [ 1 2 0 0  1 2*numElementsBlade + 1 + i ] ;
      mesh.conecCell{ i  + numElementsBlade + 1 , 1 } = [ 1 3 0 0  1     numElementsBlade + 1 + i ] ;
      % mesh.conecCell{ i  + 2*numElementsBlade +1, 1 } = [ 1 2 0 0  1   2*numElementsBlade + 1 + i ] ;
    elseif i > 1
      conecElemMatrix( i + numElementsBlade,   : ) = [ 1 2 0 0    numElementsBlade + i    numElementsBlade + i + 1  ] ;
      % conecElemMatrix( i + 2*numElementsBlade, : ) = [ 1 3 0 0  2*numElementsBlade + i  2*numElementsBlade + i + 1  ] ;
      mesh.conecCell{ i  +   numElementsBlade +1 , 1 } = [ 1 2 0 0    numElementsBlade + i   numElementsBlade + i + 1  ] ;
      % mesh.conecCell{ i  + 2*numElementsBlade +1 , 1 } = [ 1 2 0 0  2*numElementsBlade + i 2*numElementsBlade + i + 1  ] ;
    end
  end
  figure
  hold on
  plot(nodesLocalBladeZ(:,2), nodesLocalBladeZ(:,3),  'linewidth', 5 )
  plot( [ 0; nodesLocalBlade120(:,2) ], [ 0; nodesLocalBlade120(:,3) ], 'linewidth', 5 )
  plot( [ 0; nodesLocalBlade240(:,2) ], [ 0; nodesLocalBlade240(:,3) ], 'linewidth', 5 )
  %-------------------------------------
  % Static case
  % -------------------------------------
  % analysisSettings
  %---------------------------- 
  %numericalMethodSettings static case
  analysisSettings.methodName    = 'newtonRaphson'                    ;
  analysisSettings.finalTime     =   1                                ;
  analysisSettings.deltaT        =   analysisSettings.finalTime / 10 ;
  analysisSettings.stopTolDeltau =   1e-6                             ;
  analysisSettings.stopTolForces =   1e-6                             ;
  analysisSettings.stopTolIts    =   30                               ;
  analysisSettings.booleanSelfWeight = true                           ;
  %
  % otherParams
  %----------------------------
  otherParams.problemName      = strcat( 'onsasExample_windTurbine_', 'SelfWeight',...
                                         '_formulation=', num2str(formulCase), ' 3D' ) ;
  otherParams.plotsFormat      = 'vtk' ;
  % Execute ONSAS
  % ----------------------------
  % [ matUsStatic, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
  %-------------------------------------
  % Dynamic case
  % -------------------------------------
  %
  % The drag and lift section function names are:
  numGaussPoints  = 1 ;
  elements(2).elemTypeAero   = [0 0 d numGaussPoints formulCase ] ;
  elements(2).userDragCoef   = 'dragCoefS809'   ;
  elements(2).userLiftCoef   = 'liftCoefS809'   ;
  elements(2).userMomentCoef = 'momentCoefS809'   ;
  elements(3).elemTypeAero   = [0 0 -d numGaussPoints formulCase ] ;
  elements(3).userDragCoef   = 'dragCoefS809'   ;
  elements(3).userLiftCoef   = 'liftCoefS809'   ;
  elements(3).userMomentCoef = 'momentCoefS809' ;
  %
  % analysisSettings
  % -------------------------------------
  analysisSettings.finalTime         =   100     ;% This value must be manually introduced into windVelNonLinearDynamic to achive max wind vl
  analysisSettings.deltaT            =   0.5   ;
  analysisSettings.methodName        = 'alphaHHT';
  analysisSettings.alphaHHT          =  -0.05    ;
  analysisSettings.stopTolIts        =   10      ;
  analysisSettings.geometricNonLinearAero = true  ;
  % 
  % Run with different velocity cases
  %----------------------------
  % matUsExamplesDyn = [] ;
  velocitiyCases = ['windVelDynamic'] ;
  velCaseIndex = 1 ;
  % for velCaseIndex = 1 : size(velocitiyCases, 1)
    %
    % otherParams
    %----------------------------
    analysisSettings.userWindVel = velocitiyCases(velCaseIndex, :) ;
    otherParams.problemName      = strcat( 'onsasExample_windTurbine_', analysisSettings.userWindVel, '_formulation=', num2str(formulCase), ' 3D' ) ;
    otherParams.plotsFormat      = 'vtk' ;
    % Execute ONSAS
    % ----------------------------
    [ matUsDyncurrVel, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
    % fill matrix of matUs for different c cases
    % matUsExamples3DDyn = [ matUsExamples3DDyn matUsDyn3DcurrVel ] ;
  end