%md# Aerodynamic linear static cantilever beam example
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
E = 70e9 ;  nu = 0.3 ; rho = 700 ; G = E / (2 * (1+nu)) ;
% geometrical scalar parameters
l = 20 ; dext = .5 ;  b = 1e-3  ; dint  = dext - 2*b ;
A = pi * (dext^2 - dint^2) / 4  ;
J = pi * (dext^4 - dint^4) / 32 ; Iyy = J/2 ; Izz = Iyy ;
Irho = diag([J Iyy Izz],3,3);
% the number of elements of the mesh
numElements = 20 ;
%md##Numerical solution
%md### MEBI parameters
%md
%md### materials
%md Since the example contains only aeroFoone rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = 'linearElastic' ;
% materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [1 A J Iyy Izz Irho(1,1) Irho(2,2) Irho(3,3)] ;
%md The drag and lift section function names are:
numGaussPoints  = 3 ;
formulationType = 4 ;
elements(2).elemTypeAero   = [0 dext 0 numGaussPoints formulationType ];
elements(2).userDragCoef   = 'dragCoefNonLinear'   ;
elements(2).userLiftCoef   = 'liftCoefNonLinear'   ;
elements(2).userMomentCoef = 'momentCoefNonLinear' ;
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
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
% mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
%md the following case only differs in the boundary condition and the node number
for i=1:numElements
  conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   .1            ;
analysisSettings.finalTime     =   1             ;%the final time must achive 20m&s
analysisSettings.stopTolDeltau =   1e-10         ;
analysisSettings.stopTolForces =   1e-10         ;
analysisSettings.stopTolIts    =   20            ;
%md the name of the wind velocity function is: 
analysisSettings.userWindVel   = 'windVelNonLinear';
%md geometrical nonlinearity in the wind force is not taken into account in this example:
%mdloadsSettings
analysisSettings.booleanSelfWeight      = false ;
analysisSettings.geometricNonLinearAero = true ;
%md
%md## otherParams
otherParams.problemName = 'onsasExample_nonLinearStaticCantilever';
otherParams.controlDofs = [ numElements+1  4 ] ;
otherParams.plotsFormat = 'vtk' ;
%md In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%md
%md## Assamble Julia solution
%md
%mdread julia solution
xJulia         = load ( 'output/solJDiffEq_xcords.txt' ) ;
dSolJulia      = zeros( 6*size( xJulia, 2 ),1          ) ;
uYJulia        = load ( 'output/solJDiffEq_uz.txt'     ) ;
thetaZdefJulia = load ( 'output/solJDiffEq_thetaY.txt' ) ;

%md fill sol vector
dSolJulia(3:6:end) = -1 * uYJulia';
dSolJulia(6:6:end) = -1 * thetaZdefJulia';
% dSolJulia(4:6:end) = thetaYdefJulia';
% dSolJulia(4:6:end) = thetaYdefJulia';

%md create and fill solution with the same size of the vector lacating their cooridnates
dSol = zeros( (numElements + 1) * 6 ,1 );
numElemJulia = size(dSolJulia(1:6:end)) - 1;

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
betaRel = acos(dot(elements(2).elemTypeAero(1:3) , [0 0 1] ));


if isfield(elements(2), 'userDragCoef')
  c_d = feval(elements(2).userDragCoef, betaRel);
else
  c_d = 0;
end
if isfield(elements(2), 'userLiftCoef')
  c_l = feval(elements(2).userLiftCoef, betaRel);
else
  c_l = 0;
end
if isfield(elements(2), 'userMomentCoef')
  c_m = feval(elements(2).userMomentCoef, betaRel);
else
  c_m = 0;
end

%mdget wind velocity
windVel = feval(analysisSettings.userWindVel, mesh.nodesCoords(1,1), analysisSettings.finalTime) ;
%mdcaracteristicDimension
dimCaracteristic = norm(elements(2).elemTypeAero(1:3)) ;

%dynamic presure
q = 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2) ;
%loads per unit of length  
qz = q * c_l * dimCaracteristic ; 
qy = q * c_d * dimCaracteristic ; 
qm = q * c_m * dimCaracteristic ; 

%reference coordinates
xref = mesh.nodesCoords(:,1) ;
yref = mesh.nodesCoords(:,2) ;
zref = mesh.nodesCoords(:,3) ;

%Analytic x vector
sizeAnalyticX = 100 ;
xanal = linspace( 0, l , sizeAnalyticX )' ;

% linear disp
ydefAnalyticLin = qy / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4);
zdefAnalyticLin = qz / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4);
% angular disp
thetaXAnalyticLin = qm   / (2 * (Izz + Iyy) * G) * ( l^2  - ( xanal - l).^2 )  ;
thetaYAnalyticLin = -qz  / (6*E*Iyy) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;
thetaZAnalyticLin = qy   / (6*E*Izz) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3) ;

%linear disp
xdefJulia = xJulia' +  dSolJulia(1:6:end) ;
ydefJulia = dSolJulia(3:6:end)            ;
zdefJulia = dSolJulia(5:6:end)            ;
%angular disp
thetaXdefJulia = dSolJulia(2:6:end) ;
thetaYdefJulia = dSolJulia(4:6:end) ;
thetaZdefJulia = dSolJulia(6:6:end) ;

% Load numerical solution
%linear disp
xdefNum = mesh.nodesCoords(:,1) + matUs(1:6:end, end);
ydefNum = mesh.nodesCoords(:,2) + matUs(3:6:end, end);
zdefNum = mesh.nodesCoords(:,2) + matUs(5:6:end, end);
%angular disp
thetaXdefNum = matUs(2:6:end, end) ;
thetaYdefNum = matUs(4:6:end, end) ;
thetaZdefNum = matUs(6:6:end, end) ;

% Plot parameters:
lw = 5 ; ms = 5 ;
% labels parameters:
labelTitle= [' Validating solution with ' num2str(numElements) ' elements' ];
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ;    

% Plot linear displacements
figure
hold on  
grid on
% plot(xdefNum  , zdefNum,         'ro', 'linewidth', lw, 'markersize', ms + 8 );
% plot(xanal    , zdefAnalyticLin, 'r:', 'linewidth', lw, 'markersize', ms     );
% plot(xdefJulia, zdefJulia,       'r-', 'linewidth', lw, 'markersize', ms     );
plot(xdefNum    , ydefNum,        'bo' , 'linewidth', lw, 'markersize', ms+5   );
plot(xanal      , ydefAnalyticLin,'b:' , 'linewidth', lw, 'markersize', ms     );
plot(xdefJulia  , ydefJulia,      'b-' , 'linewidth', lw, 'markersize', ms     );
% legend('z_n', 'z_{lin}','z_{sa}','y_n', 'y_{lin}', 'y_{sa}', 'location', 'north')
% legend('z_n', 'z_{sa}','y_n', 'y_{sa}', 'location', 'north')
legend('y_n', 'y_{lin}', 'y_{sa}', 'location', 'north')
% legend('z_n', 'z_{lin}', 'z_{sa}', 'location', 'north')
labx=xlabel(' x (m)');    laby=ylabel('Displacements (m)');
title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
print('./output/linearDisp.png')

% Plot angular displacements
figure
hold on  
grid on
% plot(xdefNum,   rad2deg(thetaYdefNum),      'ro' , 'linewidth', lw, 'markersize', ms + 8);
% plot(xanal,     rad2deg(thetaYAnalyticLin), 'r:' , 'linewidth', lw, 'markersize', ms    );
% plot(xdefJulia, rad2deg(thetaYdefJulia),    'r-' , 'linewidth', lw, 'markersize', ms    );
plot(xdefNum,     rad2deg(thetaZdefNum),      'bo' , 'linewidth', lw,'markersize', ms+5   );
plot(xanal,       rad2deg(thetaZAnalyticLin), 'b:' , 'linewidth', lw, 'markersize', ms    );
plot(xdefJulia,   rad2deg(thetaZdefJulia),    'b-' , 'linewidth', lw, 'markersize', ms    );
% plot(xdefNum, rad2deg(thetaXdefNum), 'go' , 'linewidth', lw, 'markersize', ms+5);
% plot(xanal, rad2deg(thetaXAnalyticLin), 'g:' , 'linewidth', lw, 'markersize', ms);
% plot(xdefJulia, rad2deg(thetaXdefJulia), 'g-' , 'linewidth', lw, 'markersize', ms);
% legend('\theta y_n', '\theta y_{lin}', '\theta y_{sa}', '\theta z_n', '\theta z_{lin}', '\theta z_{sa}', '\theta x_n', '\theta x_{lin}','\theta x_{sa}', 'location','eastoutside')
% legend('\theta y_n', '\theta y_{sa}', '\theta z_n', '\theta z_{sa}', '\theta x_n', '\theta x_{sa}', 'location','eastoutside')
legend('\theta z_n', '\theta z_{lin}', '\theta z_{sa}', 'location','eastoutside')
% legend('\theta y_n', '\theta y_{lin}', '\theta y_{sa}', 'location','eastoutside')
labx=xlabel(' x (m)');    laby=ylabel('Angle (º)');
title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
print('./output/angDisp.png')

% %Plot 3D deformed
figure
hold on
grid on
plot3(xref, yref, zref,'k-' , 'linewidth', lw+300,'markersize', ms+200);
plot3(xanal, ydefAnalyticLin, zdefAnalyticLin,'b:' , 'linewidth', lw,'markersize', ms);
plot3(xdefJulia, ydefJulia, zdefJulia,'g-' , 'linewidth', lw,'markersize', ms);
plot3(xdefNum, ydefNum, zdefNum,'ro' , 'linewidth', lw,'markersize', ms+5);
% legend('Reference config','Linear defomred config','Julia DiffEqSol defomred config',  'Numerical defomred config','location','northEast')
legend('Reference config','Linear defomred config', 'DiffEq.jl defomred config',  'Numerical defomred config', 'location','northEast')
labx=xlabel('x (m)');    laby=ylabel('y(m)'); labz=zlabel('z(m)');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
view([0.5 +0.5 -1])
% axis ( [0 L max(ydefAnalyticLin) -max(ydefAnalyticLin)/10 max(zdefAnalyticLin) -max(zdefAnalyticLin)/10 ] )
grid on
print('./output/def.png','-dpng')    