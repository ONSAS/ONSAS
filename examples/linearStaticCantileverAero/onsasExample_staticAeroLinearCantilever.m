%md# Aerodynamic linear static cantilever beam example
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
addpath( genpath( [ pwd ] ) );
% material scalar parameters
E = 70e9 ;  nu = 0.3 ; rho = 700 ; G = E / (2 * (1+nu))
% geometrical scalar parameters
l = 20 ; dext = .5 ;  b = 1e-3  ; dint  = dext - 2*b ;
A = pi * (dext^2 - dint^2) / 4  ;
J = pi * (dext^4 - dint^4) / 32 ; Iyy = J/2 ; Izz = Iyy ;
Irho = diag([J Iyy Izz],3,3);
% the number of elements of the mesh
numElements = 10 ;
%md##Numerical solution
%md### MEBI parameters
%md
%md### materials
%md Since the example contains only aeroFoone rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
materials.density = rho ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [1 A J Iyy Izz Irho(1,1) Irho(2,2) Irho(3,3)] ;
%md The drag and lift section function names are:
elements(2).elemTypeAero   = [0 dext 0];
elements(2).userDragCoef   = 'dragCoefFunction'   ;
elements(2).userLiftCoef   = 'liftCoefFunction'   ;
elements(2).userMomentCoef = 'momentCoefFunction' ;
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) 0;
boundaryConds(2).loadsBaseVals = [ 0 0 0 -1 0 0 ] ;
%md the name of the wind velocity function is: 
boundaryConds(3).userWindVel    = 'windVel';
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
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1   ] ;
%md the following case only differs in the boundary condition and the node number
for i=1:numElements,
  mesh.conecCell{ i+2,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   0.5  ;
analysisSettings.finalTime     =   0.5    ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
otherParams.problemName = 'aeroLinStaticCantilever';
otherParams.controlDofs = [ numElements+1  4 ] ;
otherParams.plotsFormat = 'vtk' ;
%md In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


%md## Verification
rhoAire = 1.2;
%evaluate drag/lift and moment coefficents
betaRel = acos(dot(elements(2).elemTypeAero , [0 0 1] ));

c_d = feval(elements(2).userDragCoef,   betaRel);
c_l = feval(elements(2).userLiftCoef,   betaRel);
c_m = feval(elements(2).userMomentCoef, betaRel);

%mdget wind velocity
windVel = feval(boundaryConds(3).userWindVel, betaRel);
%mdcaracteristicDimension
dimCaracteristic = norm(elements(2).elemTypeAero);

%dynamic presure
q = 1/2 * rhoAire * (windVel(3)^2 + windVel(2)^2);
%loads per unit of length  
qz = q * c_d * dimCaracteristic; 
qy = q * c_l * dimCaracteristic; 
qm = q * c_m * dimCaracteristic; 
%reference coordinates
xref = mesh.nodesCoords(:,1);
yref = mesh.nodesCoords(:,2);
zref = mesh.nodesCoords(:,3);

%Analytic x vector
sizeAnalyticX = 100;
xanal = linspace(0,l,sizeAnalyticX)';

%Evaluate analytical solutions
% linear disp
ydefAnalytic = qy / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4);
zdefAnalytic = qz / (24*E*Izz) * (6*l^2*xanal.^2 -4*l*xanal.^3+xanal.^4);
% angular disp
thetaXAnalytic = qm / (2 * (Izz + Iyy) * G) * ( l^2  - ( xanal - l).^2 );
thetaYAnalytic = -qz / (6*E*Iyy) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3);
thetaZAnalytic = qy / (6*E*Izz) * (3* l^2 * xanal -3*l*xanal.^2+xanal.^3);

% Load numerical solution
%linear disp
xdefNum = mesh.nodesCoords(:,1) + matUs(1:6:end,end);
ydefNum = mesh.nodesCoords(:,2) + matUs(3:6:end,end);
zdefNum = mesh.nodesCoords(:,2) + matUs(5:6:end,end);
%angular disp
thetaXdefNum = matUs(2:6:end,end);
thetaYdefNum = matUs(4:6:end,end);
thetaZdefNum = matUs(6:6:end,end);

% Plot parameters:
lw = 5 ; ms = 8 ;
% labels parameters:
labelTitle= [' Validating solution with ' num2str(numElements) ' element' ];
axislw= 2; axisFontSize = 20 ; legendFontSize = 15; curveFontSize = 15;       

% Plot linear displacements
figure
hold on  
grid on
plot(xdefNum, zdefNum,'ro' , 'linewidth', lw, 'markersize', ms);
plot(xanal, zdefAnalytic,'r-' , 'linewidth', lw, 'markersize', ms);
plot(xdefNum, ydefNum, 'bo' , 'linewidth', lw,'markersize', ms);
plot(xanal, ydefAnalytic,'b-' , 'linewidth', lw, 'markersize', ms);
legend('z_n', 'z_a','y_n', 'y_a', 'location', 'north')
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
plot(xdefNum, rad2deg(thetaYdefNum), 'ro' , 'linewidth', lw, 'markersize', ms);
plot(xanal, rad2deg(thetaYAnalytic), 'r-' , 'linewidth', lw, 'markersize', ms);
plot(xdefNum, rad2deg(thetaZdefNum), 'bo' , 'linewidth', lw,'markersize', ms);
plot(xanal, rad2deg(thetaZAnalytic), 'b-' , 'linewidth', lw, 'markersize', ms);
plot(xdefNum, rad2deg(thetaXdefNum), 'go' , 'linewidth', lw, 'markersize', ms);
plot(xanal, rad2deg(thetaXAnalytic), 'g-' , 'linewidth', lw, 'markersize', ms);
legend('\theta y_n', '\theta y_a', '\theta z_n', '\theta z_a', '\theta x_n', '\theta x_a', 'location','eastoutside')
labx=xlabel(' x (m)');    laby=ylabel('Angle (º)');
title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
print('./output/angDisp.png')

%Plot 3D deformed

figure
hold on
grid on
plot3(xref, yref, zref,'k-' , 'linewidth', lw+300,'markersize', ms+200);
plot3(xanal, ydefAnalytic, zdefAnalytic,'r-' , 'linewidth', lw,'markersize', ms);
plot3(xdefNum, ydefNum, zdefNum,'bo' , 'linewidth', lw,'markersize', ms);
legend('Reference config','Numerical def config', 'Analytic def config','location','northEast')
labx=xlabel('x (m)');    laby=ylabel('y(m)'); labz=zlabel('z(m)');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ; set(labz, 'FontSize', axisFontSize) ;
view([0.5 +0.5 -1])
% axis ( [0 L max(ydefAnalytic) -max(ydefAnalytic)/10 max(zdefAnalytic) -max(zdefAnalytic)/10 ] )
grid on
print('./output/def.png','-dpng')    
