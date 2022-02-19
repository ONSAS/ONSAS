%md# Wind turbine example
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );
%
% General  problem parameters
%----------------------------
% material scalar parameters
E = 210e6 ;  nu = 0.3 ; rho = 7850 ; G = E / (2 * (1+nu)) ;
% geometrical scalar parameters
% l = 5 ; a = 0.2 ; J = 1/3 * 0.40147 * a^4 ; Iyy = a ^ 4 / 12  ; Izz = Iyy ;  
l = 2 ; d = 0.1;  
% the number of elements of the mesh for static case
numElementsBlade = 1;
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
% nodes
elements(1).elemType = 'node'  ;
% first blade aligned with -z global axis
numGaussPoints  = 4            ;
formulCase      = 4            ;
elements(2).elemType = 'frame' ;
elements(2).elemTypeGeometry = [3 d ] ;
elements(2).elemTypeAero     = [0 0 d numGaussPoints formulCase ] ;
elements(2).userLiftCoef     = 'liftCoef'                         ;
% second blade in (z,-y) quarter 
elements(3).elemType = 'frame' ;
elements(3).elemTypeGeometry = [3 d ] ;
elements(3).elemTypeAero     = [0 d 0 numGaussPoints formulCase ] ;
elements(3).userLiftCoef     = 'liftCoef'                         ;
% third blade in (z,y) quarter 
elements(4).elemType = 'frame' ;
elements(4).elemTypeGeometry = [3 d ] ;
elements(4).elemTypeAero     = [0 d 0 numGaussPoints formulCase ] ;
elements(4).userLiftCoef     = 'liftCoef'                         ;
%
% boundaryConds
%----------------------------
% The elements are submitted to two different BC settings. The first BC corresponds to a free angle in x condition 
boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
%
% initial Conditions
%----------------------------
% homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%
% mesh parameters
mesh.nodesCoords = [ 0        0              0            ; ...
                     0  l*sin( pi )        l*cos( pi )    ; ...
                     0  l*sin( pi/3  )     l*cos( pi/3 )  ; ... 
                     0  l*sin( 4*pi/3 )   -l*cos( 4*pi/3 ); ] 

mesh.conecCell         = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0 0   1 2 ] ;
mesh.conecCell{ 3, 1 } = [ 1 3 0 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 4 0 0   1 4 ] ;
%
% analysisSettings
% -------------------------------------
analysisSettings.finalTime              =   210     ;
analysisSettings.deltaT                 =   5     ;
analysisSettings.methodName             = 'alphaHHT';
analysisSettings.alphaHHT               =  -0.05    ;
analysisSettings.stopTolIts             =   30      ;
analysisSettings.geometricNonLinearAero = true      ;
analysisSettings.booleanSelfWeight      = false     ;
%
% add wind veloctiy into analysisSettings struct
analysisSettings.userWindVel = 'windVel' ;
%
% otherParams
%----------------------------
otherParams.problemName = strcat( 'onsasExample_simpleWindTurbine' ) ;
otherParams.plotsFormat = 'vtk' ;
%
% Execute ONSAS
% ----------------------------
[ matUs, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%md
%md## Verification
%mdcompute solution by the second carindal, firt the wind parameters are lodaded
rhoA = 1.225 ; c_l = feval('liftCoef', 0) ; vwind = feval('windVel', 0,0) ;
%md lift load per unit of length: 
fl = 1 / 2 * c_l * rhoA * norm(vwind) ^ 2 * d ;
%md the total moment induced in node 1 in x direction for is the sum for three blades: 
moment1x = 3 * fl * l * l / 2 ;
%md then the angular moment is:
bladeMass = rho * l * pi * d ^2 /4 ; 
Jrho =  3 * 1/3 * bladeMass  * l ^ 2 ; 
angleXnode1 = @(t)  moment1x / Jrho / 2 * t .^ 2 ;
%md numercial time vector is given by:
timeVec = linspace(0, analysisSettings.finalTime, size(matUs, 2) ) ;
%md numercial rotation angle is:
dofAngleXnode1 = 2 ;
angleXnode1Numeric = -matUs(dofAngleXnode1,:) ;
%md analytical rotation angle is:
angleXnode1Analytic = angleXnode1(timeVec) ;
%md
%md## Verification
%md
verifBoolean = norm( angleXnode1Numeric - angleXnode1Analytic )  ...
                    < ( norm( angleXnode1Numeric ) * 1e-2 ) ;
%md
%md## Plots
%md
lw = 2.0 ; ms = 10; plotfontsize = 22 ;
spanPlotTime = 2 ;
figure
plot( timeVec(1:spanPlotTime:end), angleXnode1Analytic(1:spanPlotTime:end) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( timeVec(1:spanPlotTime:end), angleXnode1Numeric(1:spanPlotTime:end), 'ko' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time(s)');   laby = ylabel('$\theta_x node 1$') ;
legend('analytic','numeric','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/verifSimpleWindTurbine.png','-dpng')