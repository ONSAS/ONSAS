% ONSAS pnonlinear pendulum

close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
global massratio;
global AMBool;
%md Case 1: the fluid is still: it is only modelled by an added inertia of the solid, the fluid propertiies are opnly defined by the massratio = rho_structure/rho_fluid
otherParams.problemName     = 'AMCase1';
%md AMBool = false because the fluid is not defined here(see case 2)
AMBool = false;
% material scalar parameters
EA = 1e8     ; A = 0.1  ; d = 2*sqrt(A/pi);
E   = EA / A ; nu = 0   ;
l0  = 3.0443 ; T = 4.13 ;
m   = 10     ; g = 9.80 ;
%md
rho = 2*m / ( A * l0 )  ;
materials.density = rho ;
%md Fluid inertia is only defined by the mass ratio:
massratio = 1;
%md Moreover, the constitutive behavior considered is the Rotated Engineering strain, thus the field `hyperElasModel` is:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and truss. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The `elemType` field is then:
elements(1).elemType = 'node' ;
elements(2).elemType     = 'frame';
%mdA rectangular $2$ section is considered with sqrt(A)xsqrt(A). However this type of section has no effect in the results, because of the inertial primacy against stiffness terms. Subsequently `elemCrossSecParams` field is:
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [ d ] ;
elements(2).massMatType  = 'consistent';
%md
%md### boundaryConds
%md
%mdThe elements are submitted to a hinged condition where onnly rotation along y is allowed. then the `boundaryConds(1)` set is:
boundaryConds(1).imposDispDofs = [ 1 2 3 5 6] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0] ;
%md### initial Conditions
%mdNon homogeneous initial conditions are not considered in this example, consequently the `initialConds`  struct is set empty:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates considering a mesh of two nodes is:
angle_init = 25; % degrees
mesh.nodesCoords = [   0                    0    l0 ; ...
                   sind(angle_init)*l0  0  l0-cosd(angle_init)*l0  ] ;
%mdThe connectivity is introduced using the _conecCell_ cell. Each entry of the cell (indexed using {}) contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md Then the entry of node $1$ is introduced:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the first MEBI parameter (Material) is set as _zero_ (since nodes dont have material). The second parameter corresponds to the Element, and a _1_ is set since `node` is the first entry of the  `elements.elemType` cell. For the BC index, we consider that node $1$ is simple fixed, then the first index of the `boundaryConds` struct is used. Finally, no specific initial conditions are set for the node (0) and at the end of the vector the number of the node is included (1).
mesh.conecCell{ 2, 1 } = [ 0 1 0 0  2   ] ;
%md and for node $2$ only the boundary condition is changed, because it is lodaded.
%md Regarding the truss elements, the first material is considered, the second type of element, and no boundary conditions are applied.
mesh.conecCell{ 3, 1 } = [ 1 2 0 0  1 2 ] ;
%md
%md### analysisSettings
%mdA Newmark algorithm is used to solve this problem with the following parameters during one period $T$:
analysisSettings.deltaT        = 0.05  ;
analysisSettings.finalTime     = 2*T  ;
analysisSettings.stopTolDeltau = 1e-12 ;
analysisSettings.stopTolForces = 1e-12 ;
analysisSettings.stopTolIts    = 30    ;
%otherParams.plots_format       = 'vtk' ;
% ------------------------------------
%md### Analysis case 2: Solution using HHT with truss element, consistent mass formulation and self weight boolean activated:
analysisSettings.booleanSelfWeight = true ;
%mdIn order to validate HHT numerical method this is executed with $alphaHHT = 0$ which necessary implies that numerical results must be identical to CASE 1, since HHT is equivalent to Newmark if AplhaHHT = 0,
analysisSettings.methodName = 'alphaHHT';
analysisSettings.alphaHHT   =  0        ;
analysisSettings.stopTolDeltau = 1e-12 ;
analysisSettings.stopTolForces = 1e-12 ;
% ------------------------------------
[matUspendulumCase1, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
% ------------------------------------
%md Case 2: the fluid has an acceeration
otherParams.problemName     = 'AMCase2';
%md AMBool is turned to true to consider fluid acceleration
AMBool = true;
%md Mesh is now a vertical pendulum
%mdThe coordinates conisdering a mesh of two nodes is:
mesh.nodesCoords = [   0   0    l0 ;...
                       l0  0    0  ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 0 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0  1 2 ] ;
% Fluid parameters
rhoFluid = rho/massratio; nuFluid = 1e-6;
AeroBoolmat = false;
%md Initially straight, motion is only driven by the added mass force with no weight
% angle_init = 0;
nameFuncVel = 'windUniform';
%md Drag and lift are ignored in this idealized example
elements(2).aeroCoefs   = {[]; []; [] }   ;
% hydro cross-section props
numGaussPoints  = 4 ;
elements(2).elemTypeAero = [0 0 -d numGaussPoints AeroBoolmat] ;
%md Analysis Settings
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
analysisSettings.geometricNonLinearAero = true;
analysisSettings.booleanSelfWeight = false ;
analysisSettings.stopTolDeltau = 1e-8 ;
analysisSettings.stopTolForces = 1e-8 ;
% ------------------------------------
[matUspendulumCase2, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
% ------------------------------------------------------------------------------
%md### extract control displacements
%mdThe mass displacement in z are:
controlDofDispZ = 6 + 5 ;
controlDispZCase1 = matUspendulumCase1( controlDofDispZ , : ) + (l0-cosd(angle_init)*l0);
controlDispZCase2 = matUspendulumCase2( controlDofDispZ , : ) ;
%mdanalogously the mass dipslacement in x are:
controlDofDispX = 6 + 1 ;
controlDispXCase1 = matUspendulumCase1( controlDofDispX , : ) + sind(angle_init)*l0;
controlDispXCase2 = matUspendulumCase2( controlDofDispX , : ) ;
%
%mdIn order to contrast the solution with the literature refrence the bounce angle measured from the vertical is computed:
angleThetaCase1= rad2deg( atan2( controlDispXCase1, l0 - controlDispZCase1 ) ) ;
angleThetaCase2= rad2deg( atan2( controlDispXCase2, l0 - controlDispZCase2 ) ) ;
%
%mdTo plot diplsacements against $t$ the time vector is:
dt = analysisSettings.deltaT;
times  = (0:length(controlDispZCase1)-1) * dt ;
% Analytical solution for Case 1
d = 2*sqrt(A/pi);
if ~isempty( massratio ) % massratio = rho_structure/rho_fluid
    AMcoef =  1+(1/massratio) ;
else AMcoef = 1;
end
T_ana_lim = 2*pi*sqrt(AMcoef*(d^2/(8*l0) + 2*l0/(3*g)));
f_ana_lim  = 1/T_ana_lim;
theta_ana = angle_init*cos(2*pi*f_ana_lim.*times);
%md Plot angle solution for case 1
figure(), hold on, grid on
plot( times, angleThetaCase1, 'rx')
hold on
plot( times, theta_ana, 'ko')
xlabel('time (s)'), ylabel('\theta (degrees)')
title(sprintf('Pendulum angle, massratio=%d', massratio))
xlabel('time (s)'), ylabel('\theta (degrees)')
title('Angle of the pendulum')
legend('ONSAS', 'analytical')
%md

if isThisOctave
  figure(), hold on, grid on
  plot( times, angleThetaCase2, 'bx')
  xlabel('time (s)'), ylabel('\theta (degrees)')

else
%md Plot angle solution for case 2
figure(), hold on, grid on
yyaxis right
plot( times, angleThetaCase2, 'bx')
xlabel('time (s)'), ylabel('\theta (degrees)')
hold on
%md Plot fluid load and pendulum angle for Case 2
a = times; f = times;
dt = times(2) - times(1);
madded = (1+1)*pi* d^2/4 * l0* rhoFluid/2; % (1+Ca) * Volume * density /2
for t = 1: length(times)-1
    acc = (feval(nameFuncVel, 0, times(t+1)) - feval(nameFuncVel, 0, times(t)))/dt ;
    a(t) = acc(1);
    f(t) = madded * acc(1);
end
yyaxis left
ylabel('fluid load x component');
plot(times(1:end-2), f(1:end-2), 'r-')
title(sprintf('Angle of a pendulum subected only to the added mass force of the swell'))
legend('added mass force', 'pendulum angle (degrees)')
end