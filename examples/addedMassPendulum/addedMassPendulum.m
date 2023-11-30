% nonlinear added mass pendulum eample

close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );

%md Case 1: the fluid is still: it is only modelled by an added inertia of the solid,
%md the fluid properties are only defined by the massratio = rho_structure/rho_fluid
otherParams.problemName     = 'addedMassPedulum_case1';

% input scalar parameters
EA = 1e8  ;  nu = 0       ; 
A  = 0.1  ;  l0  = 3.0443 ;
m  = 10   ;  g = 9.80     ;
T  = 4.13 ;

% computed scalar parameters
E   = EA / A ;              % young modulus
d   = 2 * sqrt(A/pi);       % cross-section diameter
rho_structure = 2 * m / ( A * l0 ) ;  % solid density

% materials
materials.modelName   = 'elastic-rotEngStr' ;
materials.modelParams = [ E nu ] ;
materials.density     = rho_structure ;
%md

%md### elements
%md
elements(1).elemType = 'node' ;

elements(2).elemType = 'frame';
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [ d ] ;
elements(2).massMatType  = 'consistent';

%md### boundaryConds

boundaryConds(1).imposDispDofs = [ 1 2 3 5 6] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0] ;

%md### initial Conditions
initialConds = {} ;

%md### mesh parameters
%mdThe coordinates considering a mesh of two nodes is:
angle_init = 25; % degrees
mesh.nodesCoords = [  0                    0  l0                       ; ...
                      sind(angle_init)*l0  0  l0*(1-cosd(angle_init))  ] ;

mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1  1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0  1 2 ] ;

%md### analysisSettings
analysisSettings.deltaT        = 0.01  ;
analysisSettings.finalTime     = 2*T *.1 ;
analysisSettings.stopTolDeltau = 1e-13 ;
analysisSettings.stopTolForces = 1e-13 ;
analysisSettings.stopTolIts    = 30    ;

analysisSettings.booleanSelfWeight = true ;

analysisSettings.methodName = 'alphaHHT';
analysisSettings.alphaHHT   =  0        ;
analysisSettings.stopTolDeltau = 1e-12 ;
analysisSettings.stopTolForces = 1e-12 ;

otherParams.plots_format       = 'vtk' ;

[matUspendulumCase1, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

# stop
# % =========================================================================


# otherParams.problemName     = 'addedMassPedulum_case2';

# analysisSettings.deltaT        = 0.05  ;

# %md AMBool is turned to true to consider fluid acceleration


# %md Mesh is now a vertical pendulum
# %mdThe coordinates conisdering a mesh of two nodes is:
# # mesh.nodesCoords = [   0   0    l0 ;...
# #                        l0  0    0  ] ;
# % Fluid parameters

massratio = inf;
rhoFluid = rho_structure/massratio; nuFluid = 1e-6;

# %md Initially straight, motion is only driven by the added mass force with no weight
# % angle_init = 0;
# nameFuncVel = 'windUniform';

# %md Drag and lift are ignored in this idealized example

# %md Analysis Settings
# analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
# analysisSettings.addedMassBool = true  ;
# analysisSettings.booleanSelfWeight = false ;

# % ------------------------------------
# [matUspendulumCase2, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


% ------------------------------------------------------------------------------
%md### extract control displacements
%mdThe mass displacement in z are:
controlDofDispZ = 6 + 5 ;
controlDispZCase1 = matUspendulumCase1( controlDofDispZ , : ) + (l0-cosd(angle_init)*l0);
# controlDispZCase2 = matUspendulumCase2( controlDofDispZ , : ) ;
controlDofDispX = 6 + 1 ;
controlDispXCase1 = matUspendulumCase1( controlDofDispX , : ) + sind(angle_init)*l0;
# controlDispXCase2 = matUspendulumCase2( controlDofDispX , : ) ;
%
%mdIn order to contrast the solution with the literature refrence the bounce angle measured from the vertical is computed:
angleThetaCase1= rad2deg( atan2( controlDispXCase1, l0 - controlDispZCase1 ) ) ;
# angleThetaCase2= rad2deg( atan2( controlDispXCase2, l0 - controlDispZCase2 ) ) ;
%
%mdTo plot diplsacements against $t$ the time vector is:
dt = analysisSettings.deltaT;
times  = (0:length(controlDispZCase1)-1) * dt ;
% Analytical solution for Case 1
d = 2*sqrt(A/pi);

AMcoef =  1+(1/massratio) ;

T_ana_lim = 2*pi*sqrt(AMcoef*(d^2/(8*l0) + 2*l0/(3*g)));
f_ana_lim  = 1/T_ana_lim;
theta_ana = angle_init*cos(2*pi*f_ana_lim.*times);
%md Plot angle solution for case 1
figure(), hold on, grid on
plot( times, angleThetaCase1, 'rx')
plot( times, theta_ana, 'ko')
xlabel('time (s)'), ylabel('Pendulum angle \theta (degrees)')
title(['AddedMassPendulum - case 1' sprintf(' massratio=%d', massratio)] )
legend('ONSAS', 'analytical')

stop
%md
if isThisOctave
  figure(), hold on, grid on
  plot( times, angleThetaCase2, 'bx')
  xlabel('time (s)'), ylabel('\theta (degrees)')
else
  % Plot angle solution for case 2
  figure(), hold on, grid on
  yyaxis right
  plot( times, angleThetaCase2, 'bx')
  xlabel('time (s)'), ylabel('\theta (degrees)')
  hold on
  % Plot fluid load and pendulum angle for Case 2
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