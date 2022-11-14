% ONSAS pnonlinear pendulum

close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
global massratio;
% material scalar parameters
EA = 1e8     ; A = 0.1  ; d = 2*sqrt(A/pi);
E   = EA / A ; nu = 0   ;
l0  = 3.0443 ; T = 4.13 ;
m   = 10     ; g = 9.80 ;
%md
rho = 2*m / ( A * l0 )  ;
materials.density = rho ;
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
elements(2).elemCrossSecParams{1,1} = 'circle' ;% 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ d ] ; 
%mdand the according to the literature example the element include conssitent mass matrix
elements(2).massMatType  = 'consistent';
%md
%md### boundaryConds
%md
%mdThe elements are submitted to two different BC settings. The first BC corresponds to a simple fixed condition (all 3 dipslacments dofs of the node are equal to zero during all time span) then the `boundaryConds(1)` set is:
boundaryConds(1).imposDispDofs = [ 1 2 3 5 6] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0] ;
%md### initial Conditions
%mdNon homogeneous initial conditions are not considered in this example, consequently the `initialConds`  struct is set empty:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates conisdering a mesh of two nodes is:
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
otherParams.problemName     = 'AMVal_nonlinearPendulum';
% ------------------------------------
[matUspnedulum, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md### extract control displacements
%mdThe mass displacement in z are:
controlDofDispZ = 6 + 5 ;
controlDispZ = matUspnedulum( controlDofDispZ , : ) + (l0-cosd(angle_init)*l0);
%mdanalogously the mass dipslacement in x are:
controlDofDispX = 6 + 1 ;
controlDispX = matUspnedulum( controlDofDispX , : ) + sind(angle_init)*l0; 

%mdIn order to contrast the solution with the literature refrence the bounce angle measured from the vertical is computed:
angleTheta= rad2deg( atan2( controlDispX, l0 - controlDispZ ) ) ;

%mdTo plot diplsacements against $t$ the time vector is:
dt = analysisSettings.deltaT;
times  = (0:length(controlDispZ)-1) * dt ;
% Analytical solution
d = 2*sqrt(A/pi);
if ~isempty( massratio ) % massratio = rho_structure/rho_fluid
    AMcoef =  1+(1/massratio) ;
else AMcoef = 1;
end
T_ana_lim = 2*pi*sqrt(AMcoef*(d^2/(8*l0) + 2*l0/(3*g))); 
f_ana_lim  = 1/T_ana_lim;
theta_ana = angle_init*cos(2*pi*f_ana_lim.*times);
%md Plot angle solution
figure(), hold on, grid on
plot( times, angleTheta, 'rx')
hold on
plot( times, theta_ana, 'ko')
xlabel('time (s)'), ylabel('\theta(º)')
title(sprintf('Pendulum angle, massratio=%d', massratio))
%title(sprintf('Angle of the pendulum, 1 +1/massratio=%d', AMcoef))
legend('ONSAS', 'analytical')

fftsig (angleTheta, dt)
fftsig (theta_ana, dt)
title(sprintf('FFT of the pendulum angle, massratio=%d', massratio))
legend('ONSAS', 'analytical')
%%  FFT and signal     
function fftsig (xdefNumlast, dt)
    Ns = length(xdefNumlast);
    xhat = fft(xdefNumlast(1:end),Ns); %same as xhat = fft(zmid(1:end));
    PSD = xhat.*conj(xhat)/Ns ;
    freq = 1/(dt*Ns)*(0:Ns);
    L = 1:floor(Ns/10); % Only plot 1/10th of frequencies
    figure(20)
    plot(freq(L), PSD(L))
    title('FFT'); xlabel('f(Hz)')
    hold on
end
