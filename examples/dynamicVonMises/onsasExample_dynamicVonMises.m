%md
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar auxiliar parameters are defined.
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters

%%% Parametros Estructura
rho = 7850; % kg/m3 (acero)
Lx  = .374/2;
L0  = .205 ;

Lc  = .240; %m
Ic  = .0254*.0032^3/12; %m4
Ac  = .0254*.0032; %m2
E   = 200000e6 ; %Pa (acero)
kc  = 3*E*Ic/Lc^3; %N/m
mConc   = 1.4; %kg Pandeo incipiente en 1.4
c   = 0; %kg/s (amortiguamiento por friccion juntas y arrastre pesa)
g   = 9.81; %m/s2

tf  = 1.0;
dt  = .000025; % sec

[u, normalForce, times ] = centralDiffDynVonMises(rho, Lx, L0, Lc, Ic, Ac, E, kc, mConc, c, g, tf, dt);

mb = L0*Ac*rho; %kg
Lz = sqrt( L0^2 - Lx^2 ); %m

rhoBarraMasa = mConc*.5 / (Lz*Ac);

M = [mb 0 ; 0 (mb+mConc)/2] ;


nu = .3 ;
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E nu ] ;
materials(1).density =  rho ;

materials(2).hyperElasModel  = '1DrotEngStrain' ;
materials(2).hyperElasParams = [ 1e10*E nu ] ;
materials(2).density = rhoBarraMasa ;

materials(3).hyperElasModel  = '1DrotEngStrain' ;
materials(3).hyperElasParams = [ kc*L0/Ac nu ] ;
materials(3).density = rho ;


elements(1).elemType = 'node' ;

elements(2).elemType = 'truss';
elements(2).elemCrossSecParams{1,1} = 'rectangle'            ;
elements(2).elemCrossSecParams{2,1} = [sqrt(Ac) sqrt(Ac) ] ;
elements(2).massMatType = 'lumped' ;

boundaryConds(1).imposDispDofs = [ 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 ] ;

boundaryConds(2).imposDispDofs =  [ 1 3 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'                  ;
boundaryConds(2).loadsTimeFact = @(t) 1                    ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -(mConc+mb)/2*g 0 ] ;

boundaryConds(3).imposDispDofs =  [ 1 3 ] ;
boundaryConds(3).imposDispVals =  [ 0 0 ] ;

boundaryConds(4).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(4).imposDispVals =  [ 0 0 0 ] ;

initialConds                = struct() ;

mesh.nodesCoords = [   0  0   0   ; ...
                      Lx  0  Lz   ; ...
                      Lx  0  Lz-Lz; ...
                      -L0  0  0    ] ;

mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0   2   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 3 0   3   ] ;
mesh.conecCell{ 4, 1 } = [ 0 1 4 0   4   ] ;

mesh.conecCell{ 4+1, 1 } = [ 1 2 0 0   1 2 ] ;
mesh.conecCell{ 4+2, 1 } = [ 2 2 0 0   2 3 ] ;
mesh.conecCell{ 4+3, 1 } = [ 3 2 0 0   1 4 ] ;

analysisSettings.methodName    = 'newmark' ;
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   0.0005  ;
analysisSettings.finalTime      =   1    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;

%md
%md### otherParams
otherParams.problemName = 'dynamicVonMisesTruss';
otherParams.plotsFormat = 'vtk' ;
%md
%md### Analysis case 1: NR with Rotated Eng Strain
%md In the first case ONSAS is run and the solution at the dof of interest is stored.
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

uONSAS = [ matUs(1,:) ; matUs(6+5,:) ] ;

timesONSAS    = 0:analysisSettings.deltaT:analysisSettings.finalTime ;


subplot(3,1,1)
plot(times(1:10:end),1000*u(1,1:10:end))
hold on
plot(timesONSAS,1000*uONSAS(1,:),'r' )
xlabel('t [s]'); ylabel('u_1 [mm]');
axis( [0 2 1e3*min(u(1,:))*1.1 1e3*max(u(1,:))*1.1] );
subplot(3,1,2)
hold on
plot(times(1:10:end),1000*u(2,1:10:end))
plot(timesONSAS,1000*uONSAS(2,:),'r' )
xlabel('t [s]'); ylabel('u_2 [mm]');
axis([0 2 1e3*min(u(2,:))*1.1 1e3*max(u(2,:))*1.1]);
subplot(3,1,3)
plot(times(1:10:end),normalForce(1:10:end))
xlabel('t [s]'); ylabel('Directa [N]')
axis([0 2 min(normalForce)*1.1 max(normalForce)*1.1]);
