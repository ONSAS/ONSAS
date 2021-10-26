
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );

E   = 210e9  ; nu  = 0      ;
A   = .7854/100^2 ; %fi 1cm
l0  = 2 ;
m   = 214     ; g   = 9.81   ;

rho = 2*m / ( A * l0 ) ;

materials.hyperElasModel  = 'SVK' ;
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasParams = [ lambda mu ] ;
materials.density = rho ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemTypeGeometry = [2 sqrt(A) sqrt(A) ] ;
elements(2).elemTypeParams = 1 ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 1.0     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -m*g 0 ] ;

initialConds                = struct() ;



analysisSettings.deltaT        =    0.01  ;
analysisSettings.finalTime      =   3 ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   30   ;

analysisSettings.methodName    = 'newmark' ;
analysisSettings.alphaNM      =   0.25   ;
analysisSettings.deltaNM      =   0.5   ;

%analysisSettings.methodName    = 'alphaHHT' ;
%analysisSettings.alphaHHT      =   0   ;
%alphaHHT = -0.05 ;


otherParams.problemName = 'nonlinearPendulum';
otherParams.plotsFormat = 'vtk' ;
otherParams.nodalDispDamping = 20 ;


mesh.nodesCoords = [   0  0   l0 ; ...
                      l0  0  l0  ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0  1 2 ] ;

% ------------------------------------
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

controlDof = 6 + 5 ;
controlDispsA = matUs( controlDof , : ) ;
timesVec = (0:length(controlDispsA)-1)*analysisSettings.deltaT ;

figure, hold on, grid on, MS = 10; LW = 1.5 ;
plot( timesVec, controlDispsA, 'b-o','markersize',MS,'linewidth',LW)
xlabel('time (s)'), ylabel('control displacement')
%legend('onsasA', 'onsasB', 'semi-analytic')
print('output.png','-dpng')
