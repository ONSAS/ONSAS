% ------------------------------------------------------------------------------
% example simple wind turbine
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/..' ] ;   problemName = 'windTurbine' ;

% ------------------------------------
E   =  200e9  ;
nu  =  0.3  ;
A   =  .2*.2    ;
I   =  5^4/12 ;
L   =  30   ;
rho =  8e3    ;

J   = I ;

global Jrho
Jrho = rho * diag( [ J; I; I ] ) ;

materialsParams = {[ rho 1 E nu ],[ rho 1 2*E nu ]} ;

crossSecsParams = [   A I I J ;
                    2*A I I J ] ;

% method
timeIncr   =  0.050    ;
finalTime  = 3 ;    
%~ finalTime  = 15 ;    
nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-6        ;
stopTolIts    = 30          ;
% ------------------------------------

nodalSprings = [ 1 inf 0 inf inf inf inf ] ;

nElemsPerBeam = 8 ;

auxRs = linspace(0,L, nElemsPerBeam+1 )' ;

Nodes = [ zeros(nElemsPerBeam+1,1)     auxRs                zeros(nElemsPerBeam+1,1) ; ...
          zeros(nElemsPerBeam  ,1)     auxRs(2:end)*cos(2*pi  /3)  auxRs(2:end)*sin(2*pi  /3) ; ...
          zeros(nElemsPerBeam  ,1)     auxRs(2:end)*cos(2*pi*2/3)  auxRs(2:end)*sin(2*pi*2/3) ] ;

aux = (1:(nElemsPerBeam+1))' ;
Conec = [ aux(1:(end-1)) aux(2:end) zeros(nElemsPerBeam,2) ; ...
           1             aux(2)+nElemsPerBeam 0 0 ; ...
          aux(2:(end-1))+nElemsPerBeam aux(3:end)+nElemsPerBeam zeros(nElemsPerBeam-1,2) ; ...
           1             aux(2)+2*nElemsPerBeam 0 0 ; ...
          aux(2:(end-1))+2*nElemsPerBeam aux(3:end)+2*nElemsPerBeam zeros(nElemsPerBeam-1,2) ] ; 

Conec = [ Conec [(ones(size(Conec,1)/3,1)*[ 1 1 2])   ; ...
                 (ones(size(Conec,1)/3,1)*[ 1 2 2])   ; ...
                 (ones(size(Conec,1)/3,1)*[ 2 2 2]) ] ] ;

% -------------------
nodalVariableLoads   = [ nElemsPerBeam+1  0  0  0  0  1  0 ];

controlDofs = [ 1 2 1 ] ;

loadFactorsFunc = @(t) 5e5*t*(t<1) + (10e5-5e5*t)*(t>=1)*(t<2) + 0 ;

DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;
numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

storeBoolean = 1;

plotParamsVector = [ 3 ]; sectPar = [ 12 1 1 ] ;
printFlag = 0 ;

reportBoolean = 0 ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
grid on
labx = xlabel('Time (s)');   laby = ylabel(sprintf('Displacement node: %2i dof %1i', controlDofs(1), controlDofs(2) ) ) ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print('windturbine','-dpdflatex','-tight')
cd(acdir);

figure
grid on
plot(timesVec, loadFactors,'r','linewidth',lw,'markersize',ms);
