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

nElemsPerBeam = 8 ;

auxRs = linspace(0,L, nElemsPerBeam+1 )' ;
auxRs(1) = [] ;

c1 = cos(2*pi*1/3) ; s1 = sin(2*pi*1/3) ;
c2 = cos(2*pi*2/3) ; s2 = sin(2*pi*2/3) ;

Nodes = [ 0                        0         0        ; ...
          zeros(nElemsPerBeam, 1)  auxRs*1   auxRs*0  ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c1  auxRs*s1 ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c2  auxRs*s2 ] ;

aux = (1:(nElemsPerBeam+1))' ;

Conec = [ aux(1:(end-1)) aux(2:end) zeros(nElemsPerBeam,2) ; ...
           1             aux(2)+nElemsPerBeam 0 0 ; ...
          aux(2:(end-1))+nElemsPerBeam aux(3:end)+nElemsPerBeam zeros(nElemsPerBeam-1,2) ; ...
           1             aux(2)+2*nElemsPerBeam 0 0 ; ...
          aux(2:(end-1))+2*nElemsPerBeam aux(3:end)+2*nElemsPerBeam zeros(nElemsPerBeam-1,2) ] ; 

Conec = [ Conec [(ones(size(Conec,1)/3,1)*[ 1 1 2])   ; ...
                 (ones(size(Conec,1)/3,1)*[ 1 2 2])   ; ...
                 (ones(size(Conec,1)/3,1)*[ 2 2 2]) ] ] ;





J   = I ;

global Jrho
Jrho = rho * diag( [ J; I; I ] ) ;

materialsParams = {[ rho 1 E nu ],[ rho 1 2*E nu ]} ;

crossSecsParams = [   A I I J ;
                    2*A I I J ] ;



% method
timeIncr   =  0.1    ;
finalTime  = 3 ;    
%~ finalTime  = 15 ;    
nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-6        ;
stopTolIts    = 30          ;
% ------------------------------------

nodalSprings = [ 1 inf 0 inf inf inf inf ] ;


% -------------------
nodalVariableLoads   = [ nElemsPerBeam+1  0  0  0  0  1  0 ];

controlDofs = [ 1 2 1 ] ;

%~ loadFactorsFunc = @(t) 5e5*t*(t<1) + (10e5-5e5*t)*(t>=1)*(t<2) + 0 ;
loadFactorsFunc = @(t) 5e5*sin( 2*pi * t / ( finalTime/4 ) ) ;

alphahht    =  -0.05               ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphahht ] ;

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
