% ------------------------------------
% TEST example pendulum
% ------------------------------------

clear all, close all

% --- general data ---
inputONSASversion = '0.1.10';

dirOnsas = [ pwd '/..' ] ;
problemName = 'simplePendulumTightTol' ;
% ------------------------------------

Es = 10e11 ;
nu      =  0    ;

A  = 0.1 ;
l0 = 3.0443     ;

rho    = 2*10 / ( A * l0 ) ;

nodalDispDamping = 0.000 ;

booleanConsistentMassMat = 0 ;

% method
timeIncr   =  0.1    ;
finalTime  =  4.13                ;
nLoadSteps = finalTime/timeIncr ;


DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;
%~ alphaHHT = 0 ;
alphaHHT = -0.05 ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForcesTight = 1e-6           ;

stopTolForcesTight = 1e-6           ;
stopTolForcesLoose = 1e+0           ;

stopTolIts    = 30              ;
% ------------------------------------


% --- structural properties ---
materialsParams = {[rho 1 Es nu ]} ;

crossSecsParams = [ A 10 10 10 ] ;

global Jrho
Jrho = diag( [ 20; 10; 10 ] ) ;


nodalSprings = [ 1  inf  0  inf  0  inf 0  ...
               ];

Nodes = [    0  0  0   ; ...
             l0 0  0 ] ;

Conec = [ 1 2 0 0 1 1 1 ] ; 

% -------------------
nodalConstantLoads   = [ 2  0  0  0  0  -rho*A*l0*0.5*9.8  0 ];
% -------------------

controlDofs = [ 2 1 1 ] ;

% analysis parameters

%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForcesTight stopTolIts DeltaNW AlphaNW] ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForcesTight stopTolIts alphaHHT ] ;

%~ plotParamsVector = [ 1 ];
sectPar = [12 .3 .3 ] 
plotParamsVector = [ 3 40 ];
%~ plotParamsVector = [2 5 ];
%~ plotParamsVector = [ 3 20 ];
printFlag = 0 ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

angs = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;

close all

problemName = 'simplePendulumLooseTol' ;

%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForcesLoose stopTolIts DeltaNW AlphaNW] ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForcesLoose stopTolIts alphaHHT ] ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;


angs2 = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;

lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
plot(timesVec, angs,'b--o','linewidth',lw,'markersize',ms);
grid on, hold on
plot(timesVec, angs2,'r-x','linewidth',lw,'markersize',ms);

labx=xlabel('time (s)'); laby=ylabel('Angle (degrees)'); labz=zlabel('z') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize*0.8 ) ;
%~ set(tit, "FontSize", plotfontsize) ;
set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize) ;
legend('tight tolerances','loose tolerances','location','southeast')

cd(dirOnsas); cd(outputDir ); 
print( [ problemName(1:end-8) '_Bathe'  ] ,'-dpdflatex','-tight') ;
cd(acdir)
