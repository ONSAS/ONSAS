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

% 2*10 = rho*A*l0 ;

rho    = 2*10 / ( A * l0 ) ;

nodalDamping = 0.000 ;
%~ nodalMass    = 

% method
timeIncr   =  0.1    ;
finalTime  =  4.13                ;
%~ finalTime  = 0.3                 ;
nLoadSteps = finalTime/timeIncr ;
DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-7           ;
%~ stopTolDeltau = 1e-2           ; 
%~ stopTolForces = 1e-2           ;
stopTolIts    = 30              ;
% ------------------------------------



% --- structural properties ---
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;

secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [ 1  inf  0  inf  0  inf 0  ...
               ];

Nodes = [    0  0  0   ; ...
             l0 0  0 ] ;
             %~ 0  0  -l0 ] ;

Conec = [ 1 2 0 0 1 1 1 ] ; 


% -------------------
nodalConstantLoads   = [ 2  0  0  0  0  -rho*A*l0*0.5*9.8  0 ];
% or
%~ nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
%~ userLoadsFilename = 'myLoadSpringMass' ;

%~ nodalDamping = cres ;
% -------------------

controlDofInfo = [ 2 1 1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
%~ u0    = 0 ;
%~ udot0 = 7.725   ; % m/sec

%~ nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
%~ nonHomogeneousInitialCondUdot0 = [ 2 1 udot0 ] ;

numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%~ plotParamsVector = [2 5 ];
plotParamsVector = [ 3 20 ];
printflag = 0 ;

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;


angs = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;

close all

stopTolForces = 1e0           ;
problemName = 'simplePendulumLooseTol' ;

numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;

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

cd(dirOnsas);
cd(outputDir )
print( [ problemName(1:end-8) '_Bathe'  ] ,'-dpdflatex','-tight') ;
cd(acdir)
