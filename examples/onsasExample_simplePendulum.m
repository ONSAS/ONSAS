% ------------------------------------
% TEST example simple pendulum
% ------------------------------------

clear all, close all

% --- general data ---
inputONSASversion = '0.1.10';
dirOnsas = [ pwd '/..' ] ;
problemName = 'simplePendulumTrussHHT' ;
% ------------------------------------

% -- scalar params -----
Es  = 10e11  ;
nu  = 0      ;
A   = 0.1    ;
l0  = 3.0443 ;
m   = 10     ;
g   = 9.81   ;

nodalDispDamping = 0 ;
rho = 2*m / ( A * l0 ) ;


% ----- geometry -----------------
Nodes = [     0 0  l0   ; ...
             l0 0  l0 ] ;

Conec = { [ 0 1 0 0 1  1   ] ; ...
          [ 0 1 1 0 2  2   ] ; ...
          [ 1 2 0 1 0  1 2 ] } ;


%  --- MELCS params ---

materialsParams = { [ rho 3 Es nu ] } ;

elementsParams = { 1 ; [ 2 0 ] } ;

loadsParams  = {[ 1 0  0  0  0  0  -rho*A*l0*0.5*g  0 ]};

crossSecsParams = { [2 sqrt(A) sqrt(A) ] } ;

springsParams = {[  inf  0  inf  0  inf 0 ] ; ...
                 [  0    0  inf  0  0   0 ] };

% -------------------


% method

timeIncr   =  0.05    ;
%~ timeIncr   =  0.1    ;

finalTime  =  4.13*10                ;
nLoadSteps = finalTime/timeIncr ;


alphaHHT = -0.05 ;
%~ DeltaNW    =  0.5               ;
%~ AlphaNW    =  0.25              ;
%~ alphaHHT = 0 ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForcesTight = 1e-6           ;

stopTolForcesTight = 1e-6           ;
stopTolForcesLoose = 1e+0           ;

stopTolIts    = 30              ;
% ------------------------------------

reportBoolean = 0 ;

controlDofs = [ 2 1 1 ] ;

% analysis parameters

%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForcesTight stopTolIts DeltaNW AlphaNW] ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForcesTight stopTolIts alphaHHT ] ;

%~ plotParamsVector = [ 1 ];
sectPar = [12 .3 .3 ] ;
plotParamsVector = [ 3 40 ];
%~ plotParamsVector = [2 5 ];
%~ plotParamsVector = [ 3 20 ];
printFlag = 0 ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

uNum = PenduloNL_HHT( l0, A, Es, m, nodalDispDamping, g, timeIncr, -alphaHHT, finalTime, stopTolForcesTight, stopTolIts, 1, 0 ) ;

figure
hold on, grid on
plot(timesVec, controlDisps, 'b-o')
plot(timesVec, uNum(1,:)-l0 ,'r-x')
ylabel('control displacement')
xlabel('time (s)')
legend('onsas','semi-analytic')
print(['../../salida_dt_' sprintf('%05.3f',timeIncr) '.png'],'-dpng')
%~ stop

%~ angs2 = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;

%~ lw  = 2   ; ms  = 5.5 ;
%~ lw2 = 3.2 ; ms2 = 23 ;
%~ plotfontsize = 22 ;

%~ figure
%~ plot(timesVec, angs,'b--o','linewidth',lw,'markersize',ms);
%~ grid on, hold on
%~ plot(timesVec, angs2,'r-x','linewidth',lw,'markersize',ms);

%~ labx=xlabel('time (s)'); laby=ylabel('Angle (degrees)'); labz=zlabel('z') ;
%~ set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize*0.8 ) ;
%~ set(tit, "FontSize", plotfontsize) ;
%~ set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize) ;
%~ legend('tight tolerances','loose tolerances','location','southeast')

%~ cd(dirOnsas); cd(outputDir ); 
%~ print( [ problemName(1:end-8) '_Bathe'  ] ,'-dpdflatex','-tight') ;
%~ cd(acdir)






%~ angs = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;

%~ close all

%~ problemName = 'simplePendulumLooseTol' ;

%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForcesLoose stopTolIts DeltaNW AlphaNW] ;
%~ numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForcesLoose stopTolIts alphaHHT ] ;

%~ acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;
