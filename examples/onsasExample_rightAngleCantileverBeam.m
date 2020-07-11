% ------------------------------------------------------------------------------
% example Right-angle cantilever
%
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/..' ] ;   problemName = 'rightAngleCantileverBeam' ;

% ------------------------------------
E   =  1e6  ;
nu  = -0.5  ;
A   =  1    ;
I   =  1e-3 ;
L   =  10   ;
rho =  1    ;

% inconsistent
J   = I ;
global Jrho
Jrho = diag( [ 20; 10; 10 ] ) ;

materialsParams = {[ rho 1 E nu ]} ;

crossSecsParams = [ A I I J ] ;

% method
timeIncr   =  0.050    ;
finalTime  = .5    ;
%~ finalTime  = 15 ;    
nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-7        ;
stopTolIts    = 30          ;
% ------------------------------------

nodalSprings = [ 1 inf inf inf inf inf inf ] ;

%~ nElemsPerBeam = 4 ;
nElemsPerBeam = 10 ;
auxCoords = linspace(0,L,nElemsPerBeam+1)' ;

Nodes = [ zeros(nElemsPerBeam+1,1)       auxCoords  zeros(nElemsPerBeam+1,1) ; ...
          -auxCoords(2:end) ones(nElemsPerBeam,1)*L       zeros(nElemsPerBeam,1) ] ;

aux = (1:(2*nElemsPerBeam+1))' ;
Conec = [ aux(1:(end-1)) aux(2:end) zeros(2*nElemsPerBeam,2) ...
          ones(2*nElemsPerBeam,2) 2*ones(2*nElemsPerBeam, 1) ] ; 

% -------------------
nodalVariableLoads   = [ nElemsPerBeam+1  0  0  0  0  1  0 ];

controlDofs = [ nElemsPerBeam+1 3  1 ] ;

loadFactorsFunc = @(t) 50*t*(t<1) + (100-50*t)*(t>=1)*(t<2) + 0 ;
DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;
numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%~ alphaHHT = 0 ;
%~ alphaHHT = -0.05 ;
%~ numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

%~ plotParamsVector = [0 ];
plotParamsVector = [ 3 ]; sectPar = [ 12 .25 .25 ] ;
printFlag = 0 ;

reportBoolean = 0 ;

storeBoolean = 0 ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
%~ plot(controlDisps,'b--o','linewidth',lw,'markersize',ms);
plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
grid on
labx = xlabel('Time (s)');   laby = ylabel(sprintf('Displacement node: %2i dof %1i', controlDofs(1), controlDofs(2) ) ) ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print('rightAngle','-dpdflatex','-tight')
%~ print('rightAngle','-dpdflatex')
%~ print('rightAngle.png','-dpng')
cd(acdir);


