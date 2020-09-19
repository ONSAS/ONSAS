% ------------------------------------------------------------------------------
% example: simple model of a wind turbine with tubular blades
% ------------------------------------------------------------------------------

clear all, close all

problemName = 'windTurbine' ;

% ------------------------------------
E         = 200e9 ;
nu        = 0.3   ;
phiExt    = 1     ;
thickness = 0.005  ; 

A   = pi * (phiExt^2/4 - (phiExt-thickness)^2/4 ) ;
I   = pi * (phiExt^4   - (phiExt-thickness)^4)   / 64 ;  J = I ;
L   =  30   ;
rho =  8e3    ;

nElemsPerBeam = 8 ;

% ----------------------------------------------------------------------

auxRs = linspace(0,L, nElemsPerBeam+1 )'; auxRs(1) = [] ;

c1 = cos(2*pi*1/3) ; s1 = sin(2*pi*1/3) ;
c2 = cos(2*pi*2/3) ; s2 = sin(2*pi*2/3) ;

Nodes = [ 0                        0         0        ; ...
          zeros(nElemsPerBeam, 1)  auxRs*1   auxRs*0  ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c1  auxRs*s1 ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c2  auxRs*s2 ] ;

Conec = cell(nElemsPerBeam*3,1) ;

for j=1:3
  i = 1 ;
  Conec{   1+(j-1)*nElemsPerBeam , 1 } = [ 1 1 0 1 0   1  ... 
                                                         2+(j-1)*nElemsPerBeam ] ;
  for i=2:nElemsPerBeam
    Conec{ i+(j-1)*nElemsPerBeam , 1 } = [ 1 1 0 1 0   i+(j-1)*nElemsPerBeam ...
                                                         i+(j-1)*nElemsPerBeam+1 ] ;
  end
end

Conec{3*nElemsPerBeam+1,1} = [ 0 2 0 0 1  1               ] ;
Conec{3*nElemsPerBeam+2,1} = [ 0 2 1 0 0  nElemsPerBeam+1 ] ;
% ----------------------------------------------------------------------

% ======================================================================
% --- MELCS parameters ---

materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

% --- Material parameters ---
materialsParams{1,1} = [ rho 2 E nu ] ;

% --- Element parameters ---
elementsParams{1,1} = [ 3 ] ;
elementsParams{2,1} = [ 1 ] ;

% --- Load parameters ---
loadsParams{1,1} = [ 1 1  0 0  0 0  1 0 ] ;

% --- CrossSection parameters ---
crossSecsParams{1,1} =  [ A I I J ] ; 
% ----------------------------------------------------------------------
% --- springsAndSupports parameters ---
springsParams{1, 1} = [ inf 1e4  inf inf inf inf ] ;


% method
timeIncr   =  0.1    ;
%~ finalTime  = 30*timeIncr ;    
finalTime  = 10 ;    
%~ finalTime  = 15 ;    
nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-6        ;
stopTolIts    = 30          ;
% ------------------------------------

controlDofs = [ 1 2 1 ] ;

loadFactorsFunc = @(t) 1e5*sin( 2*pi * t / ( finalTime ) ) ;

alphahht    =  -0.05               ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphahht ] ;

storeBoolean = 1;

plotParamsVector = [ 3 ]; sectPar = [ 12 1 1 ] ;
printFlag = 0 ;

reportBoolean = 0 ;



% --------------------------------------
run( [ pwd '/../ONSAS.m' ] ) ;
% --------------------------------------

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
