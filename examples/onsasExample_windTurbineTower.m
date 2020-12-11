% ------------------------------------------------------------------------------
% example: simple model of a wind turbine with tubular blades
% ------------------------------------------------------------------------------

clear all, close all

dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName = 'windTurbine' ;

% ------------------------------------
E         = 200e9 ;
nu        = 0.3   ;
phiExt    = 3     ;
thickness = 0.01  ; 

A   = pi * (phiExt^2/4 - (phiExt-thickness)^2/4 ) ;
I   = pi * (phiExt^4   - (phiExt-thickness)^4)   / 64 ;  J = I ;
L   =  30   ;
rho =  8e3    ;

deepn = 0.1;

nElemsPerBeam = 8 ;

% ----------------------------------------------------------------------

auxRs = linspace(0,L, nElemsPerBeam+1 )'; auxRs(1) = [] ;

c1 = cos(2*pi*1/3) ; s1 = sin(2*pi*1/3) ;
c2 = cos(2*pi*2/3) ; s2 = sin(2*pi*2/3) ;

Nodes = [ 0                        0         0        ; ...
          zeros(nElemsPerBeam, 1)  auxRs*1   auxRs*0  ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c1  auxRs*s1 ; ...
          zeros(nElemsPerBeam, 1)  auxRs*c2  auxRs*s2 ; ...
	  -deepn*L 0 0   ; ...
	  -deepn*L -L 0 ; ...
	  -deepn*L -2*L 0 ; ...
	  -deepn*L -2.5*L 0 ; ...
	  -deepn*L -3*L 0 ...
	  ] ;

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

%~ Conec{3*nElemsPerBeam+1,1} = [ 0 2 0 0 1  1               ] ;
Conec{3*nElemsPerBeam+1,1} = [ 0 2 0 0 1  nElemsPerBeam*3+1+1+4 ] ;

Conec{3*nElemsPerBeam+2,1} = [ 0 2 1 0 0  nElemsPerBeam+1 ] ;

Conec{ 3*nElemsPerBeam+3 , 1 } = [ 1 1 0 2 0   1 3*nElemsPerBeam+2 ] ;

Conec{ 3*nElemsPerBeam+3+1 , 1 } = [ 1 1 0 1 0   3*nElemsPerBeam+2 3*nElemsPerBeam+3] ;
Conec{ 3*nElemsPerBeam+3+2 , 1 } = [ 1 1 0 1 0   3*nElemsPerBeam+3 3*nElemsPerBeam+4] ;
Conec{ 3*nElemsPerBeam+3+3 , 1 } = [ 1 1 0 1 0   3*nElemsPerBeam+4 3*nElemsPerBeam+5] ;
Conec{ 3*nElemsPerBeam+3+4 , 1 } = [ 1 1 0 1 0   3*nElemsPerBeam+5 3*nElemsPerBeam+6] ;

Conec{ 3*nElemsPerBeam+3+5 , 1 } = [ 0 2 2 0 2   1 ] ;
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
loadsParams{2,1} = [ 1 0  -1e6 0  0 0  0 0 ] ;

% --- CrossSection parameters ---
crossSecsParams{1,1} =  [ 2 phiExt*.4 phiExt*.4 ] ; 
crossSecsParams{2,1} =  [ 1 A 0.01*J I I ] ; 
% ----------------------------------------------------------------------
% --- springsAndSupports parameters ---
%~ springsParams{1, 1} = [ inf 1e4  inf inf inf inf ] ;
springsParams{1, 1} = [ inf inf inf inf inf inf ] ;
springsParams{2, 1} = [ 0 0 0 0 inf 0 ] ;


% method
timeIncr   =  0.08    ;
%~ finalTime  = 30*timeIncr ;    
finalTime  = 5 ;    
%~ finalTime  = 15 ;    
nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-6        ;
stopTolIts    = 30          ;
% ------------------------------------

controlDofs = [ 1 2 1 ] ;

loadFactorsFunc = @(t) 1e6*sin( 2*pi * t / ( finalTime ) ) ;

alphahht    =  -0.05               ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphahht ] ;

storeBoolean = 1;

plotParamsVector = [ 3 ] ;
printFlag        =   0   ;

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
