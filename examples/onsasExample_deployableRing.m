% -----------------------------
% deployable ring example
% -----------------------------

clear all, close all

% --- scalar parameters ---
E  = 200e3 ; nu = 0.3 ;  R = 120 ;  wy = .6 ;   wz = 6 ;  halfNElem = 10 ;

problemName = 'deployableRing' ;

% --- MELCS ---

materialsParams = { [ 0 1 E nu ] } ;

elementsParams  = { 1; 3 } ;

loadsParams     = { [ 1 1  0 1 0 0 0 0 ]} ;

A  = wy*wz      ;
It = wz*wy^3/3 ; % approximation. improve computation in https://github.com/ONSAS/ONSAS/issues/134

Iy = wy*wz^3/12 ; Iz = wz*wy^3/12 ;

crossSecsParams = {[ 1 A It Iy Iz ]} ;

springsParams   = {[ inf  inf  inf  inf  inf  inf ] ; ...
                   [ 0    0    inf  inf  inf  inf ] } ;

%~ Nelem = 2*32 ;
%~ Nelem = 2*5 ; 
%~ Nelem = 2*20 ;
Nelem = 2*halfNElem ;

Np = Nelem ; % number points (the last is the same as the first)

% --- geometry ---
deltaAngle = 2 * pi / Nelem ;
angs       = linspace( 0, 2*pi-deltaAngle , Nelem )' ;

Nodes = [ R*cos( angs )  R*sin( angs ) zeros(size(angs)) ] ;

Conec = {[ 0 1 0 0 1   1          ] ; ...
         [ 0 1 1 0 2   Nelem/2+1  ] } ;
for i=1:(Nelem-1)
  Conec{i+2, 1} = [ 1 2 0 1 0   i i+1 ] ;
end
Conec{Nelem+2, 1} = [ 1 2 0 1 0   Nelem 1 ] ;


controlDofs = [ Nelem/2+1 2 1 ] ;

% --- analysis parameters ---
stopTolIts      = 30     ;
stopTolDeltau   = 1.0e-10 ;
stopTolForces   = 1.0e-10 ;
targetLoadFactr = 3*8e2  ;
nLoadSteps      = 180 ; incremArcLen     = .5    ;
%~ nLoadSteps      = 3000 ; incremArcLen     = .8    ;

plotParamsVector = [ 3 50 ] ;

numericalMethodParams = [ 2 ...
 stopTolDeltau stopTolForces stopTolIts ...
 targetLoadFactr nLoadSteps  incremArcLen ] ; 

%~ numericalMethodParams = [ 1 ...
 %~ stopTolDeltau stopTolForces stopTolIts ...
 %~ targetLoadFactr nLoadSteps   ] ; 

% --- ONSAS execution ---
run( [  pwd  '/../ONSAS.m' ] ) ;

return
lw = 3.0 ; ms = 11 ; plotfontsize = 22 ;

figure
plot( controlDisps , loadFactors, 'b-.' , 'linewidth', lw,'markersize',ms )
grid on

labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
%~ legend('analytic','NRAL','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print( [ 'ring' ] ,'-dpdflatex','-tight') ;
cd(acdir);
