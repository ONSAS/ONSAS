% ------------------------------------------------------------------------------ 
% ------      ONSAS example file: cantilever with nodal moment example    ------
% ------------------------------------------------------------------------------
clear all, close all

dirOnsas = [ pwd '/..' ]; addpath( dirOnsas );
problemName = 'uniformCurvatureCantilever' ;

l = 10   ;    ty = .1 ;  tz = .1 ;  Nelem = 20 ;

Nodes = [ (0:(Nelem))'*l/Nelem zeros(Nelem+1,2) ] ;

auxconec = [ (ones(Nelem,1)*[ 1 2 0 1 0]) (1:(Nelem))' (2:(Nelem+1))' ] ;

Conec = cell(2+Nelem,1) ;

Conec{1, 1} = [ 0 1 0 0 1                     1       ] ; % fixed node
Conec{2, 1} = [ 0 1 1 0 0                     Nelem+1 ] ; % loaded node
for i=1:Nelem
  Conec{2+i, 1} =  auxconec(i,:) ;
end


% ======================================================================
% --- MELCS parameters ---

E = 200e9 ;  nu = 0.3 ;  rho = 0 ;
materialsParams = { [ rho 1 E nu] } ;

elementsParams  = { 1; 3} ;

loadsParams   = {[ 1 1   0 0 0 -1 0 0 ]} ;

% --- cross section ---
A = ty*tz ;     Iy = ty*tz^3/12 ;    Iz = tz*ty^3/12 ;     It = 1 ;
% crossSecsParams = {[ 1 A It Iy Iz ]} ;
crossSecsParams = {[ 2 ty tz ]}      ;

springsParams    = {[ inf  inf  inf  inf  inf  inf ]} ;

storeBoolean    = 1 ;

controlDofs     = [ Nelem+1  4  -1 ] ;

stopTolIts      = 30      ;
stopTolDeltau   = 0       ;
stopTolForces   = 1.0e-10 ;
targetLoadFactr = E * Iy * 2 * pi / l ;  % moment of curvature 2pi/l
nLoadSteps      = 10 ;

plotParamsVector = [ 3 ] ; printFlag = 2 ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 

analyticSolFlag = 1 ; analyticCheckTolerance = 1e-4 ; analyticFunc = @(w) w * l / ( E * Iy ) ;

% --- ONSAS execution ---
ONSAS
return

lw = 3.5; ms = 11; plotfontsize = 22 ;

xs = Nodes(:,1) ;
u  = matUs(1:2:end,11) ;
ux11=u(1:3:end);
uz11=u(3:3:end);

u  = matUs(1:2:end,6) ;
ux6=u(1:3:end);
uz6=u(3:3:end);

figure
plot( xs, 0*xs ,'b-.' , 'linewidth', lw,'markersize',ms )
hold on, grid on
axis equal
plot( xs+ux6, uz6 ,'k' , 'linewidth', lw,'markersize',ms )
plot( xs+ux11, uz11 ,'r' , 'linewidth', lw,'markersize',ms )
%~ plot( controlDispsNRAL, loadFactorsNRAL,'r-s' , 'linewidth', lw,'markersize',ms )
%~ plot( controlDispsNR, loadFactorsNR,'k-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('x (m)');   laby = ylabel('z (m)') ;
%~ legend('analytic','NRAL','NR','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
print( [ 'unif' ] ,'-dpdflatex','-tight') ;
cd(acdir);
