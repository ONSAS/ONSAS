% ------------------------------------------------------------------------------ 
% ------      ONSAS example file: cantilever with nodal moment example    ------
% ------------------------------------------------------------------------------
clear all, close all

dirOnsas = [ pwd '/..' ] ;      problemName = 'cantileverNodalMoment' ;

l = 10   ;    b = .1 ;  h = .2 ;  Nelem = 20 ;

E = 200e9 ;  nu = 0.3 ;  rho = 0 ;
materialsParams = { [ rho 1 E nu] } ;

% --- cross section ---
A = b*h ;     Iy = b*h^3/12 ;    Iz = h*b^3/12 ;     It = 1 ;

sectPar = [ b h ]; 

storeBoolean = 1 ;

crossSecsParams = [ A Iy Iz It ] ;

nodalSprings    = [ 1  inf  inf  inf  inf  inf  inf ] ;

Nodes = [ (0:(Nelem))'*l/Nelem zeros(Nelem+1,2) ] ;

Conec = [ (1:(Nelem))' (2:(Nelem+1))'  zeros(Nelem,2) (ones(Nelem,1)*[ 1 1 2]) ] ;

%~ booleanCSTangs = 1 ;

nodalVariableLoads   = [ Nelem+1  0 0 0 -1 0 0 ] ;

controlDofs = [ Nelem+1  4  -1 ] ;

stopTolIts      = 30      ;
stopTolDeltau   = 0       ;
stopTolForces   = 1.0e-10 ;
targetLoadFactr = E * Iy / ( l / ( 2 * pi ) ) ;  % curvradius corresponding to perimeter = l
%~ nLoadSteps      = 2 ;
nLoadSteps      = 10 ;

%~ plotParamsVector = [ 2 5 ] ;    plotsViewAxis = [ 0 -1 0 ] ;
plotParamsVector = [ 3 ] ; sectPar = [ 12 b h ] ;
printFlag = 2 ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 

analyticSolFlag = 1 ; analyticCheckTolerance = 1e-4 ; analyticFunc = @(w) w * l / ( E * Iy ) ;

% --- ONSAS execution ---
acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;



lw = 3.5 ; ms = 11 ; plotfontsize = 22 ;

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

