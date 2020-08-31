% -----------------------------
% deployable ring example
% -----------------------------

clear all, close all

E  = 200e3 ;   nu = 0.3   ;  R = 120   ;

dirOnsas = [ pwd '/..' ] ;
inputONSASversion = '0.1.10';
problemName = 'deployableRing' ;

materialsParams = cell(1,1) ;
materialsParams{1} = [0 1 E nu] ;

b = .6;   h = 6 ;

sectPar = [12 b h ]; 

A  = b*h      ; It = h*b^3/3 ;
Iy = b*h^3/12 ; Iz = h*b^3/12 ;

crossSecsParams = [ A Iy Iz It ] ;

% Defino los apoyos de la estructura

%~ Nelem = 2*32 ;  Np = Nelem +1;
%~ Nelem = 2*5 ;  Np = Nelem +1;
%~ Nelem = 2*20 ;  Np = Nelem +1;
Nelem = 2*22 ;  Np = Nelem +1;

al = 2*pi/Nelem ;
Nodes = zeros(Np,3);
for i=1:(Np)
  Nodes(i,1) = R * cos( pi+ (i-1)*al );
  Nodes(i,2) = R * sin( pi+ (i-1)*al );
end

nodalSprings = [ 1         inf  inf  inf  inf  inf  inf ; ...
                 Nelem/2+1   0    0  inf  inf  inf  inf ] ;

Conec = [ (1:(Nelem))' (2:(Nelem+1))' zeros(Nelem,2)  (ones(Nelem,1)*[ 1 1 2]) ] ;
Conec ( end,2) = 1;
Nodes(end,:)=[];

nodalVariableLoads   = [ Nelem/2+1  0 1  0 0 0 0 ] ;

controlDofs = [ Nelem/2+1 2 1 ] ;

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-10 ;
stopTolForces  = 1.0e-10 ;
targetLoadFactr = 3*8e2 ;
%~ nLoadSteps      = 180 ; incremArcLen     = .5    ;
nLoadSteps      = 3000 ; incremArcLen     = .8    ;

%~ targetLoadFactr = 3*8e1 ;
%~ nLoadSteps      = 15 ; incremArcLen     = .75    ;

plotParamsVector = [ 3 50 ] ;
%~ plotParamsVector = [ 2  4 ] ;
plotsViewAxis = [ 2 -1 1];     printflag = 0;

numericalMethodParams = [ 2 ...
 stopTolDeltau stopTolForces stopTolIts ...
 targetLoadFactr nLoadSteps  incremArcLen ] ; 


acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;

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
