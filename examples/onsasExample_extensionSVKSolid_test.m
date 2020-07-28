%% Example uniaxialSolid
% Linear elastic solid submitted to uniaxial loading. 
% Geometry given by $L_x$, $L_y$ and $L_z$, tension $p$ applied on 
% face $x=L_x$.
%
% Analytical solution to be compared with numerical:
% $$ u_x(x=L_x,y,z) = \frac{p L_x}{E} $$
%%

% uncomment to delete variables and close windows
clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'extensionSVKSolidManual' ;

%% Structural properties

% tension applied and x, y, z dimensions
p = 1 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;

% an 8-node mesh is considered with its connectivity matrix
Nodes = [ 0    0    0 ; ...
          0    0   Lz ; ...
          0   Ly   Lz ; ...
          0   Ly    0 ; ...
          Lx   0    0 ; ...
          Lx   0   Lz ; ...
          Lx  Ly   Lz ; ...
          Lx  Ly    0 ] ;

Conec = [ 1 4 2 6 1 1 3 ; ...
          6 2 3 4 1 1 3 ; ...
          4 3 6 7 1 1 3 ; ...
          4 1 5 6 1 1 3 ; ...
          4 6 5 8 1 1 3 ; ...
          4 7 6 8 1 1 3 ] ;

% Material and geometry properties
E = 1 ; nu = 0.3 ;
  
materialsParams = cell(1,1) ;  
materialsParams{1} = [ 0 6 E nu ] ;

crossSecsParams = [ 0 0 0 0 ] ;

% Displacement boundary conditions and springs
nodalSprings = [ 1 inf 0  inf 0   inf 0 ; ...
                 2 inf 0  inf 0   0   0 ; ...
                 3 inf 0  0   0   0   0 ; ...
                 4 inf 0  0   0   inf 0 ; ...
                 5 0   0  inf 0   inf 0 ; ...
                 6 0   0  inf 0   0   0 ; ...
                 8 0   0  0   0   inf 0 ] ;

%% Loading parameters
nodalForce = p * Ly * Lz / 6 ;
nodalVariableLoads = [ (5:8)' nodalForce*[1 2 1 2]' zeros(4,5) ] ;

%% Analysis parameters
stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-12 ;
stopTolForces    = 1.0e-12 ;
targetLoadFactr  = 2    ;

%~ nLoadSteps       = 2    ;
nLoadSteps       = 10    ;

controlDofs = [ 7 1 1 ] ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

% Analytic sol
analyticSolFlag = 2 ;
analyticCheckTolerance = 1e-8 ;
analyticFunc = @(w) E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) )

%% Output parameters
plotParamsVector = [ 0 ] ;
printflag = 2 ;

reportBoolean = 0;

%% run ONSAS
acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

controlDispsCase1 = controlDisps ;
analyticValsCase1 = analyticVals ;
loadFactorsCase1  = loadFactors  ;

close all

% --------------------------------------------------------

problemName = 'extensionSVKGMSHAndComplexStep' ;

consMatFlag = 1 ;

[ nodesMat, conecMat ] = meshFileReader( 'geometry_extensionSVK.msh' ) ;

suppsMat = [ inf 0  0 	0   0 	0 ; ...
             0 	 0  inf 0   0   0 ; ...
             0 	 0  0   0   inf 0 ] ;

% Loads matrix: 		Is defined by the corresponding load label. First entry is a boolean to assign load in Global or Local axis. (Recommendation: Global axis). 
%										Global axis -> 1, local axis -> 0. 
%										The structure of the matrix is: [ 1/0 Fx Mx Fy My Fz Mz ]

%~ loadsMat = [0   0 0 0 0 p 0 ] ;  % --- local loading ---
loadsMat = [1   p 0 0 0 0 0 ] ; % --- global loading ---

[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat ) ;

clear nodesMat conecMat loadsMat suppsMat

plotParamsVector = [ 3 ] ;

% run ONSAS
acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

% --- plots ---

lw = 2.0 ; ms = 10 ; plotfontsize = 22 ;

figure
plot( controlDisps, analyticVals ,'b-o' , 'linewidth', lw,'markersize',ms )
grid on, hold on
plot( controlDispsCase1, loadFactorsCase1  ,'k-s' , 'linewidth', lw,'markersize',ms)
plot( controlDisps, loadFactors  ,'r-x' , 'linewidth', lw,'markersize',ms)

%~ figure
%~ semilogy(controlDisps, abs( analyticVals-loadFactors) )

% ---------------
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic Sol','numerical Sol 1','numerical Sol 2','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
%~ print( [ 'plotsExtensionSVK'  ] ,'-depslatex') ;

cd(dirOnsas); cd(outputDir);
print( [ 'plotsExtensionSVK' ] ,'-dpdflatex','-tight') ;
cd(acdir);
