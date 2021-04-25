%% Example uniaxialSolid
% Elastic solid submitted to uniaxial loading. 
% Geometry given by $L_x$, $L_y$ and $L_z$, tension $p$ applied on 
% face $x=L_x$.

clear all, close all

%% set ONSAS.m directory
dirOnsas = [ pwd '/../../src' ] ; addpath( dirOnsas );
otherParams.problemName = 'uniaxialExtension_Manual' ;

%% Structural properties
E = 1 ; nu = 0.3 ;
p = 3 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;


lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasModel = {'SVK'} ;
materials.hyperElasParams = { [ lambda mu ] } ;


elements.elemType = { 'triangle', 'tetrahedron' } ;
elements.elemTypeParams = { [];[] } ;
elements.elemTypeGeometry = { [];[] } ;

boundaryConds.loadsCoordSys = {'global'; [] ; [] ; [] } ;
boundaryConds.loadsTimeFact = { @(t) t ; [] ; [] ; []} ;
boundaryConds.loadsBaseVals = { [p 0 0 0 0 0 ] ; [] ; [] ; [] } ;
boundaryConds.imposDispDofs = { [] ; [1] ; [3] ; [5] } ;
boundaryConds.imposDispVals = { [] ; [0] ; [0] ; [0] } ;


initialConds = struct();

% tension applied and x, y, z dimensions

% an 8-node mesh is considered with its connectivity matrix
mesh.nodesCoords = [ 0    0    0 ; ...
                     0    0   Lz ; ...
                     0   Ly   Lz ; ...
                     0   Ly    0 ; ...
                     Lx   0    0 ; ...
                     Lx   0   Lz ; ...
                     Lx  Ly   Lz ; ...
                     Lx  Ly    0 ] ;

mesh.conecCell = {[ 0 1 1 0    5 8 6   ]; ... % loaded face
                  [ 0 1 1 0    6 8 7   ]; ... % loaded face
                  [ 0 1 2 0    4 1 2   ]; ... % x=0 supp face
                  [ 0 1 2 0    4 2 3   ]; ... % x=0 supp face
                  [ 0 1 3 0    6 2 1   ]; ... % y=0 supp face
                  [ 0 1 3 0    6 1 5   ]; ... % y=0 supp face
                  [ 0 1 4 0    1 4 5   ]; ... % z=0 supp face
                  [ 0 1 4 0    4 8 5   ]; ... % z=0 supp face
                  [ 1 2 0 0    1 4 2 6 ]; ... % tetrahedron
                  [ 1 2 0 0    6 2 3 4 ]; ... % tetrahedron
                  [ 1 2 0 0    4 3 6 7 ]; ... % tetrahedron
                  [ 1 2 0 0    4 1 5 6 ]; ... % tetrahedron
                  [ 1 2 0 0    4 6 5 8 ]; ... % tetrahedron
                  [ 1 2 0 0    4 7 6 8 ]  ... % tetrahedron
                } ;
       
% ----------------------------------------------------------------------

%% --- Analysis parameters ---
analysisSettings.methodName = 'newtonRaphson' ;
analysisSettings.stopTolIts       = 30      ;
analysisSettings.stopTolDeltau    = 1.0e-12 ;
analysisSettings.stopTolForces    = 1.0e-12 ;
analysisSettings.finalTime        = 1       ;
analysisSettings.deltaT           = .1      ;

%~ controlDofs = [ 7 1 1 ] ;

%~ storeBoolean = 1 ;

%% Output parameters
otherParams.plotParamsVector = [ 3 ] ;
%~ printflag = 2 ;

%~ % --- Analytic sol ---
%~ analyticSolFlag        = 2 ;
%~ analyticCheckTolerance = 1e-8 ;
%~ analyticFunc           = @(w) 1/p * E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) ) ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


return
Conec = {[ 0 1 1 0 0   5 8 6   ]; ... % loaded face
         [ 0 1 1 0 0   6 8 7   ]; ... % loaded face
         [ 0 1 0 0 1   4 1 2   ]; ... % x=0 supp face
         [ 0 1 0 0 1   4 2 3   ]; ... % x=0 supp face
         [ 0 1 0 0 2   6 2 1   ]; ... % y=0 supp face
         [ 0 1 0 0 2   6 1 5   ]; ... % y=0 supp face
         [ 0 1 0 0 3   1 4 5   ]; ... % z=0 supp face
         [ 0 1 0 0 3   4 8 5   ]; ... % z=0 supp face
         [ 1 2 0 0 0   1 4 2 6 ]; ... % tetrahedron
         [ 1 2 0 0 0   6 2 3 4 ]; ... % tetrahedron
         [ 1 2 0 0 0   4 3 6 7 ]; ... % tetrahedron
         [ 1 2 0 0 0   4 1 5 6 ]; ... % tetrahedron
         [ 1 2 0 0 0   4 6 5 8 ]; ... % tetrahedron
         [ 1 2 0 0 0   4 7 6 8 ]  ... % tetrahedron
        } ;

iniMatUs = matUs ;
storeBoolean = 0 ;

ONSAS

% --------------------------------------------------------

clear iniMatUs


controlDispsValsCase1         = controlDisps  ;
loadFactorAnalyticalValsCase1 = analyticVals  ;
loadFactorNumericalValsCase1  = numericalVals ;

close all

% --------------------------------------------------------
% solid model using gmsh mesh, local tension load and complex step 
% --------------------------------------------------------

problemName = 'uniaxialExtension_GMSH_ComplexStep' ;

[ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;

loadsParams{1,1}    = [ 0 1  0 0 0 0 p 0 ] ; % local coords appliend tension

elementsParams{2,1} = [ 4 1 ] ; % complex step constitutive tensor

plotParamsVector = [ 0 ] ;
analyticSolFlag        = 0 ;

% run ONSAS
ONSAS

controlDispsValsCase2         = controlDisps  ;
loadFactorNumericalValsCase2  = numericalVals ;


% --------------------------------------------------------
% truss element model
% --------------------------------------------------------

problemName = 'uniaxialExtension_truss' ;

Nodes = [ 0    0    0 ; ...
          Lx   0    0   ...
        ] ;

Conec = {[ 0 1 0 0 1   1   ] ; ... % fixed node
         [ 0 1 1 0 2   2   ] ; ... % loaded node
         [ 1 2 0 1 0   1 2 ]   ... % truss element
        } ;

% ======================================================================
% --- MELCS parameters ---

materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

% --- Material parameters ---
E = 1 ; nu = 0.3 ;
materialsParams{1} = [ 0 2 E nu ] ;

% --- Element parameters ---
elementsParams = { 1  ; [ 2 0 ]} ;

% --- Load parameters ---
loadsParams{1,1} = [ 1 1  p 0 0 0 0 0 ] ;

% --- CrossSection parameters ---
crossSecsParams = { [ 2 Ly Lz] } ; %

% ----------------------------------------------------------------------
% --- springsAndSupports parameters ---
springsParams{1, 1} = [ inf 0  inf 0   inf 0 ] ;
springsParams{2, 1} = [ 0   0  inf 0   inf 0 ] ;

% ======================================================================

plotParamsVector       = [ 0 ] ;


controlDofs = [ 2 1 1 ] ;

%% run ONSAS
ONSAS

controlDispsValsCase3         = controlDisps  ;
%~ loadFactorNumericalValsCase3  = numericalVals .* (1+controlDisps) / Lx ;
loadFactorNumericalValsCase3  = numericalVals ;


% ----------------------------------------------------------------------
% --- plots ---
lw = 2.0 ; ms = 10 ; plotfontsize = 22 ;

figure, grid on, hold on

plot( controlDispsValsCase1, ...
      loadFactorAnalyticalValsCase1 ,'b-o' , 'linewidth', lw,'markersize',ms )

plot( controlDispsValsCase1, ...
      loadFactorNumericalValsCase1  ,'k-s' , 'linewidth', lw,'markersize',ms)

plot( controlDispsValsCase2, ...
      loadFactorNumericalValsCase2  ,'r-x' , 'linewidth', lw,'markersize',ms)

plot( controlDispsValsCase3, ...
      loadFactorNumericalValsCase3  ,'g--' , 'linewidth', lw,'markersize',ms)

labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic Sol','numerical Sol 1','numerical Sol 2','numerical Sol 3','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
%~ print( [ 'plotsExtensionSVK' ] ,'-dpdflatex','-tight') ;
print( [ 'plotsExtensionSVK.png' ] ,'-dpng') ;

% ----------------------------------------------------------------------
