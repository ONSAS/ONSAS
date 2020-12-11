%% Example uniaxialSolid
% Linear elastic solid submitted to uniaxial loading.
% Geometry given by $L_x$, $L_y$ and $L_z$, tension $p$ applied on
% face $x=L_x$.
% We use Neo-Hookean constitutive law. The solid is in a compression regime.
%%

% uncomment to delete variables and close windows
clear all, close all

%% General data
dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName = 'uniaxialCompression_Manual' ;

%% Structural properties

% compression applied and x, y, z dimensions
p = -3 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;

% an 8-node mesh is considered with its connectivity matrix
Nodes = [ 0    0    0 ; ...
          0    0   Lz ; ...
          0   Ly   Lz ; ...
          0   Ly    0 ; ...
          Lx   0    0 ; ...
          Lx   0   Lz ; ...
          Lx  Ly   Lz ; ...
          Lx  Ly    0 ] ;

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


% ======================================================================
% --- MELCS parameters ---

materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

% --- Material parameters ---
E = 1 ; nu = 0.3 ;
materialsParams{1} = [ 0 3 E nu ] ;

% --- Element parameters ---
elementsParams{1,1} = [ 5   ] ;
elementsParams{2,1} = [ 4 1 ] ; # the second index indicates that we take the complex-step computation expression

% --- Load parameters ---
loadsParams{1,1} = [ 1 1  p 0 0 0 0 0 ] ;

% --- CrossSection parameters ---

% ----------------------------------------------------------------------
% --- springsAndSupports parameters ---
springsParams{1, 1} = [ inf 0  0   0   0   0 ] ;
springsParams{2, 1} = [ 0   0  inf 0   0   0 ] ;
springsParams{3, 1} = [ 0   0  0   0   inf 0 ] ;

% ======================================================================


%% --- Analysis parameters ---
stopTolIts       = 30      ;
stopTolDeltau    = 1.0e-12 ;
stopTolForces    = 1.0e-12 ;
targetLoadFactr  = 1       ;
nLoadSteps       = 10      ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ;

controlDofs = [ 7 1 1 ] ;

%% Output parameters
plotParamsVector = [ 3 ] ;
printflag = 2 ;

reportBoolean = 0;

% --- Analytic sol ---
analyticSolFlag        = 0 ; % not available

%% run ONSAS
acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;
% --------------------------------------------------------
