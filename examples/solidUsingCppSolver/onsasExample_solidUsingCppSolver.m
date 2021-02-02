% --------------------------------------------------
% solid analysis example for testing C++ solver
% --------------------------------------------------

clear all, close all

%% General data
dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path
problemName = 'solidUsingCppSolver_case1' ;

%% Structural properties

% tension applied and x, y, z dimensions
p = 3 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;


% ======================================================================
% --- MELCS parameters ---

% tension applied and x, y, z dimensions
p = 3 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;

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

% --- Material parameters ---
E = 1 ; nu = 0.3 ;
materialsParams = {[ 0 2 E nu ]} ;

% --- Element parameters ---
elementsParams = { [ 5   ] ; ...
                   [ 4 2 ] } ; % analytic constitutive tensor

% --- Load parameters ---
loadsParams = {[ 1 1  p 0 0 0 0 0 ]} ;  % global coords tension applied

% --- CrossSection parameters ---
crossSecsParams = cell(1,1) ; %

% --- springsAndSupports parameters ---
springsParams = {[ inf 0  0   0   0   0 ] ; ...
                 [ 0   0  inf 0   0   0 ] ; ...
                 [ 0   0  0   0   inf 0 ] } ;

% ----------------------------------------------------------------------

%% --- Analysis parameters ---
stopTolIts       = 30      ;
stopTolDeltau    = 1.0e-12 ;
stopTolForces    = 1.0e-12 ;
targetLoadFactr  = 1       ;
nLoadSteps       = 2      ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ;

controlDofs = [ 7 1 1 ] ;

%% Output parameters
plotParamsVector = [ 3 ] ;
printflag = 2 ;

% --- Analytic sol ---
analyticSolFlag        = 0 ;
analyticCheckTolerance = 1e-8 ;
analyticFunc           = @(w) 1/p * E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) ) ;

%% run ONSAS
ONSAS


controlDispsValsCase1         = controlDisps  ;
%~ loadFactorAnalyticalValsCase1 = analyticVals  ;
%~ loadFactorNumericalValsCase1  = numericalVals ;

% ==============================================================================

problemName = 'solidUsingCppSolver_case2' ;

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

cppSolverBoolean = 1 ;
%% run ONSAS
ONSAS

% ==============================================================================


problemName = 'solidUsingCppSolver_case3' ;

[ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;

cppSolverBoolean = 0 ;

%% run ONSAS
ONSAS
% ==============================================================================


problemName = 'solidUsingCppSolver_case4' ;

[ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;
cppSolverBoolean = 1 ;

%% run ONSAS
ONSAS

% ==============================================================================


problemName = 'solidUsingCppSolver_case5' ;

[ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension_dense.msh' ) ;

%% run ONSAS
ONSAS


%~ problemName = 'solidUsingCppSolver_case6' ;

%~ [ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension_dense.msh' ) ;
%~ cppSolverBoolean = 0 ;

%~ %% run ONSAS
%~ ONSAS
