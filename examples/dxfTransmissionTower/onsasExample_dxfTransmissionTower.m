%% Three dimensional truss tower
% In this example dxf reading feature is shown. A simple 3D truss structure is 
% considered for a nonlinear elastic analysis.
%
%%

clear all, close all

dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

% Problem name in order to write LaTex report 
problemName = 'dxfTransmissionTower' ; 

% scalar auxiliar parameters
E = 210e9 ;  A = 2.5e-1 ; nu = 0 ;  rho = 0 ; 

% ----------------------------------------------------------------------
% MELCS parameters
% Materials
% structure of cell: [ density model YoungModulus Poisson coefficient ]
materialsParams = {[ rho 3 E nu ]} ;
% Elements
% Structure of the cell: { elem_1 ; ... ; elem_n }
elementsParams  = { 1; 2} ;
% Loads
% Structure of the cell: { [ loadLabel global/local Fx Mx Fy My Fz Mz ] }
loadsParams     = { [ 1 1   0 0 0   0 -15e4 0] ;
										[ 1 1   0 0 1e2 0 -5e4  0 ] } ;
% Cross-Sections
crossSecsParams = { [ 2 .2 .2 ] } ;
% Springs
% structure of the cell: { [ ux thetax uy thetay uz thetaz ] }
springsParams   = { [ inf  0  inf  0  inf   0 ]} ;

% Nodes and Conectivity matrix from .dxf file
addpath( [ pwd '/../sources' ] ) ;
[ Nodes, Conec ] = meshFileReader( 'geometry_transmissionTower.dxf' ) ;

% analysis parameters
stopTolDeltau    = 1.0e-8 ;    stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 1e2    ;    nLoadSteps       = 10     ;
stopTolIts       = 30     ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

stabilityAnalysisBoolean = 1 ;

% Output parameters
plotParamsVector = [ 3 ] ;

% run ONSAS
run( [  pwd  '/../ONSAS.m' ] ) ;
