%% Three dimensional truss tower
% In this example dxf reading feature is shown. A simple 3D truss structure is 
% considered for a nonlinear elastic analysis.
%
%%

clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ; addpath( dirOnsas );

% Problem name in order to write LaTex report 
problemName = 'dxf3DTrussTower' ; 

% scalar auxiliar parameters
E = 210e9 ;  A = 2.5e-4 ; nu = 0 ;  rho = 0 ; 

% ----------------------------------------------------------------------
% MELCS parameters
% Materials
% structure of cell: [ density model YoungModulus Poisson coefficient ]
materialsParams = {[ rho 1 E nu ]} ;
% Elements
% Structure of the cell: { elem_1 ; ... ; elem_n }
elementsParams  = { 1; 2} ;
% Loads
% Structure of the cell: { [ loadLabel global/local Fx Mx Fy My Fz Mz ] }
loadsParams     = { [ 1 1   0 0 10e4 0 -10e4 0] ;
										[ 2 1   0 0 15e4 0 -5e4 0] } ;
% Cross-Sections
crossSecsParams = { [ 2 .1 .1 ] } ;
% Springs
% structure of the cell: { [ ux thetax uy thetay uz thetaz ] }
springsParams   = { [ inf  0  inf  0  inf   0 ]} ;

% Nodes and Conectivity matrix from .dxf file
addpath( [ dirOnsas '/sources' ] );
[ Nodes, Conec ] = meshFileReader( 'torre.dxf' ) ;

% Output parameters
plotParamsVector = [ 3 ] ;

% run ONSAS
ONSAS
