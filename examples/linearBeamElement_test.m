% ======================================================================
% Linear Beam Element example

% ===========================
% first case: linear analysis

close all, clear all        ;

dirOnsas = [ pwd '/../src' ] ; % set ONSAS.m directory
addpath( dirOnsas )             ; % add ONSAS directory to path

problemName = 'linearBeamElement_test' ; %

% ----------------------------------------------------------------------
% scalar auxiliar parameters
E = 210e6 ; L = 5 ; nu = 0.3 ;  rho = 0 ; b = 0.3 ;

% ----------------------------------------------------------------------
% MELCS parameters

% Materials
materialsParams = {[ rho 1 E nu ]} ;

% Elements
elementsParams  = { 1; 5 } ;

% Loads
loadsParams   = { [ 1 1   0 0 0  0 -5 0] } ;

% Cross-Sections
crossSecsParams = { [ 2 b b ] } ;

% springs parameters
springsParams = { [ inf 0 inf 0 inf 0  ] ; ...
                  [ 0   0 inf 0 inf 0  ] } ;
% ----------------------------------------------------------------------
% mesh parameters

% node coordinates matrix
Nodes = [ 0		0	0		; ...
					L 	0 0 	; ...
          2*L 0 0 	] ;

% connectivity cell
%						M	E	L	C	S
Conec = { [ 0 1 0 0 1  1   ] ; ... % fixed node
          [ 0 1 1 0 0  2   ] ; ... % loaded node
          [ 0 1 0 0 2  3   ] ; ... % fixed node
          [ 1 2 0 1 0  1 2 ] ; ... % beam element
          [ 1 2 0 1 0  2 3 ] } ;   % beam element

plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;

ONSAS
% ----------------------------------------------------------------------

