% ======================================================================
% 

clear all, close all
dirOnsas    = [ pwd '/..' ] ; addpath( dirOnsas );
problemName = '1DHeatTransfer' ;

nelem = 10 ;

%~ aux_toyHeatTransferExample( 1, nelem, 1 ) ;

rho = 1 ; c = 1; k = 1;  E = 1; nu = 0 ; A = 1 ;

% Materials
materialsParams = {[ rho c k 1 E nu ]} ;

% Elements
elementsParams  = { 1; [2 1 1] } ;

% Loads
loadsParams     = { [ 1 1   0 0 0 0 0 0    1] } ;

% Cross-Sections
crossSecsParams = { [ 2 sqrt(A) sqrt(A) ] } ;

% Springs
springsParams   = { [ inf  0  inf  0  inf   0   ]  ...
                    [ 0    0  inf  0    0   0 ] } ;

% ----------------------------------------------------------------------
% nodes coordinates matrix and connectivity cell

% node coordinates matrix
Nodes = [   0  0  0  ; ...
            1  0  0  ] ;

% connectivity cell
Conec = { [ 0 1 1   1   ] ; ... %  node
          [ 0 1 1 0 0  2   ] ; ... % node
          [ 1 2 0 1 0  1 2 ] } ;   % truss element

controlDofs      = [ 2 5 -1 ] ; % [ node nodaldof scalefactor ]
plotParamsVector = [ 3 ] ;
reportBoolean    = 1     ;

analyticSolFlag = 2    ;
analyticFunc    = @(w) 2 * E * A * sin(ang1*pi/180)^2 * w / L ;

ONSAS
