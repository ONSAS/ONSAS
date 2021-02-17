% ------------------------------------------------------------------------------
% example simple transmission line problem
% ------------------------------------------------------------------------------

clear all, close all

dirOnsas = [ pwd '/../../src' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName = 'simpleTransmissionLine' ;

% Geometric properties cable:
%-------------------- 
dc  = 0.1 ; lc  = 406.5;

%Material properties cable :
%--------------------
Ec   = 70e9    ;
nuc  = 0.3      ; 
rhoc = 2700	    ;
Gc	= Ec/(2*(1+nuc));	

%Elem CrosSec and elem Params Cells:
%--------------------
%1 cable 2 isolator chain
materialsParams = {[rhoc 1 Ec nuc ] } ;

elementsParams  = { 1; 3} ;

%Conec Elem 2 type secs	
crossSecsParams = {[3 dc] } ;	

loadsParams    = { [ 1  1    0 0 0 0 -1 0 ] } ;

springsParams  = { [ inf inf  inf inf  inf inf ] } ;

%Nodes and Conec matrix:
%-----------------------

%Nelem Cable
NelemC  = 10 ;
NnodesC = NelemC+1;

%             x                       y                       z
Nodes = [ (0:(NelemC))'*lc/NelemC    zeros(NnodesC,1)         zeros(NnodesC,1) ] ;

Conec = {} ;

             % M E L C S
Conec{1,1} = [ 0 1 0 0 1  1 ] ;
Conec{2,1} = [ 0 1 0 0 1  NnodesC ] ;
Conec{3,1} = [ 0 1 1 0 0  round(NnodesC/2) ] ;

for i=1:NelemC
  Conec{i+3,1} = [ 1 2 0 1 0  i i+1 ] ; 
end


%Parameters and tolerances:
%-----------------------
% times
timeIncr   =  1e-2       ;
finalTime  = timeIncr*10 ;

nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 1e-4        ; 
stopTolForces = 1e-4        ;
stopTolIts    = 30          ;
alphaHHT      = -0.05       ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;
    
%Loads and supports:
%Loads---------------------

loadFactorsFunc = @(t)(  t );

% Damping
nodalDispDamping = .1 ;
%Booleans control and plot:
%-----------------------
%controlDof
controlDofs = [ NelemC/2 5 1 ] ;
%Plots
plotParamsVector = [ 0 ];

% RunONSAS
ONSAS;
