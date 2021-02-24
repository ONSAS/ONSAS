% =======================================
% sript for Heat transfer code validation
% =======================================


% =========   case 1 =========
% 1D problem
% boundary conditions: dirichlet on both ends
% initial temperature profile
% https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/
% paramters table 2.8 from Bofang's book

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions

global rho       = 2485  ;

timeIncr  = 1/3   ;
Tfinal    = 24*10    ;
cSpHe     = 0.967 ;
kCond     = 9.37  ;
Lx        = 1     ;
Ly        = 1     ;
Lz        = 1      ;
Tdiri     = 0      ;

diriDofs = [];
robiDofs = [] ;

nPlots = inf ;

problemName = 'dirichlet3D' ;

initialTempFunc = 'uniformTemp' ;
internalHeatFunc = 'hydrationHeat' ;
%~ initialTempFunc = 'myInitialTemp' ;
%~ internalHeatFunc = 'myInternalHeat' ;

ndivs     = [ 12 12 12 ];

hConv = [];
hConvAir =  80 ;
hConvSol = 800 ;

Tamb     = 23 ;
qInpLeft = [];
qInpRight = [];
Tdiri = 0 ;

qInp = 0 ;

diriFacesAndVals = [  ] ;
neumFacesAndVals = [ 1 0 ; 3 0 ] ;
robiFacesAndVals = [ 2 hConvSol ;  4 hConvSol ; 5 hConvSol ; 6 hConvAir ] ;

[Ts3D, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, internalHeatFunc, ...
  diriFacesAndVals, neumFacesAndVals, robiFacesAndVals  );

% evolution plot
figure, grid on, hold on
MS = 15 ; LW = 2 ;
%~ indplot = round( ndivs(1)/2 ) ;
indplot = round( prod((ndivs+1))/2 ) ;

% numerical solution plot
plot( times(1:5:end), Ts3D(indplot,1:5:end), 'b-s',  'markersize', MS,'linewidth',LW )
xlabel('time (s)')
ylabel( ['Temp at ' num2str(indplot) ' (C)' ] )
