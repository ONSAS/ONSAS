% =======================================
% sript for Heat transfer code validation
% =======================================


% =========   case 1 =========
% 1D problem
% boundary conditions: dirichlet on both ends
% initial temperature profile
% https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions

timeIncr  = 0.0001 ;
Tfinal    = 0.005   ;
rho       = 1.     ;
cSpHe     = 1.     ;
kCond     = 4      ;
Lx        = 1      ;
Ly        = .25     ;
Lz        = .25     ;
Tdiri     = 0      ;

diriDofs = [];
robiDofs = [] ;

nPlots = inf ;

problemName = 'dirichlet3D' ;

initialTempFunc = 'myInitialTemp' ;

ndivs     = [ 200 1 1 ];

hConv = [];
Tamb = [];
qInpLeft = [];
qInpRight = [];
Tdiri = 0 ;

qInp = 0;

diriFaces = [ Tdiri 1 2 ];
neumFaces = [ qInp  3 4 5 6 ] ;
robiFaces = [  ] ;
  
[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, ...
  diriFaces, neumFaces, robiFaces  );


% evolution plot
figure, grid on, hold on
MS = 20 ; LW = 1.5 ;
indplot = round( ndivs(1)/2 ) ;

% numerical solution plot
plot( times, Ts(indplot,:), 'b',  'markersize', MS,'linewidth',LW )

% analytic solution computation
xsAnly = NodesCoord(indplot,1) ;
alpha = kCond / ( rho * cSpHe ) ;
tanali = exp(-(  pi)^2 * alpha * times ) * sin(pi * xsAnly) ...
                  + exp(-(3*pi)^2 * alpha * times ) * 0.5 * sin( 3 * pi * xsAnly ) ;

% analytic solution plot
plot( times, tanali, 'r',  'markersize', MS,'linewidth',LW )
xlabel('t'), ylabel('Temp')
% ----------------------------------
