% =======================================
% sript for Heat transfer code validation
% =======================================


% =========   case 1 =========
% 3D problem
% boundary conditions: dirichlet on both ends
% initial temperature profile
% https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions

timeIncr  = 0.0001 ;
Tfinal    = 0.01   ;
rho       = 1.     ;
cSpHe     = 1.     ;
kCond     = 4      ;
Lx        = 1      ;
Ly        = .5     ;
Lz        = .5     ;
Tdiri     = 0      ;

diriDofs = [];
robiDofs = [] ;

nPlots = inf ;

problemName = 'dirichlet3D' ;

initialTempFunc = 'myInitialTemp' ;

ndivs     = [ 20 2 2 ];

hConv = [];
Tamb = [];
qInpLeft = [];
qInpRight = [];
Tdiri = 0 ;

qInp = 0;

diriFaces = [ Tdiri 1 2 ];
neumFaces = [ qInp  3 4 5 6 ] ;
robiFaces = [  ] ;
  
[Ts3D, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], ...
  diriFaces, neumFaces, robiFaces  );


% =========   case 2 =========
% 1D problem
% boundary conditions: dirichlet on both ends
% initial temperature profile
% https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/

problemName = 'dirichlet1D' ;
diriDofs = [ 1 ndivs(1)+1 ] ;
nPlots = 0 ;

[ Ts1D, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 Lx Ly*Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], ...
  diriFaces, neumFaces, robiFaces  );


% evolution plot
figure, grid on, hold on
MS = 15 ; LW = 2 ;
indplot = round( ndivs(1)/2 ) ;

% numerical solution plot
plot( times(1:5:end), Ts3D(indplot,1:5:end), 'b-s',  'markersize', MS,'linewidth',LW )

% analytic solution computation
xsAnly = NodesCoord(indplot,1) ;
alpha = kCond / ( rho * cSpHe ) ;
tanali = exp(-(  pi)^2 * alpha * times ) * sin(pi * xsAnly) ...
                  + exp(-(3*pi)^2 * alpha * times ) * 0.5 * sin( 3 * pi * xsAnly ) ;

% analytic solution plot
plot( times(1:5:end), tanali(1:5:end), 'r-o',  'markersize', MS,'linewidth',LW )
xlabel('t'), ylabel('Temp')


plot( times(1:5:end), Ts1D(indplot,1:5:end), 'g-x',  'markersize', MS,'linewidth',LW )

legend('num3D','analytic','num1D','location','southeast')

print('pngs/valid3D.png','-dpng')
% ----------------------------------
