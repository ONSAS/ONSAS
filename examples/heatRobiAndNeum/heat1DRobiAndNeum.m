% =========   caso de validacion 2 =========
%  diri no homogeneo ( ver fundamentos en https://www.math.uzh.ch/li/index.php?file&key1=25297 )

Tdiri     = 1  ;
%~ Tfinal = .2 ;
Tfinal = .02 ;
ndivs = [ 50 ] ;
diriDofs = [ 1 ndivs+1 ]
problemName = 'nonHomogDirichlet' ;

[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 L Area ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], [], [], [] );


% numerical solution plot
figure
plot( times, Ts(indplot,:), 'b',  'markersize', MS,'linewidth',LW )



% =========   caso de validacion 3 =========
%  conv and neum

diriDofs = [ ];
robiDofs = [1 ] ;
qInpRight = 4 ;

 Tfinal = 4.5;
Tfinal = 0.5;

Tamb      = 5 ;

hConv = 2;

problemName = 'robinAndNeumann' ;

[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 L Area ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, [], [], [], [] );






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
