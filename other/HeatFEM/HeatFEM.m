% ---------------------------
% Heat Transfer FEM solver
% ---------------------------

function Ts = HeatFEM( ...
  timeIncr, Tfinal, ...
  materialParams, ...
  geometryParams, ...
  meshParams, ...
  boundaryCondParams,
  plotBoolean, nCurves, anlyBoolean, problemName )

% ==============================================
% seteos previos

% variables auxiliares para futuro en ONSAS
close all

% material params
rho   = materialParams(1) ;  
cSpHe = materialParams(2) ;
kCond = materialParams(3) ;

%
geometryType = geometryParams(1) ;

if geometryType == 1 % unidimensional
  L = geometryParams(2); A = geometryParams(3);
elseif geometryType == 2 % solid 
  Lx = geometryParams(2);
  Ly = geometryParams(3);
  Lz = geometryParams(4);
end

% mesh creation
if geometryType == 2
  ndivs = meshParams (1:3) ;
end

[NodesCoord, Conec] = createSolidRegularMesh( ndivs, Lx, Ly, Lz )


plotTemp( NodesCoord, Conec,0,zeros(size(NodesCoord,1),1) )

%~ hConv,Tamb, Tdiri diriDofs, robiDofs, qInpLeft, qInpRight



% defino/renombro variables de entrada
dt = timeIncr ;

if exist('nt')==0, nTimes = Tfinal / dt + 1 ; end

tiempos = linspace( 0, Tfinal, nTimes )' ;

if isempty( nCurves ),
  nCurves = 4 ;
elseif nCurves == inf,
  nCurves = nTimes ;
end

alpha = kCond / ( rho * cSpHe ) ;

return


nnodes = nelem +1 ;
nnodesAnly = 100 ;

% discretizacion puntos para ploteos
xs     = linspace(0, L, nnodes     )' ;
xsAnly = linspace(0, L, nnodesAnly )' ;

% grados de libertad de neumann
neumdofs = 1:nnodes ;
neumdofs( diriDofs ) = [] ;

% calculo largo de cada elemento
lelem  = L/nelem ;
% ==============================================



% ==============================================
% construccion matrices FEM


% local elemental diffussion equation
Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;

MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;

% volumetric external heat source
bQhe = 0.5 * Area * lelem * [ 1 ; 1 ] * 0 ;

% initial temperature
if initialTempFlag == 1
  T0 = sin(pi*xs) + 0.5*sin(3*pi*xs) ;
end
%~ T0 = 15.0 * ( 1 + .1* sin( pi * xs / L ) ) ;
%T0 = 0.5*Tamb * ( 1 + .1* sin( pi * xs / L ) ) ;


% ------------------------
% Ts matriz para almacenar temperatuas en cada tiempo
Ts = T0 ;

% matrices/vectors assembly
KdiffG = zeros( nnodes, nnodes ) ;
MintEG = zeros( nnodes, nnodes ) ;
MrobiG = zeros( nnodes, nnodes ) ;
QhG    = zeros( nnodes, 1      ) ;

for i = 1 : nelem
  nodeselem = [ i i+1 ] ;

  elemDofs = nodes2dofs ( nodeselem, 1 ) ;

  KdiffG( elemDofs , elemDofs ) = ...
  KdiffG( elemDofs , elemDofs ) + Kdiffe  ;

  MintEG( elemDofs , elemDofs ) = ...
  MintEG( elemDofs , elemDofs ) + MintEe  ;

  QhG   ( elemDofs            ) = ...
  QhG   ( elemDofs            ) + bQhe  ;

end

if  ~isempty( robiDofs )
  MrobiG ( robiDofs, robiDofs ) = hConv ;
  KdiffG = KdiffG + MrobiG ;
end


% ------------------------

CDD = MintEG( diriDofs, diriDofs ) ;
CND = MintEG( neumdofs, diriDofs ) ;
CNN = MintEG( neumdofs, neumdofs ) ;

KdiffGNN = KdiffG(neumdofs, neumdofs ) ;
KdiffGND = KdiffG(neumdofs, diriDofs ) ;


Matrix = ( KdiffGNN * dt + CNN ) ;

qext = zeros( nnodes, 1 ) ;

% input fluxes
if ~isempty( qInpLeft  ), qext(   1) = qInpLeft ; end
if ~isempty( qInpRight ), qext( end) = qInpRight; end

if ~isempty( diriDofs )
  qext( neumdofs ) = qext( neumdofs ) - KdiffGND * ones( length( diriDofs ), 1 ) * Tdiri ;
end

if ~isempty( robiDofs )
  qext( robiDofs ) = qext( robiDofs ) + hConv * Tamb ;
end

qext = qext + QhG ;
% ==============================================


% =========================================
% simulacion

if plotBoolean
  figure, hold on, grid on
end

% plot settings
MS = 10 ; LW = 1.5 ;

fext = zeros( size(Ts) ) ;

% ------------------------
% time loop

indsPlotStep = round( nTimes / nCurves ) 

for ind = 1:nTimes %ind es el indice de tiempo que se esta hallando
  t = dt*(ind-1) ;

  fext ( neumdofs, ind ) = qext ( neumdofs    ) ;
  
  if ind > 1
    f = (    fext( neumdofs, ind ) * dt  ...
         + MintEG( neumdofs, : ) * Ts( :, ind-1 )
        );

    Tip1 = Matrix \ f ;
    Ts( neumdofs, ind ) = Tip1 ;

    if ~isempty( diriDofs ),
      Ts( diriDofs, ind ) = Tdiri   ;
    end

  end % if i not zero

  if anlyBoolean
    TsAnly(:,ind) = exp(-(  pi)^2 * alpha * t ) * sin(pi * xsAnly) ...
                  + exp(-(3*pi)^2 * alpha * t ) * 0.5 * sin( 3 * pi * xsAnly ) ;
  end

  % --- plots ---
  if plotBoolean
    if ind == 1 || mod( ind, indsPlotStep ) == 0
      plot( xs, Ts(:, ind), 'b-o', 'markersize', MS,'linewidth',LW );
      if anlyBoolean
        plot( xsAnly, TsAnly(:, ind), 'r--'  , 'markersize', MS,'linewidth',LW );
      end
    end
  end
  % ---------------


end % for times
% ------------------------

if plotBoolean
  %~ axis equal
  if anlyBoolean
    legend('numerical','analytical')
  end
  print( ['pngs/' problemName '.png'],'-dpng')
end

% ploteos en puntos específicos
if plotBoolean
  figure, hold on, grid on
  plot(tiempos, Ts(round(nnodes/2),:), 'b-o','linewidth',LW,'markersize',MS)
  %~ plot(tiempos, Ts(end,:), 'r-x','linewidth',LW,'markersize',MS)
  %~ legend('izq', 'der')
  xlabel('t'), ylabel('Temp')
end

M = zeros( nnodes, nnodes ) ;
K = KdiffGNN ;
C = CNN ;

us    = Ts ;
udots = [] ;

save -mat 'mats/outputMatrices.mat'  K C M fext timeIncr us udots xs






