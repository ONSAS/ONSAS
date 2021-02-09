% ---------------------------
% Heat Transfer FEM prototype solver
% ---------------------------

function [Ts, NodesCoord, tiempos] = HeatFEM( ...
  timeIncr, Tfinal, ...
  materialParams, ...
  geometryParams, ...
  meshParams, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc )

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
  L = geometryParams(2); Area = geometryParams(3);
elseif geometryType == 2 % solid 
  Lx = geometryParams(2);
  Ly = geometryParams(3);
  Lz = geometryParams(4);
end

% mesh creation
if geometryType == 1
  ndivs = meshParams (1) ;
elseif geometryType == 2
  ndivs = meshParams (1:3) ;
end

if geometryType == 1
  NodesCoord = linspace(0, L, ndivs(1)+1 )' ;
  Conec = [ (1:(ndivs(1)))' (2:(ndivs(1)+1))' ] ;
elseif geometryType == 2
  [NodesCoord, Conec] = createSolidRegularMesh( ndivs, Lx, Ly, Lz )
end



% defino/renombro variables de entrada
dt = timeIncr ;

if exist('nt')==0, nTimes = Tfinal / dt + 1 ; end

tiempos = linspace( 0, Tfinal, nTimes )' ;

if isempty( nPlots ),
  nPlots = 4 ;
  plotBoolean = 1,
elseif nPlots == 0,
  plotBoolean = 0,
else
  plotBoolean = 1,
  if nPlots == inf,
    nPlots = nTimes ;
  end
end


nnodes = size(NodesCoord,1);
nelem  = size(Conec,1) ;


% grados de libertad de neumann
neumdofs = 1:nnodes ;
neumdofs( diriDofs ) = [] ;

% ==============================================



% ==============================================
% construccion matrices FEM

if geometryType == 1
  % calculo largo de cada elemento
  lelem  = L/nelem ;
  
  % local elemental diffussion equation
  Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;
  
  MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;
  
  % volumetric external heat source
  bQhe = 0.5 * Area * lelem * [ 1 ; 1 ] * 0 ;

  % initial temperature
  T0 = feval( initialTempFunc, NodesCoord ) ;
end

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

fext = zeros( size(Ts) ) ;

% ------------------------
% time loop

indsPlotStep = round( nTimes / nPlots ) 

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

  % --- plots ---
  if plotBoolean
    if ind == 1 || mod( ind, indsPlotStep ) == 0
      plotTemp( NodesCoord, Conec,0, Ts(:,ind), geometryType )
    end
  end
  % ---------------

end % for times
% ------------------------

if plotBoolean
  print( [ problemName '.png'],'-dpng')
end

M = zeros( nnodes, nnodes ) ;
K = KdiffGNN ;
C = CNN ;
xs = NodesCoord(:,1) ;

us    = Ts ;
udots = [] ;

save -mat 'outputMatrices.mat'  K C M fext timeIncr us udots xs

