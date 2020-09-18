% ---------------------------
% 1D Heat Transfer Problem
% ---------------------------

function onsasExample_1DHeatTransfer( onsasDir, scalarParams, plotBoolean )

close all, problemName = '1DHeatTransfer' ;

% --- scalar parameters ---

% heat transfer system 
if nargin < 2

  timeIncr   = 0.001 ;
  Tfinal = 1 ;
  rho    = 1 ;
  cSpHe  = 1 ;
  kCond  = .5 ;
  L      = 2.75 ; % domain [0,L]
  Area   = 1 ;
  hConv  = 10 ;
  nelem  = 10  ;

  anlyBoolean = 0 ;              wx = 1 ;

  Tamb        = 20       ;
  diriDofs    = [ ] ;  robiDofs    = [ 1 nelem+1 ] ;  

else
  %~ k        = scalarParams(1) ;
  %~ c        = scalarParams(2) ;
  %~ m        = scalarParams(3) ;
  %~ omegaBar = scalarParams(4) ;
  %~ p0       = scalarParams(5) ;
  %~ u0       = scalarParams(6) ;
  %~ timeIncr = scalarParams(7) ;  
end

if nargin < 3
  plotBoolean = 1 ;
  close all
end

dt = timeIncr ;

if exist('nt')==0, nt = Tfinal / dt ; end

if exist('nCurves') == 0, nCurves = 15 ; end

alpha = kCond / ( rho * cSpHe ) ;

nnodes = nelem +1 ;

nnodesAnly = 100 ;

xs     = linspace(0, L, nnodes     )' ;
xsAnly = linspace(0, L, nnodesAnly )' ;


neumdofs = 1:nnodes ;
neumdofs( diriDofs ) = [] ;

lelem  = L/nelem ;

% local elemental diffussion equation
Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;

MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;

bQhe = 0.5 * Area * lelem * [ 1 ; 1 ] * 0 ;

% initial temperature
%~ T0     = sin(pi*xs*wx) + 0.5*sin(3*pi*xs*wx) ;
T0     = 0.5*Tamb * ones( size( xs ) ) ;
T0     = 0.5*Tamb * ( 1 + .1* sin( pi * xs / L ) ) ;
% ---

Ts     = T0 ;

   
% ------------------------
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
end

KdiffG = KdiffG + MrobiG ;

% ------------------------

CDD = MintEG(diriDofs, diriDofs) ;
CND = MintEG(neumdofs, diriDofs) ;
CNN = MintEG(neumdofs, neumdofs) ;

qext = zeros( nnodes, 1 ) ;
%~ if exist( 'qentrIzq' ) ~= 0, qext(  1) = qentrIzq ; end
%~ if exist( 'qentrDer' ) ~= 0, qext(end) = qentrDer ; end

if ~isempty( robiDofs )
  qext( robiDofs ) = qext( robiDofs ) + hConv * Tamb ;
end

Tamb
T0
hConv

%~ return


qext = qext + QhG 


if plotBoolean
  figure
  hold on, grid on
end

MS = 10 ; 
LW = 1.5 ; 

%~ nt

fext = zeros( size(Ts) ) ;

% ------------------------
% time loop
for i=0:nt
  i
  t = i*dt 

  if i~=0
    fext ( :,i) = qext ( neumdofs    ) ;
  
    f = (    qext ( neumdofs    ) * dt ...
          - KdiffG( neumdofs, : ) * Ts( :, i ) * dt ...
          + MintEG( neumdofs, : ) * Ts( :, i ) 
        ) 
    
    Tip1 = CNN \ f 
    Ts( neumdofs, i+1 ) = Tip1 ;
    
    if ~isempty( diriDofs ),
      Ts( diridofs, i+1 ) = Tdiri   ;
    end
  end  

  %~ if anlyBoolean
    %~ if caseNum == 1 || caseNum == 2
      %~ TsAnly (:,i+1) = exp(-(  pi* wx*alpha )^2 * t ) *       sin(     pi * xsAnly * wx ) ...
                     %~ + exp(-(3*pi* wx*alpha )^2 * t ) * 0.5 * sin( 3 * pi * xsAnly * wx ) ;
    %~ end
  %~ end
  
  % --- plots ---
  if plotBoolean
    if mod(i, round(nt/nCurves) )==0
      plot( xs    , Ts    (:, i+1), 'b-o', 'markersize', MS,'linewidth',LW );  
      
      if anlyBoolean
        plot( xsAnly, TsAnly(:, i+1), 'r--'  , 'markersize', MS,'linewidth',LW );  
      end
    end
  end
    % ---------------
  
end
% ------------------------

if plotBoolean
  axis equal
  print( '../../1DheatTransfer.png','-dpng')

  %~ figure, plot( Ts(2,:) )
end

KdiffGNN = KdiffG(neumdofs, neumdofs ) ;

%~ save -mat auxVars.mat KdiffGNN CNN Ts

M = zeros( nnodes, nnodes ) ;
K = KdiffG ;
C = MintEG ;

us = [] ;
udots = [] ;

save -mat outputMatrices.mat  K C M fext timeIncr us udots
