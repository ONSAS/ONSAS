% ---------------------------
% 1D Heat Transfer Problem
% ---------------------------

function onsasExample_1DHeatTransfer( onsasDir, scalarParams, plotBoolean )

close all, problemName = '1DHeatTransfer' ;

% --- scalar parameters ---

% heat transfer system 
if nargin < 2

  timeIncr  = 0.001 ;
  Tfinal    = .004 ;
  rho       = 1 ;
  cSpHe     = 1 ;
  kCond  = 2  ;
  L      = 3.4 ; % domain [0,L]
  Area   = 1   ;
  hConv  = 2  ;
  nelem  = 10   ;

  anlyBoolean = 0 ;              wx = 1 ;

  Tamb        = 8       ;
  T0val = 10 ;
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
Kke = kCond * Area * lelem * 1/(lelem^2) * [ 1 -1 ; -1 1 ] ;

%~ Ce  = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;
Ce  = rho * cSpHe * Area * lelem / 2 * [ 1 0 ; 0 1 ] ;

bQhe = 0.5 * Area * lelem * [ 1 ; 1 ] * 1 ;

% initial temperature
%~ T0     = sin(pi*xs*wx) + 0.5*sin(3*pi*xs*wx) ;
T0     = T0val * ones( size( xs ) ) ;
%~ T0     = T0val * ( 1 + .0* sin( pi * xs / L ) ) ;
% ---

Ts     = T0 ;

   
% ------------------------
% matrices/vectors assembly
KkG    = zeros( nnodes, nnodes ) ;
CG     = zeros( nnodes, nnodes ) ;
KcG    = zeros( nnodes, nnodes ) ;
QhG    = zeros( nnodes, 1      ) ;

for i = 1 : nelem
  nodeselem = [ i i+1 ] ;
  
  %~ elemDofs = nodes2dofs ( nodeselem, 1 ) ;
  elemDofs = nodeselem ;
  
  KkG( elemDofs , elemDofs ) = ...
  KkG( elemDofs , elemDofs ) + Kke  ; 
  
  CG( elemDofs , elemDofs ) = ...
  CG( elemDofs , elemDofs ) + Ce  ; 

  QhG   ( elemDofs            ) = ...
  QhG   ( elemDofs            ) + bQhe  ; 

end

if  ~isempty( robiDofs )
  
  KcG ( robiDofs, robiDofs ) = hConv * Area * eye(length(robiDofs)) ;
end

KG  = KkG + KcG ;

% ------------------------

if isempty( diriDofs )
  systMatrix = CG + KG * dt ;
else
  error('add case');
end

qext = zeros( nnodes, 1 ) ;
%~ if exist( 'qentrIzq' ) ~= 0, qext(  1) = qentrIzq ; end
%~ if exist( 'qentrDer' ) ~= 0, qext(end) = qentrDer ; end

% convection
%~ KcG
qext = qext + KcG * Tamb * ( ones(nnodes,1) ) ;

%~ Tamb
%~ T0
%~ hConv

%~ qext 

qext = qext + QhG  ;

if plotBoolean
  figure
  hold on, grid on
end

MS = 10 ; 
LW = 1.5 ; 

%~ nt

fext = qext ( neumdofs    ) ;

% ------------------------
% time loop
for i=1:nt % nt increments
  
  Ti = Ts( :, i ) ;  

  ip1 = i+1  ;
  tp1 = ip1 * dt ;

  fext ( :,ip1) = qext ( neumdofs    ) ;
 
  %~ qext ( neumdofs    ) * dt
  
  %~ a = hConv * ( Tamb - Ts(:,i) )
  
  %~ qext ( neumdofs    ) * dt ...
          %~ - KdiffG( neumdofs, : ) * Ts( :, i ) * dt
  
  qconv = qext ( neumdofs  ) * dt ;

  qabs = CG * Ti ;

  systIndTerm = ( qconv ...
                  + qabs ) ;
    
  Tip1 = systMatrix \ systIndTerm  ;

  Ts( neumdofs, ip1 ) = Tip1 ;
  
  analyIndTeConv = hConv * Tamb * dt ;
  analyIndTeAbs  = rho * cSpHe * lelem * .5 * Ti;
  
  analyIndTe = hConv * Tamb * dt + rho * cSpHe * lelem * .5 * Ti ;
  
  analyTp1 = analyIndTe / ( rho * cSpHe * lelem * .5 + hConv * dt ) ;
    

  %~ if anlyBoolean
    %~ if caseNum == 1 || caseNum == 2
      %~ TsAnly (:,i+1) = exp(-(  pi* wx*alpha )^2 * t ) *       sin(     pi * xsAnly * wx ) ...
                     %~ + exp(-(3*pi* wx*alpha )^2 * t ) * 0.5 * sin( 3 * pi * xsAnly * wx ) ;
    %~ end
  %~ end
  
  % --- plots ---
  if plotBoolean
    if i == 1
      plot( xs    , Ts    (:, 1), 'g-o', 'markersize', MS,'linewidth',LW );
    end
    
    if mod(i, round(nt/nCurves) )==0
      plot( xs    , Ts    (:, ip1), 'b-o', 'markersize', MS,'linewidth',LW );  
      
      
      if anlyBoolean
        plot( xsAnly, TsAnly(:, ip1), 'r--'  , 'markersize', MS,'linewidth',LW );  
      end
    end
  end
    % ---------------
  
end
% ------------------------

%~ return

if plotBoolean
  ylabel('Temperature (^\circC)')
  xlabel('height (m)')
  print( '../../1DheatTransfer.png','-dpng')
  %~ figure, plot( Ts(2,:) )
end


M = zeros( nnodes, nnodes ) ;
K = KG ;
C = CG ;

us = [] ;
udots = [] ;

fext
K
C
M
save -mat outputMatrices.mat  K C M fext timeIncr us udots
