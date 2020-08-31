% ---------------------------
% Toy Heat Transfer Problem
% ---------------------------

function aux_toyHeatTransferExample(caseNum, nelem, plotBoolean )

addpath('../sources/')

dt     = 0.001 ;
rho    = 1 ;
cSpHe  = 1 ;
kCond  = 1 ;
Ltot   = 1 ; % domain [0,1]
Area   = 1 ;

%~ nt     = 3 ; % defined if stopped before Tfinal


switch caseNum

case 1  % diri-diri conds

  Tfinal      = 0.5   ;
  %~ nelem       = 5     ;
  diridofs    = [ 1 nelem+1 ] ;
  Tdiri       = 0 ;
  anlyBoolean = 1;
  
case 2 % diri-neum conds

  %~ nelem  = 10     ;
  Tfinal = 1   ;
  diridofs = [ 1 ] ;
  Tdiri  = 0 ;
  anlyBoolean = 0;
  qentr = -0 ;

case 3 % diri-robin conds

  nelem  = 10     ;
  Tfinal = 1   ;
  diridofs = [ 1 ] ;
  Tdiri  = 0 ;
  anlyBoolean = 0;
  qentr = -0 ;
   
end



if exist('nt')==0
  nt     = Tfinal / dt ;
end

alpha = kCond / ( rho * cSpHe ) ;

nnodes = nelem +1 ;

nnodesAnly = 100 ;

xs     = linspace(0, Ltot, nnodes     )' ;
xsAnly = linspace(0, Ltot, nnodesAnly )' ;


neumdofs = 1:nnodes ;
neumdofs( diridofs ) = [] ;

lelem  = Ltot/nelem ;

% local elemental diffussion equation
Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;

MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;

% initial temperature
T0     = sin(pi*xs) + 0.5*sin(3*pi*xs) ;
Ts     = T0 ;

   
% ------------------------
% matrices assembly
KdiffG = zeros( nnodes, nnodes ) ;
MintEG = zeros( nnodes, nnodes ) ;

for i = 1 : nelem
  nodeselem = [ i i+1 ] ;
  
  elemDofs = nodes2dofs ( nodeselem, 1 ) ;
  
  KdiffG( elemDofs , elemDofs ) = ...
  KdiffG( elemDofs , elemDofs ) + Kdiffe  ; 
  
  MintEG( elemDofs , elemDofs ) = ...
  MintEG( elemDofs , elemDofs ) + MintEe  ; 

end
% ------------------------

CDD = MintEG(diridofs, diridofs) ;
CND = MintEG(neumdofs, diridofs) ;
CNN = MintEG(neumdofs, neumdofs) ;

qext = zeros( nnodes, 1 ) ;
if caseNum == 2
  qext(end) = qentr ;
end

if plotBoolean
figure
hold on, grid on
end

MS = 10 ; 
LW = 1.5 ; 

% ------------------------
% time loop
for i=0:nt
  t = i*dt ;
  
  if i~=0
    f = (    qext ( neumdofs    ) * dt ...
          - KdiffG( neumdofs, : ) * Ts( :, i ) * dt ...
          + MintEG( neumdofs, : ) * Ts( :, i ) 
        ) ;
        
    Ts( neumdofs, i+1 ) = CNN \ f ;
    Ts( diridofs, i+1 ) = Tdiri   ;
  end  

  if anlyBoolean
    TsAnly (:,i+1) = exp(-(  pi*alpha)^2 * t ) *       sin(     pi * xsAnly ) ...
                   + exp(-(3*pi*alpha)^2 * t ) * 0.5 * sin( 3 * pi * xsAnly ) ;
  end
  
  % --- plots ---
  if plotBoolean
  
    if mod(i, round(nt/10) )==0
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
  figure
  plot( Ts(2,:) )
end

KdiffGNN = KdiffG(neumdofs, neumdofs ) ;

save -mat auxVars.mat KdiffGNN CNN Ts
