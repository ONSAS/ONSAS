% ---------------------------
% Toy Heat Transfer Problem
% ---------------------------

clear all, close all

addpath('../sources/')

nelem  = 5 ;
Tfinal = .5 ;
dt     = .001 ;
rho    = 1 ;
cSpHe  = 1 ;
kCond  = 1 ;
Ltot   = 1 ; % domain [0,1]
Area   = 1 ;

%~ nt     = 3 ; % defined if stopped before Tfinal

if exist('nt')==0
  nt     = Tfinal / dt ;
end

alpha = kCond / ( rho * cSpHe ) ;

nnodes = nelem +1 ;

nnodesAnly = 100 ;

xs     = linspace(0, Ltot, nnodes     )' ;
xsAnly = linspace(0, Ltot, nnodesAnly )' ;

diridofs = [ 1 nnodes ] ;
Tdiri  = 0 ;

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

figure
hold on, grid on

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

  TsAnly (:,i+1) = exp(-(  pi*alpha)^2 * t ) *       sin(     pi * xsAnly ) ...
                 + exp(-(3*pi*alpha)^2 * t ) * 0.5 * sin( 3 * pi * xsAnly ) ;
  
  if mod(i, 10)==0
    plot( xs    , Ts    (:, i+1), 'b-o', 'markersize', MS,'linewidth',LW );  
    plot( xsAnly, TsAnly(:, i+1), 'r--'  , 'markersize', MS,'linewidth',LW );  
  end
end
% ------------------------

figure
plot( Ts(2,:) )


KdiffGNN = KdiffG(neumdofs, neumdofs ) ;

save -mat auxVars.mat KdiffGNN CNN
