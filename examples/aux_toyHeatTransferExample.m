% ---------------------------
% Toy Heat Transfer Problem
% ---------------------------

clear all, close all

addpath('../sources/')

nelem  = 1 ;
Tfinal = .0003 ;
dt     = .0001 ;
Area   = 1 ;
rho    = 1 ;
cSpHe  = 1 ;
kCond  = 1 ;
Ltot   = 2 ;

Temp0  = 10 ;

diridofs = 1 ;
Tdiri  = 10 ;


nnodes = nelem +1 ;

neumdofs = 1:nnodes ;

neumdofs( diridofs ) = [] ;

nt     = Tfinal / dt ;

lelem  = Ltot/nelem ;

Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;

MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;

T0     = ones( nnodes, 1)*Temp0 ;

Ts     = T0 ;

KdiffG = zeros( nnodes, nnodes ) ;
MintEG = zeros( nnodes, nnodes ) ;

% assembly
for i = 1 : nelem
  nodeselem = [ i i+1 ] ;
  
  elemDofs = nodes2dofs ( nodeselem, 1 ) ;
  
  KdiffG( elemDofs , elemDofs ) = ...
    KdiffG( elemDofs , elemDofs ) + Kdiffe  ; 
  
  MintEG( elemDofs , elemDofs ) = MintEG( elemDofs , elemDofs ) + MintEe  ; 
end

CDD = MintEG(diridofs, diridofs) ;
CND = MintEG(neumdofs, diridofs) ;
CNN = MintEG(neumdofs, neumdofs) ;

qext = zeros(nnodes,1) ;

%~ qext (1) = 1*Area ;
qext (end) = -1e6*Area ;

figure
hold on
plot( Ts(:,1) )  


for i=1:nt
  qext ( neumdofs    ) * dt
  - KdiffG( neumdofs, : ) * Ts( :, i ) * dt
  + MintEG( neumdofs, : ) * Ts( :, i ) 
  
  f = (    qext ( neumdofs    ) * dt ...
        - KdiffG( neumdofs, : ) * Ts( :, i ) * dt ...
        + MintEG( neumdofs, : ) * Ts( :, i ) 
      )
  
  Ts( neumdofs, i+1 ) = CNN \ f ;
  Ts( diridofs, i+1 ) = Tdiri ;
  
  plot( Ts(:,i+1) );  

end
