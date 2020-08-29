% ---------------------------
% Toy Heat Transfer Problem
% ---------------------------

clear all, close all

addpath('../sources/')

nelem  = 100 ;
Tfinal = 1 ;
dt     = .0001 ;
Area   = 1 ;
rho    = 1 ;
cSpHe  = 1 ;
kCond  = 1 ;
Ltot   = 1 ;

%~ nt     = Tfinal / dt ;
nt     = 10 ;

nnodes = nelem +1 ;

xs = linspace(0,Ltot,nnodes )' ;
T0     = sin(pi*xs) + 0.5*sin(3*pi*xs) ;

diridofs = [ 1 nnodes ] ;
Tdiri  = 0 ;

neumdofs = 1:nnodes ;

neumdofs( diridofs ) = [] ;


lelem  = Ltot/nelem ;

Kdiffe = kCond * Area / lelem * [ 1 -1 ; -1 1 ] ;

MintEe = rho * cSpHe * Area * lelem / 6 * [ 2 1 ; 1 2 ] ;

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
%~ qext (end) = -1e1*Area ;

figure
hold on
plot( Ts(:,1) )  

for i=1:nt
  
  f = (    qext ( neumdofs    ) * dt ...
        - KdiffG( neumdofs, : ) * Ts( :, i ) * dt ...
        + MintEG( neumdofs, : ) * Ts( :, i ) 
      ) ;
      
  eig(CNN)
  Ts( neumdofs, i+1 ) = CNN \ f ;
  Ts( diridofs, i+1 ) = Tdiri ;
  
  plot( Ts(:,i+1) );  

end
