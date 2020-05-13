%
% functions for computation of Boundary Conditions Degrees of Freedom..
%

function [neumdofs, diridofs, KS] = computeNeumDiriDofs(nnodes, Conec, nelems, nodalSprings )

neumDofs = zeros( 6*nnodes, 1 ) ; % maximum possible vector

for elem = 1:nelems
  
  aux = nodes2dofs( Conec( elem, 1:4), 6)' ;
  
  switch Conec( elem, 7)
  case 1 % truss
    neumdofs ( aux(1:2:11) ) = aux(1:2:11) ;
  case 2 % beam
    neumdofs ( aux(1:11) ) = aux(1:11) ;
  case 3 % solid
    neumdofs ( aux(1:2:(6*4-1) ) ) = aux(1:2:(6*4-1)) ;
  end  
end
% ----------------------
fixeddofs = [] ;
KS      = sparse( 6*nnodes, 6*nnodes );  

for i=1:size(nodalSprings,1)
  aux = nodes2dofs ( nodalSprings (i,1) , 6 ) ;
  for k=1:6
    %
    if nodalSprings(i,k+1) == inf,
      fixeddofs = [ fixeddofs; aux(k) ] ;
    elseif nodalSprings(i,k+1) > 0,
      KS( aux(k), aux(k) ) = KS( aux(k), aux(k) ) + nodalSprings(i,k+1) ;
    end
  end
end

diridofs = unique( fixeddofs ) ; % remove repeated dofs
%~ diridofs = [ diridofs ; releasesDofs] ;


neumdofs( diridofs ) = 0 ;

neumdofs = unique( neumdofs ) ;
if neumdofs(1) == 0,
  neumdofs(1)=[];
end
% -------------------------------------------------------------
