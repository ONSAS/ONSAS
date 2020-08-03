%
% functions for computation of Boundary Conditions Degrees of Freedom..
%

function [neumDofs, diridofs, KS] = computeBCDofs( nnodes, Conec, nelems, nodalSprings, elementsParamsMat )


neumDofs = zeros( 6*nnodes, 1 ) ; % maximum possible vector

% loop for construction of vector of dofs
for elemNum = 1:size(elementsParamsMat,1)

  elements = find( Conec( :, 4+2 ) == elemNum ) ;

  if length( elements ) > 0
    elemType = elementsParamsMat(elemNum,1)  
    
    [numNodes, dofsStep] = elementTypeInfo ( elemType ) ;
    
    nodes    = Conec( elements, 1:numNodes) ;
    dofs     = nodes2dofs( nodes, 6)'       ;
    dofs     = dofs(1:dofsStep:end)         ;  
  
    neumDofs ( dofs ) = dofs ;
  end
end
% ----------------------------------------------------------------------


% ----------------------
fixeddofs = [] ;
KS        = sparse( 6*nnodes, 6*nnodes );  

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


neumDofs( diridofs ) = 0 ;
neumDofs = unique( neumDofs ) ;
if neumDofs(1) == 0,
  neumDofs(1)=[];
end
% -------------------------------------------------------------
