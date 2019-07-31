% function for testing experimental feature.

function fextAngSpr = loadsAngleSpring( Nodes, Conec, bendStiff )

n = size(Nodes,1);

fextAngSpr = zeros( 6*n,1);

for i = 1:n
  %
  
  if bendStiff(i) > 0  
    %
    if Conec(i-1  ,2) ~= Conec(i,1);
      error('conectivity error for angle springs.');
    end
    
    EI = bendStiff (i);
    
    indxim1 = Conec(i-1,1) ;
    indxi   = Conec(i-1,2) ;
    indxip1 = Conec(i  ,2) ;
    
    ti   = Nodes(indxi  ,:) - Nodes(indxim1,:) ;
    tip1 = Nodes(indxip1,:) - Nodes(indxi  ,:) ;
  
    li   = norm( ti  ) ;
    lip1 = norm( tip1) ;
    
    dtds    = ( tip1/lip1 - ti/li ) / ( lip1/2 + li/2 ) ;
  
    Fiim1 = +dtds *  EI  / norm( cross( ti  , dtds) )  ;
    Fim1i = -dtds *  EI  / norm( cross( ti  , dtds) )  ;
    
    Fiip1 = +dtds *  EI  / norm( cross( tip1, dtds) )  ;
    Fip1i = -dtds *  EI  / norm( cross( tip1, dtds) )  ;
  
    aux = nodes2dofs( indxim1,6) ;
    fextAngSpr( aux(1:2:6) ) = Fim1i ; 
  
    aux = nodes2dofs( indxi  ,6) ;
    fextAngSpr( aux(1:2:6) ) = Fiim1 + Fiip1 ; 
  
    aux = nodes2dofs( indxip1,6) ;
    fextAngSpr( aux(1:2:6) ) = Fip1i ; 

  end    
end
