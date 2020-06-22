function normal = computeNormalForces(Conec, Stresses, crossSecsParams )

nElems = size(Conec,1) ;

normal = Stresses(:,1) .* crossSecsParams( Conec(:, 6 ) , 1) ... 
  .* (( Conec(:,6)==1) + (Conec(:,6)==2) ) ;
