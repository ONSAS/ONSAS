function v = mat2voigt( Tensor, factor )
  
  %~ if norm( Tensor - Tensor' ) > 1e-10
    %~ Tensor
    %~ norm( Tensor - Tensor' )
    %~ error('tensor not symmetric')
  %~ end
  
  v = [ Tensor(1,1)  Tensor(2,2)  Tensor(3,3)  factor*Tensor(2,3) factor*Tensor(1,3) factor*Tensor(1,2) ]' ;
