function [S, ConsMat] = cosseratSVK( consParams, Egreen, consMatFlag )

young   = consParams(1) ;
nu      = consParams(2) ;

lambda  = young * nu / ( (1 + nu) * (1 - 2*nu) ) ;
shear   = young      / ( 2 * (1 + nu) )          ;

S       = lambda * trace(Egreen) * eye(3)  +  2 * shear * Egreen ;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;
  
elseif consMatFlag == 1 % complex-step computation expression

  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratSVK', consParams, Egreen ) ;

elseif consMatFlag == 2  % analytical expression

  ConsMat = zeros(6,6);
  ConsMat (1,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ 1-nu , nu   , nu   ] ; 
  ConsMat (2,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , 1-nu , nu   ] ;
  ConsMat (3,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , nu   , 1-nu ] ;
  ConsMat (4,4  ) = shear ;
  ConsMat (5,5  ) = shear ;
  ConsMat (6,6  ) = shear ;

end
