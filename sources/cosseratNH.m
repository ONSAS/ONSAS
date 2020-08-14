function [S, ConsMat] = cosseratNH( consParams, Egreen, consMatFlag)

young   = consParams(1) ;
nu      = consParams(2) ;

lambda  = young * nu / ( (1 + nu) * (1 - 2*nu) ) ;
shear   = young      / ( 2 * (1 + nu) )          ;

C       = 2*Egreen + eye(3);  % Egreen = 1/2 (C - I)
invC    = inv(C);
detC    = det(C); % TODO use analyDet ?
J       = sqrt(detC);
S       = shear * (eye(3) - invC) + lambda * log(J) * invC;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;

elseif consMatFlag == 1 % complex-step computation expression

  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratNH', consParams, Egreen ) ;

else
  error("the analytical expression for the Neo-Hookean constitutive law is not available")

end
