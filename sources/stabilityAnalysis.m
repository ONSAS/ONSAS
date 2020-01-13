
% --- nonlinear buckling analysis as in section 6.8.2 from Bathe, FEM Procedures 2nd edition. ---

function [ factor_crit, nKeigpos, nKeigneg] = stabilityAnalysis ( KTtm1red, KTtred, currLoadFactor, nextLoadFactor );  

  [a,b] = eig( KTtred ) ;
  Keigvals = diag(b) ; 
  nKeigpos = length( find(Keigvals >  0 ) ) ;
  nKeigneg = length( find(Keigvals <= 0 ) ) ;

  [vecgamma, gammas ] = eig( KTtred, KTtm1red ) ;
  
  gammas = diag( gammas);
 
  if length( find( gammas >  0 ) ) > 0,
  
    gamma_crit  = min ( gammas ( find( gammas >  0 ) ) ) ;
    if gamma_crit ~= 1 
      lambda_crit = 1 / ( 1 - gamma_crit )  ;               
      factor_crit = currLoadFactor + lambda_crit * (nextLoadFactor - currLoadFactor) ;
    else
      factor_crit = 0 ;
    end
  else
    factor_crit = 0;
  end
