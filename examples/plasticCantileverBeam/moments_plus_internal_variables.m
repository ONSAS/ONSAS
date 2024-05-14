
% algorithm for the update of the internal variables for elastoplasticity with hardening
% and for the computation of the moment in the bulk Mn1

% displacements in time n + 1, dpn1 v1, v2, theta1, theta2, alpha, xd
% plastic curvature in time n / kappa_plas_n

function [kappa_plas_n1, xin11val, Mn1] = moments_plus_internal_variables( v1, v2, theta1, theta2 , xd, alpha, xin1, kappa_plas_n, Mc, My, kh1, kh2, E, Iy, l)

% integration points
x = [0 l/2 l] ;

Bv1 = -6/l^2*(1-2*x/l) ;
Bv2 =  6/l^2*(1-2*x/l) ;

Bt1 = -2/l*(2-3*x/l) ;
Bt2 = -2/l*(1-3*x/l) ;

G_bar = -(1+3*(1-2*xd/l)*(1-2*x/l))/l ;

% smooth curvature
kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alpha ;



kappa_plas_test = kappa_plas_n ;


Mn1_test = E*Iy*(kappa_bar - kappa_plas_test) ;

alpha

% ----------------------------------------
% softening 
if ( alpha >0) || (Mn1_test > Mc)
  hola

% ----------------------------------------
% hardening 
else

  % hardening function
  if xin1 <= (My - Mc)/kh1
    q = -kh1*xin1 ;
  else
    q = -(My - Mc)*(1-kh2/kh1) - kh2*xin1 ;
  end

  % yield function test
  phi_test = abs(Mn1_test)- (Mc - q) ;

  % 
  if phi_test <= 0  % 
  
      kappa_plas_n1 = kappa_plas_n ;
      xin11val = xin1 ;
      Mn1 = Mn1_test ;
  
  else
  
      % calculation of gamma_n1
      if xin1 + phi_test/(kh1 + E*Iy) <= (My - Mc)/kh1
  
          gamma_n1 = phi_test/(kh1 + E*Iy) ;
  
      else
  
          gamma_n1 = phi_test/(kh2 + E*Iy) ;
  
      end
  
      kappa_plas_n1 = kappa_plas_n + gamma_n1.*sign(Mn1_test) ;
      xin11val = xin1 + gamma_n1 ;
  
      Mn1 = E*Iy*(kappa_bar - kappa_plas_n1) ;
  
  end

end