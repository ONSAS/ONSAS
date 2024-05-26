% algorithm for the update of the internal variables for elastoplasticity with hardening
% and for the computation of the moment in the bulk Mn1

% displacements in time n + 1, dpn1 v1, v2, theta1, theta2, alpha, xd
% plastic curvature in time n / kappa_plas_n

function [kappa_plas_n1, xin11val, xin21val, alfan1, Mn1] = softHinge1DOF_semiAnalyticSol( v1, v2, theta1, theta2 , xd, alfan, xin1, xin2, kappa_plas_n, Mc, My, Mu, kh1, kh2, Ks, E, Iy, l)

% integration points
x = [0 l/2 l] ;
wp = [1/3 4/3 1/3]*l*0.5 ;

npi = length(x) ;

tM = 0 ;

Bv1 = -6/l^2*(1-2*x/l) ;
Bv2 =  6/l^2*(1-2*x/l) ;

Bt1 = -2/l*(2-3*x/l) ;
Bt2 = -2/l*(1-3*x/l) ;

G_bar = -(1+3*(1-2*xd/l)*(1-2*x/l))/l ;

kappa_plas_test = kappa_plas_n ;

% smooth curvature
kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alfan ;

Mn1 = E*Iy*(kappa_bar - kappa_plas_test) ;

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

% hardening

  % hardening function
  if xin1 <= (My - Mc)/kh1
    q = -kh1*xin1 ;
  else
    q = -(My - Mc)*(1-kh2/kh1) - kh2*xin1 ;
  end

  % yield function test
  phi_test = abs(Mn1)- (Mc - q) ;

  if phi_test <= 0
  
      kappa_plas_n1 = kappa_plas_n ;
      xin11val = xin1 ;
  
  else
  
      % calculation of gamma_n1
      if xin1 + phi_test/(kh1 + E*Iy) <= (My - Mc)/kh1
  
          gamma_n1 = phi_test/(kh1 + E*Iy) ;
  
      else
  
          gamma_n1 = phi_test/(kh2 + E*Iy) ;
  
      end
  
      kappa_plas_n1 = kappa_plas_n + gamma_n1.*sign(Mn1) ;
      xin11val = xin1 + gamma_n1 ;
  
  end
    
  alfan1      = alfan ;
  xin21val    = xin2 ;

  % smooth curvature
    kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alfan ;

    Mn1 = E*Iy*(kappa_bar - kappa_plas_n1) ;

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

% softening

for ii = 1:npi

if abs(Mn1(ii)) >= Mu

for ip = 1:npi

tM = tM - G_bar(ip)*Mn1(ip)*wp(ip) ;

end

qfailxpi = min(-Ks*xin2, Mu) ;

phifailxpi = abs(tM)-(Mu-qfailxpi) ;

if phifailxpi <= 0

    alfan1      = alfan ;
    xin21val    = xin2 ;

else
    
    if  xin2 <= -Mu/Ks

        gamma2 = phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks) ;

    else

        gamma2 = abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)) ;
    
    end

    alfan1      = alfan     + gamma2*sign(tM) ;
    xin21val    = xin2      + gamma2 ;

end

kappa_plas_n1 = kappa_plas_n ;
xin11val = xin1 ;

% smooth curvature

kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alfan1 ;
 
Mn1 = E*Iy*(kappa_bar - kappa_plas_n1) ;

end

end