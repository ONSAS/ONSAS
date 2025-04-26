% =========================================================================

% Euler-Bernoulli element with embeded discontinuity

% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan
% =========================================================================

% Solution of local equations

% plastic hardening
% the standard trial-corrector (return mapping) algorithm is used

function [kp_np1, xi1_np1, Cep_np1, alfan1, xin21, xdn] = localsystem( E, Iy, xpi, xi1_n, kp_n, My, Mc, kh1, kh2, Ms, xd, alfan, xin2, tM, l, Mu, Ks, SH_boole_np1)

if SH_boole_np1 == false

kp_np1  = kp_n ;
xi1_np1 = xi1_n ;

npi       = length(xpi) ;
qs        = zeros(npi,1) ;
phi_np1   = zeros(npi,1) ;
Cep_np1   = zeros(npi,1) ;
alfan1 = 0 ;
xin21 = 0 ;
xdn = xd ;

for ip = 1:npi

  % yield criterion
  if xi1_np1(ip) <= (My-Mc)/kh1

    qs(ip) = -kh1*xi1_np1(ip) ;

  else

    qs(ip) = -(My-Mc)*(1-kh2/kh1)-kh2*xi1_np1(ip) ;

  end

  phi_np1(ip) = abs(Ms(ip)) - (Mc - qs(ip)) ;


  % gamma_val values calculations (gamma_val derivative is the plastic multiplier)
  % the new values of internal variables are computed
  if phi_np1(ip) <= 0 % elastic increment
    
      gamma_val = 0 ;
      
  else

    if (xi1_n(ip) + phi_np1(ip)/(kh1+E*Iy)) <= (My-Mc)/kh1
        
        gamma_val = phi_np1(ip)/(kh1+E*Iy) ;
    
    else
        
        gamma_val = phi_np1(ip)/(kh2+E*Iy) ;
    
    end

    kp_np1(ip)  = kp_n(ip)  + gamma_val*sign(Ms(ip)) ;
    xi1_np1(ip) = xi1_n(ip) + gamma_val ;
  
  end

  % elastoplastic tangent bending modulus
  
  if gamma_val == 0
      
      Cep_np1(ip) = E*Iy ;

  elseif gamma_val > 0 && xi1_np1(ip) <= (My-Mc)/kh1
    
      Cep_np1(ip) = E*Iy*kh1/(E*Iy + kh1) ;

  elseif gamma_val > 0 && xi1_np1(ip) > (My-Mc)/kh1
    
      Cep_np1(ip) = E*Iy*kh2/(E*Iy + kh2) ;

  end

end

else

% plastic softening at the discontinuity
% the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity
% softening criterion (failure function) at integration points

Cep_np1 = ones(3,1)*E*Iy ;

kp_np1  = kp_n ;
xi1_np1 = xi1_n ;
xdn = xd ;

qfailxpi = min(-Ks*xin2, Mu) ; % test

phifailxpi = abs(tM)-(Mu-qfailxpi) ;

if phifailxpi <= 0

    alfan1 = alfan ;
    xin21 = xin2 ;

else

    integral_conG = (4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2) ;
    
    gamma_tent = phifailxpi / ( integral_conG + Ks) ;

    if  (xin2 + gamma_tent) <= -Mu/Ks

        gamma2 = gamma_tent ;

    else

        gamma2 = abs(tM) / integral_conG ;
    
    end

    alfan1      = alfan     + gamma2*sign(tM) ;
    xin21       = xin2      + gamma2 ;

end

end