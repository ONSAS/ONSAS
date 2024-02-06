% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan

% =========================================================================

function [Kfd, Kfalfa, Khd, Khalfa, kpn1xpi, xin11xpi, xin21xpi, M1xpi, tM, xd, Fi, alfan1] = integrand_plastic(jj, xpi, xd, l, uvector, vvector, thetavector, alfan, xin1, kpn, kpn1, E, Iy, My, Mc, kh1, kh2, A, Ks, xin2, Mu, Cep,tM)

% elastoplasticity with hardening
% the usual trial-corrector (return mapping) algorithm
% is used at each of the Gauss–Lobatto integration points.

Bu = [-1/l 1/l] ;

N = bendingInterFuns (xpi, l, 2) ;
Bv = [N(1) N(3)] ;
Btheta = [N(2) N(4)] ;

Bd = [ Bu  0 0 0 0    ; ...
       0 0 Bv  Btheta ] ;

Ghat = -1/l * ( 1 + 3*(1-2*xd/l)*(1-2*xpi/l) ) ;

% curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)

khat = Bv*vvector + Btheta*thetavector + Ghat*alfan ;

% khat2 = dirac(xd)*alfan ;
% kn = khat + khat2 ;
ken = khat - kpn(jj) ;

% moment

Mxpi = E*Iy*ken ;

% yield criterion
if xin1(jj) <= (My-Mc)/kh1
  qxpi = -kh1*xin1(jj) ;
        
else
  qxpi = -(My-Mc)*(1-kh2/kh1)-kh2*xin1(jj) ;

end

phixpi = abs(Mxpi) - (Mc - qxpi) ;

% test values

% kpn1test = kpn ;
% xin11test = xin1 ;
phitest =  phixpi ;

% gamma values calculations (gamma derivative is the plastic multiplier)
% the new values of internal variables are computed

if phitest <= 0 % elastic increment

    gamma    = 0 ;
    kpn1xpi  = kpn(jj) ;
    xin11xpi = xin1(jj) ;
    M1xpi    = Mxpi ;

else

    if xin1(jj) + phitest/(kh1+E*Iy)<=(My-Mc)/kh1

        gamma = phitest/(kh1+E*Iy) ;

    else

        gamma = phitest/(kh2+E*Iy) ;
    
    end
    
    kpn1xpi     = kpn(jj) + gamma*sign(Mxpi) ;
    xin11xpi    = xin1(jj) + gamma ;
    M1xpi       = E*Iy*(khat-kpn1xpi) ;

end

% elastoplastic tangent bending modulus

if      gamma == 0
        Cep = E*Iy ;

elseif  gamma > 0 && xin11xpi <= (My-Mc)/kh1
        Cep = E*Iy*kh1/(E*Iy + kh1) ;

elseif  gamma > 0 && xin11xpi > (My-Mc)/kh1
        Cep = E*Iy*kh2/(E*Iy + kh2) ;

end

% stiffness matrices

Kfd     = Bd'*[E*A 0; 0 Cep]*Bd ;

Kfalfa  = Bd'*[E*A 0; 0 Cep]*[0 Ghat]' ;

Khd     = [0 Ghat]*[E*A 0; 0 Cep]*Bd ;

Khalfa  = Ghat*Cep*Ghat ;

epsilon = Bu*uvector ;

Fi      = Bd' * [E*A*epsilon; M1xpi] ;

% plastic softening at the discontinuity
% the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity

% test values

% alfan1test = alfan ;
% xin2test = xin2 ;

% softening criterion (failure function) at integration points

qfailxpi = min(-Ks*xin2(jj), Mu) ;

phifailxpi = abs(tM)-(Mu-qfailxpi) ;

if phifailxpi <= 0
    alfan1 = alfan ;
    xin21xpi = xin2(jj) ;
else

    if  xin2(jj)<=-Mu/Ks

        gamma2 = phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks) ;

    else

        gamma2 = abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)) ;
    
    end
    
    alfan1      = alfan + gamma2*sign(tM) ;
    xin21xpi    = xin2(jj) + gamma2 ;
    xd          = xpi ;       

end