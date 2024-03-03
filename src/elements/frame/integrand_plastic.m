% Copyright 2024, ONSAS Authors (see documentation)
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan

% =========================================================================

function [soft_hinge_boolean, Kfd, Kfalfa, Khd, Khalfa, kpn1xpi, xin11xpi, M1xpi, xd, Fi] ...
  = integrand_plastic(soft_hinge_boolean, jj, xpi, xd, l, A, uvector, vvector, ...
                      thetavector, alfan, xin1, kpn, E, Iy, My, Mc, kh1, kh2, Cep)

% elastoplasticity with hardening
% the usual trial-corrector (return mapping) algorithm
% is used at each of the Gauss–Lobatto integration points
Bu = [-1/l 1/l] ;

N = bendingInterFuns (xpi, l, 2) ;

Bv = [N(1) N(3)] ;
Btheta = [N(2) N(4)] ;

Bd = [ Bu  0 0 0 0    ; ...
       0 0 Bv  Btheta ] ;

Ghat = -1/l * ( 1 + 3*(1-2*xd/l)*(1-2*xpi/l) ) ;

% curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)

if soft_hinge_boolean == true

thetavector(2) = thetavector(1) + alfan ;
vvector(2) = vvector(1) + xd*thetavector(1) + (l-xd)*(alfan + thetavector(1)) ;

end

khatxpi = Bv*vvector + Btheta*thetavector + Ghat*alfan ;

% khat2 = dirac(xd)*alfan ;
% kn = khat + khat2 ;
kenxpi = khatxpi - kpn(jj) ;

% moment
Mxpi = E*Iy*kenxpi ;

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

        kpn1xpi     = kpn(jj)   + gamma*sign(Mxpi) ;
        xin11xpi    = xin1(jj)  + gamma ;

    else

        gamma = phitest/(kh2+E*Iy) ;

        kpn1xpi     = kpn(jj)   + gamma*sign(Mxpi) ;
        xin11xpi    = xin1(jj)  + gamma ;
    
    end

    if soft_hinge_boolean == true && xpi == xd

        kpn1xpi  = kpn(jj) ;
        xin11xpi = xin1(jj) ;

    end

    % curvatures (time n + 1) / khat1 (continuous part of the curvature), khat2 (localized part of the curvature)
    
    if soft_hinge_boolean == true

        thetavector(2) = thetavector(1) + alfan ;
        vvector(2) = vvector(1) + xd*thetavector(1) + (l-xd)*(alfan + thetavector(1)) ;

    end
    
    khat1xpi = Bv*vvector + Btheta*thetavector + Ghat*alfan ;

    M1xpi = E*Iy*(khat1xpi-kpn1xpi) ;

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

end