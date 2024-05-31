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

% plastic hardening
% the standard trial-corrector (return mapping) algorithm is used

function [kp_np1, xi1_np1, Cep_np1] = plastic_hardening_step( E, Iy, xpi, xi1_n, kp_n, My, Mc, kh1, kh2, Ms)

kp_np1  = kp_n ;
xi1_np1 = xi1_n ;

npi       = length(xpi) ;
qs        = zeros(npi,1) ;
phis_test = zeros(npi,1) ;
Cep_np1   = zeros(npi,1) ;
Mnp1      = zeros(npi,1) ;


% disp(' /\  /\  /\  /\  /\  HARDENING /\  /\  /\  /\  /\') ;

for ip = 1:npi
  
  % yield criterion
  if xi1_n(ip) <= (My-Mc)/kh1

    qs(ip) = -kh1*xi1_n(ip) ;

  else

    qs(ip) = -(My-Mc)*(1-kh2/kh1)-kh2*xi1_n(ip) ;

  end

  phitest = abs(Ms(ip)) - (Mc - qs(ip)) ;
  phis_test(ip) = phitest ;


  % gamma values calculations (gamma derivative is the plastic multiplier)
  % the new values of internal variables are computed
  if phis_test(ip) <= 0 % elastic increment
    
      gamma = 0 ;
      
  else

    if (xi1_n(ip) + phis_test(ip)/(kh1+E*Iy)) <= (My-Mc)/kh1
        
        gamma = phis_test(ip)/(kh1+E*Iy) ;
    
    else
        
        gamma = phis_test(ip)/(kh2+E*Iy) ;
    
    end

    kp_np1(ip)  = kp_n(ip)  + gamma*sign(Ms(ip)) ;
    xi1_np1(ip) = xi1_n(ip) + gamma ;
  
  end

  % elastoplastic tangent bending modulus
  
  if gamma == 0
      
      Cep_np1(ip) = E*Iy ;

  elseif gamma > 0 && xi1_np1(ip) <= (My-Mc)/kh1
    
      Cep_np1(ip) = E*Iy*kh1/(E*Iy + kh1) ;

  elseif gamma > 0 && xi1_np1(ip) > (My-Mc)/kh1
    
      Cep_np1(ip) = E*Iy*kh2/(E*Iy + kh2) ;

  end

end

     % AA = kh1*xi1_n(1)+ Mc ;
     % BB = (xi1_n(1) + phis_test(1)/(kh1+E*Iy))*kh1 + Mc ;
     % CC = kh1*xi1_np1(1) + Mc ;