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

% plastic softening at the discontinuity
% the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity
% softening criterion (failure function) at integration points

function [alfan1, xin21, xd] ...
  = soft_hinge(xd, alfan, xin2, tM, l, E, Iy, Mu, Ks)

qfailxpi = min(-Ks*xin2, Mu) ;

phifailxpi = abs(tM)-(Mu-qfailxpi) ;

if phifailxpi <= 0

    alfan1 = alfan ;
    xin21 = xin2 ;

else
    
    if  xin2 <= -Mu/Ks

        gamma2 = phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks) ;

    else

        gamma2 = abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)) ;
    
    end

    alfan1      = alfan     + gamma2*sign(tM) ;
    
    xin21       = xin2      + gamma2 ;

end