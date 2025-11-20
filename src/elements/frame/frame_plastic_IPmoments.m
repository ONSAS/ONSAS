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

function [Mnp1, tM, Ghats] = frame_plastic_IPmoments(E, Iy, vvector, thetavector, npi, xpi, xd, l, alfa, kp, wpi)

  tM = 0;
  Ghats   = zeros(npi, 1);
  Mnp1    = zeros(npi, 1);

  for ip = 1:npi

    N = bendingInterFuns (xpi(ip), l, 2);

    Bv = [N(1) N(3)];
    Btheta = [N(2) N(4)];

    Ghats(ip) = -1 / l * (1 + 3 * (1 - 2 * xd / l) * (1 - 2 * xpi(ip) / l));

    khatxpi = Bv * vvector + Btheta * thetavector + Ghats(ip) * alfa;
    kenxpi = khatxpi - kp(ip);

    % moments Mn at integration points, corresponding to time n + 1
    Mnp1(ip) = E * Iy * kenxpi;

    % tM calculated with the moments Mn corresponding to time n + 1
    tM = tM - Ghats(ip) * Mnp1(ip) * wpi(ip);

  end

end
