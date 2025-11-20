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

function [Kfd, Kfalfa, Khd, Khalfa, Fint] = framePlasticMatrices(E, Ks, A, l, uvector, npi, xpi, wpi, Mnp1, Cep_np1, Ghats, alfa)

  Kfd    = zeros(6, 6);
  Kfalfa = zeros(6, 1);
  Khd    = zeros(1, 6);
  Khalfa = 0;

  Fint   = zeros(6, 1);

  Bu = [-1 / l 1 / l];

  for ip = 1:npi

    N = bendingInterFuns (xpi(ip), l, 2);

    Bv = [N(1) N(3)];
    Btheta = [N(2) N(4)];

    Bd = [Bu   0 0 0 0; ...
          0  0 Bv  Btheta];

    Kfdj     = Bd' * [E * A 0; 0 Cep_np1(ip)] * Bd;

    Kfalfaj  = Bd' * [E * A 0; 0 Cep_np1(ip)] * [0 Ghats(ip)]';

    Khdj     = [0 Ghats(ip)] * [E * A 0; 0 Cep_np1(ip)] * Bd;

    Khalfaj  = Ghats(ip) * Cep_np1(ip) * Ghats(ip);

    epsilon = Bu * uvector;

    Fi      = Bd' * [E * A * epsilon; Mnp1(ip)];

    % stiffness matrices / integration (Gauss-Lobatto)
    Kfd    = Kfd    + Kfdj    * wpi(ip);
    Kfalfa = Kfalfa + Kfalfaj * wpi(ip);
    Khd    = Khd    + Khdj    * wpi(ip);
    Khalfa = Khalfa + Khalfaj * wpi(ip);

    % internal forces / integration (Gauss-Lobatto)
    Fint = Fint + Fi * wpi(ip);

  end

  Khalfa = Khalfa + Ks;

end
