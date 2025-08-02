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
%
function [Kel, fintLocCoord, Kb] = localShellTriangle(x02, x03, y03, E, nu, h, Ul, flag_OPT)

  % Calculate the area of the triangle
  area = x02 * y03 / 2;

  % membrane stiffness
  if flag_OPT == 0
	aux1 = h *  E / (1 - nu^2);
	aux2 = nu * aux1;
	Dm = [[aux1, aux2, 0]; [aux2, aux1, 0]; [0, 0, aux1 * (1 - nu) / 2]];
	Bm = cstB(x02, x03, y03);
	Km = area * Bm' * Dm * Bm;
	im = [1, 2, 7, 8, 13, 14];
	% membrane forces (constant)
	Ulm = Ul(im);
	N = Dm * Bm * Ulm;
  else
    Km = opt_membrane_element() ;
  end

  % bending stiffness
  aux1 = E * h^3 / (12 * (1 - nu^2));
  aux2 = nu * aux1;
  Db = [[aux1, aux2, 0]; [aux2, aux1, 0]; [0, 0, aux1 * (1 - nu) / 2]];
  ib = [3, 4, 5, 9, 10, 11, 15, 16, 17];
  int_point = [[1 ./ 6. 1 ./ 6.]; [2 ./ 3., 1 ./ 6.]; [1 ./ 6., 2 ./ 3.]];
  Kb = zeros(9, 9);
  wgt = area / 3.0;
  for ipt = 1:3
    psi = int_point(ipt, 1);
    eta = int_point(ipt, 2);
    Bb = dktB(psi, eta, x02, x03, y03);
    Kb = Kb + wgt * Bb' * Db * Bb;
  end
  % plate moments (calculated at the center of the element)
  Ulb = Ul(ib);
  Bb = dktB(1 ./ 3., 1 ./ 3., x02, x03, y03);
  curv = Bb * Ulb;
  M = Db * curv;

  % assembling the stiffness matrix of the shell element in local coordinates
  Kel = zeros(18, 18);
  Kel(im, im) = Km;
  Kel(ib, ib) = Kb;

  k_dr = min(min(abs(Kb))) * 1.e-4;
  % k_dr = max(k_dr, 1.e-4);
  % k_dr=0;
  Kel(6, 6) = k_dr;
  Kel(12, 12) = k_dr;
  Kel(18, 18) = k_dr;

  Kel = real(Kel);

  % returning moments in local element coordiantes
  fintLocCoord = zeros(1, 3);
  fintLocCoord = fintLocCoord + M;  % Store the moments in fintLocCoord
end
