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
function [T, x02, x03, y03] = edgeLocalAxisShellTriangle(p1, p2, p3)
  % Calculates the matrix for transformation of basis between global and local axis;
  % p1, p2 and p3 are the position vector for the nodes in global coordinates;
  % the local x axis is paralel to the side connecting nodes 1 and 2
  % the local z axis is normal to the element plane;
  % the origin of local axis is located at node 1

  p12 = p2 - p1;
  p13 = p3 - p1;

  au_zl =  cross(p12, p13);
  u_zl = au_zl / norm(au_zl);

  x02 = norm(p12);
  u_xl = p12 / x02;

  u_yl = cross(u_zl, u_xl);

  T = [u_xl u_yl u_zl]';

  x03 = dot(u_xl, p13);
  y03 = dot(u_yl, p13);
end