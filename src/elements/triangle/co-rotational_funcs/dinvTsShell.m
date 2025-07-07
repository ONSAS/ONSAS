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
function [dinv_Ts] = dinvTsShell(t, v)

  % t - local thetas (deformed reference frame)
  % v - moments

  a = norm(t) ;

  if a == 0
      dinv_Ts = -0.5 * skew(v);
  else
      
      eta = (2 * sin(a) - a * (1+cos(a))) / (2 * a^2 * sin(a));
      mu  = (a * (a + sin(a)) - 8 * sin(a/2)^2) / (4 * a^4 * sin(a/2)^2);
      I = eye(3);
      dinv_Ts = eta * (t * v' -2 * v * t' + (t' * v) * I) + mu * (skew(t) * skew(t)) * (v * t') -0.5 * skew(v);
  end