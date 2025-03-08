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
function C_d = dragCoefFunctionLA(betaRel, Re)
  % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
  % Emulation of drag extrated from the reference above
  if deg2rad(-210) < betaRel && betaRel < deg2rad(210)
    omega_drag = 0.8 * 2 * pi / deg2rad(180);
    C_d = +1.2 - 0.6 * cos(omega_drag * betaRel);
  else
    betaRel;
    error ('The angle must be smaller than 210 to use this aerodynamic coefficent');
  end

end
