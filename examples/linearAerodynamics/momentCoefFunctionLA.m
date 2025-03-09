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
function cMoment = momentCoefFunctionLA(betaRel, Re)
  % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
  % Emulation of aerodinamic moment extracted from the reference above
  if max(betaRel) < deg2rad(210)
    betaRel = -rad2deg(betaRel);
    cMoment = betaRel / 25  .* (betaRel < 25) +  (-(betaRel - 25) / 100 + 1) .* (betaRel > 25) + 0.15;
  end
end
