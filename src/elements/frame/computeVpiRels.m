% Copyright 2025, ONSAS Authors (see documentation)
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
% This function returns the relative projected velocity in local coordinates

function [VpiRel, VpiRelPerp, VrelG] = computeVpiRels(udotFlow, udotFrame, Rroof, Rr, L2, L3)
  % the relative velocity in global cooridantes is:
  VrelG = udotFlow - udotFrame;
  % then the projection (in t2,t3 plane) of the relative flow velocity in deformed coordinates is:
  VpiRel = L2 * Rroof' * Rr' * VrelG;
  % the perpendicular flow relative velocity projection in deformed coordinates is:
  VpiRelPerp = L3 * VpiRel;
end
