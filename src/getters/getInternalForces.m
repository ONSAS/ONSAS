% Copyright 2023, ONSAS Authors (see documentation)
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


% function that computes one value of each internal force for each element.
% the output is a nelemx6 matrix where the columns represent Fx Mx Fy My Fz Mz
% if the forces distribution is nonlinear, then the average is computed 

function Nx = getInternalForces( cellFint, mesh )

    for 
    cellFint


