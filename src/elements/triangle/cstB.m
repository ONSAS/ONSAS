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
function [ B ] = cstB(x02, x03, y03)
    %calculate the strain-displacement matrix for the constant stress triangular element (CST)
    %x02, x03 and y03 are the local coordinates of the nodes 2 and 3 

    area02 = x02*y03;

    %( with y01 = y02 = x01 = 0)
    %bi = (yj - yk )/2A
    b1 = - y03/area02;
    b2 =   y03/area02;
    b3 =  0;

    %ci = (xk - xj)/2A
    c1 = (x03 - x02)/area02;
    c2 = - x03/area02;
    c3 = x02/area02;

    B = [   [b1,  0, b2,  0, b3,  0];
            [ 0, c1,  0, c2,  0, c3];
            [c1, b1, c2, b2, c3, b3]; ];

end
