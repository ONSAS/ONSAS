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

%
function [Gt] = matrixGt_jv(q1g, q2g, Rr, l)

    qg  = (q1g + q2g) / 2;

    % local coords
    q  = Rr' *  qg;
    q1 = Rr' * q1g;
    q2 = Rr' * q2g;

    nu   = q(1) / q(2);
    nu11 = q1(1) / q(2);
    nu12 = q1(2) / q(2);

    % nu21_jv = q2(1) / q(2);
    % nu22_jv = q2(2) / q(2);

    nu21 = 2 * nu - nu11;
    nu22 = 2 - nu12;

    % [ nu21_jv nu21 nu22_jv nu22 ]

    Gt = [  0   0    nu/l  nu12/2  -nu11/2  0  0    0    -nu/l     nu22/2  -nu21/2  0 ;...
            0   0    1/l     0        0     0  0    0    -1/l       0        0      0 ;...
            0  -1/l  0       0        0     0  0  1/l       0       0        0      0   ];



end