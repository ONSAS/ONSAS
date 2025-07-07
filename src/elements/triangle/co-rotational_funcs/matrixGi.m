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

function [G1, G2, G3] = matrixGi(a1, a2, a3, r1, r2, r3, e1_flag, origin_flag) % ok

  a21 = a2 - a1;
  a31 = a3 - a1;
  a32 = a3 - a2;
  a13 = a1 - a3;

  % if origin_flag == 0
    % Eq. (27) of 10.1016/j.cma.2006.10.006
    
    d = norm(a21);
    v = norm(cross(a21, a31));
    c = a1' * r1 + a2' * r2 + a3' * r3;

    G1 = zeros(6, 3);
    G1(3, 1) = a32(1) / v;
    G1(3, 2) = a32(2) / v;

    G2 = zeros(6, 3);
    G2(3, 1) = a13(1) / v;
    G2(3, 2) = a13(2) / v;
    
    G3 = zeros(6, 3);
    G3(3, 1) = a21(1) / v;
    G3(3, 2) = a21(2) / v;
    % G3(3, 1) = a32(1) / v;
    % G3(3, 2) = a32(2) / v;

    % e1_flag
    if e1_flag == 0
      %
      G1(1, 3) =  -r1(2) / c;
      G1(2, 3) =   r1(1) / c;
      %
      G2(1, 3) =  -r2(2) / c;
      G2(2, 3) =   r2(1) / c;
      %
      G3(1, 3) =  -r3(2) / c;
      G3(2, 3) =   r3(1) / c;
    
    else
      G1(1, 3) = 0;
      G1(2, 3) = -1/d;
      %
      G2(1, 3) = 0;
      G2(2, 3) = 1/d;
      %
      % G3(1, 3) = 0;
      % G3(2, 3) = 0;    
    end

  % elseif origin_flag == 1

  %   G1 = zeros(6,3);
  %   G1(3,1) = a32(1)/(a3(2)*a2(1));
  %   G1(3,2) = 1/a2(1);
  %   G1(2,3) = -1/(a2(1));
        
  %   G2 = zeros(6,3);
  %   G2(3,1) = -a3(1)/(a3(2)*a2(1));
  %   G2(3,2) = -1/a2(1);
  %   G2(2,3) = 1/(a2(1));

  %   G3 = zeros(6,3);
  %   G3(3,1) = -1/(a3(2));

  % end

end