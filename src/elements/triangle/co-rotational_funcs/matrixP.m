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

function [P, A] = matrixP(a1, a2, a3, G1, G2, G3, flag_second_mod) % ok
  % Eq. (26) of 10.1016/j.cma.2006.10.006
  
  if flag_second_mod == 1
    % Matrix A
    Ai = zeros(5, 3);
    Ai(3, 1) = 1;
    Ai(4, 2) = 1;
    Ai(5, 3) = 1;
    % Identity matrix
    I = zeros(5, 6);
    I(1, 1) = 1;
    I(2, 2) = 1;
    I(3, 4) = 1;
    I(4, 5) = 1;
    I(5, 6) = 1;

    % Projector Matrix
    P = zeros(15, 18);
    % Node 1
    A1 = Ai;
    A1(1, 3) = -a1(2);
    A1(2, 3) =  a1(1);
    % Node 2
    A2 = Ai;
    A2(1, 3) = -a2(2);
    A2(2, 3) =   a2(1);
    % Node 3
    A3 = Ai;
    A3(1, 3) = -a3(2);
    A3(2, 3) =   a3(1);
  else
    % Matrix A
    Ai = zeros(6, 3);
    Ai(4, 1) = 1;
    Ai(5, 2) = 1;
    Ai(6, 3) = 1;
    % Identity matrix
    I = eye(6);

    % Projector Matrix
    P = zeros(18, 18);
    % Node 1
    A1 = Ai;
    A1(1, 3) = -a1(2);
    A1(2, 3) =  a1(1);
    %
    A1(3, 1) =  a1(2);
    A1(3, 2) = -a1(1);
    % Node 2
    A2 = Ai;
    A2(1, 3) = -a2(2);
    A2(2, 3) =   a2(1);
    %
    A2(3, 1) =  a2(2);
    A2(3, 2) = -a2(1);
    % Node 3
    A3 = Ai;
    A3(1, 3) = -a3(2);
    A3(2, 3) =   a3(1);
    %
    A3(3, 1) =  a3(2);
    A3(3, 2) = -a3(1);
  end
  P1 = [I - A1 * G1', -A1 * G2', -A1 * G3'];
  P2 = [-A2 * G1', I - A2 * G2', -A2 * G3'];
  P3 = [-A3 * G1', -A3 * G2', I - A3 * G3'];
  
  % A1_jv = [-skew(a1) ; eye(3)]
  % A2_jv = [-skew(a2) ; eye(3)]
  % A3_jv = [-skew(a3) ; eye(3)]

  P = [P1; P2; P3];
  A = [A1; A2; A3];
  G = [G1; G2; G3];

  % PP=P*P;
  % PP-P
  % stop
  
  % eye
  % A'*G
  % G'*A

  %zeros
  % P*A
  % stop
end