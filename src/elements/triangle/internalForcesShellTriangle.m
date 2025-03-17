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
% Implementation of a triangular finite element with 6 dfos (3 translations and 3 rotations) per node for the analysis of linear elastic isotropic shells with constant thickness.
% The element is formed by the superposition of a plate element (DKT) and a plane stress element (CST) with
% with addition to artificial drilling (rotation about the axis normal to the element plane) stiffness.
%
function [fs, ks, fintLocCoord] = internalForcesShellTriangle(elemCoords, elemDisps, modelName, modelParams, thickness)

  % material and geometric parameters
  young_modulus = modelParams(1);
  poisson_ratio = modelParams(2);
  h = thickness;

  elemCoords = elemCoords';

  r1g = elemCoords(1:3);
  r2g = elemCoords(4:6);
  r3g = elemCoords(7:9);
  rog = (r1g + r2g + r3g) / 3;

  Ug = switchToTypeIndexing(elemDisps);

  u1g = Ug(1:3);
  q1  = Ug(4:6);
  u2g = Ug(7:9);
  q2  = Ug(10:12);
  u3g = Ug(13:15);
  q3  = Ug(16:18);

  p1g = r1g + u1g;
  p2g = r2g + u2g;
  p3g = r3g + u3g;
  pog = (p1g + p2g + p3g) / 3;

  % eq. (35) of 10.1016/j.cma.2006.10.006
  R1g = globalRotationMatrix(q1);
  R2g = globalRotationMatrix(q2);
  R3g = globalRotationMatrix(q3);

  [To, x02, x03, y03] = edgeLocalAxisShellTriangle(r1g, r2g, r3g);
  [Tr, x02, x03, y03] = edgeLocalAxisShellTriangle(p1g, p2g, p3g);

  Ro = To';
  Rr = Tr';

  r1o = To * (r1g - rog);
  r2o = To * (r2g - rog);
  r3o = To * (r3g - rog);

  % eq. (1) of 10.1016/j.cma.2006.10.006
  u1def = Tr * (p1g - pog) - r1o;
  u2def = Tr * (p2g - pog) - r2o;
  u3def = Tr * (p3g - pog) - r3o;

  % eq. (2) of 10.1016/j.cma.2006.10.006
  R1def = Tr * R1g * Ro;
  R2def = Tr * R2g * Ro;
  R3def = Tr * R3g * Ro;

  % eq. (13) of 10.1016/j.cma.2006.10.006
  v1def = rotationVector(R1def);
  v2def = rotationVector(R2def);
  v3def = rotationVector(R3def);

  % eq. (12) of 10.1016/j.cma.2006.10.006
  pl = zeros(15, 1);
  pl(1:2) = u1def(1:2);
  pl(3:5) = v1def;
  pl(6:7) = u2def(1:2);
  pl(8:10) = v2def;
  pl(11:12) = u3def(1:2);
  pl(13:15) = v3def;

  pl_full = zeros(18, 1);
  index_full = [1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18];
  pl_full(index_full) = pl;

  % calculating the stiffness matrix and internal force vector of the shell element in local coordinates
  [Kl_full, fintLocCoord] = localShellTriangle(x02, x03, y03, young_modulus, poisson_ratio, h, pl_full);

  % reducing to 15 dofs
  Kl = Kl_full(index_full, index_full);

  % local internal force vector
  fl = Kl * pl;

  % eq. (15) of 10.1016/j.cma.2006.10.006
  Ta1 = matrixTa(R1def);
  Ta2 = matrixTa(R2def);
  Ta3 = matrixTa(R3def);

  % eq. (18) of 10.1016/j.cma.2006.10.006
  fa = fl;
  fa(3:5) = Ta1' * fl(3:5);
  fa(8:10) = Ta2' * fl(8:10);
  fa(13:15) = Ta3' * fl(13:15);

  % eq. (19) of 10.1016/j.cma.2006.10.006
  Ba = eye(15);
  Ba(3:5, 3:5) = Ta1;
  Ba(8:10, 8:10) = Ta2;
  Ba(13:15, 13:15) = Ta3;

  % eq. (21) of 10.1016/j.cma.2006.10.006
  Kh = zeros(15, 15);
  Kh(3:5, 3:5) = matrixKhi(R1def, fl(3:5));
  Kh(8:10, 8:10) = matrixKhi(R2def, fl(8:10));
  Kh(13:15, 13:15) = matrixKhi(R3def, fl(13:15));

  % eq. (20) of 10.1016/j.cma.2006.10.006
  % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  Ka = Ba' * Kl * Ba + Kh;

  % eq. (7) of 10.1016/j.cma.2006.10.006
  a1 = u1def + r1o;
  a2 = u2def + r2o;
  a3 = u3def + r3o;

  % eq. (27) of 10.1016/j.cma.2006.10.006
  [G1, G2, G3] = matrixGi(a1, a2, a3, r1o, r2o, r3o);
  G = [G1; G2; G3];

  % eq. (26) of 10.1016/j.cma.2006.10.006
  P = matrixP(a1, a2, a3, G1, G2, G3);

  % eq. (25) of 10.1016/j.cma.2006.10.006
  E = blkdiag(Rr, Rr, Rr, Rr, Rr, Rr);

  % eq. (30) of 10.1016/j.cma.2006.10.006
  n = P' * fa; % eq. (31) of 10.1016/j.cma.2006.10.006
  [F1, F2] = matrixF(n);

  % eq. (29) of 10.1016/j.cma.2006.10.006
  % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  fg = E * n;
  Kg = E * (P' * Ka * P  - G * F1' * P - F2 * G') * E';

  % eq. (39) of 10.1016/j.cma.2006.10.006
  Bm = eye(18);
  Bm(4:6, 4:6) = matrixTm(q1);
  Bm(10:12, 10:12) = matrixTm(q2);
  Bm(16:18, 16:18) = matrixTm(q3);

  % eq. (40) of 10.1016/j.cma.2006.10.006
  Kk = zeros(18, 18);
  Kk(4:6, 4:6) = matrixKki(q1, fg(4:6));
  Kk(10:12, 10:12) = matrixKki(q2, fg(10:12));
  Kk(16:18, 16:18) = matrixKki(q3, fg(16:18));

  % eq.(38) of 10.1016/j.cma.2006.10.006
  % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  fm = Bm' * fg;
  Km = Bm' * Kg * Bm + Kk;

  % shifting lines and columns to onsas convention of dofs ordering
  ks = {switchToNodalIndexing(Km)};
  fs = {switchToNodalIndexing(fm)};

end

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

function [B] = cstB(x02, x03, y03)
  % calculate the strain-displacement matrix for the constant stress triangular element (CST)
  % x02, x03 and y03 are the local coordinates of the nodes 2 and 3

  area02 = x02 * y03;

  % ( with y01 = y02 = x01 = 0)
  % bi = (yj - yk )/2A
  b1 = -y03 / area02;
  b2 =   y03 / area02;
  b3 =  0;

  % ci = (xk - xj)/2A
  c1 = (x03 - x02) / area02;
  c2 = -x03 / area02;
  c3 = x02 / area02;

  B = [[b1,  0, b2,  0, b3,  0]
       [0, c1,  0, c2,  0, c3]
       [c1, b1, c2, b2, c3, b3]];

end

function [R] = globalRotationMatrix(q)
  % Eq. (35) of 10.1016/j.cma.2006.10.006
  q1 = q(1);
  q2 = q(2);
  q3 = q(3);
  q0 = sqrt(1.0 - q(1)^2 - q(2)^2 - q(3)^2);

  R  = 2 * [[(q0^2 + q1^2 - 0.5),  (q1 * q2 - q0 * q3),    (q1 * q3 + q0 * q2)]
            [(q1 * q2 + q0 * q3),  (q0^2 + q2^2 - 0.5),    (q2 * q3 - q0 * q1)]
            [(q1 * q3 - q0 * q2),  (q2 * q3 + q0 * q1),    (q0^2 + q3^2  - 0.5)]];
end

function [v] = rotationVector(R)
  % Eq. (13) of 10.1016/j.cma.2006.10.006
  v = zeros(3, 1);
  v(1) = .5 * (R(3, 2) - R(2, 3));
  v(2) = .5 * (R(1, 3) - R(3, 1));
  v(3) = .5 * (R(2, 1) - R(1, 2));
end

function [Ta] = matrixTa(R)
  % Eq. (15) of 10.1016/j.cma.2006.10.006
  Ta = 0.5 * [[R(2, 2) + R(3, 3), -R(1, 2),        -R(1, 3)]
              [-R(2, 1),           R(1, 1) + R(3, 3), -R(2, 3)]
              [-R(3, 1),          -R(3, 2),         R(1, 1) + R(2, 2)]
             ];
end

function [Tm] = matrixTm(q)
  % Eq. (37) of 10.1016/j.cma.2006.10.006
  q1 = q(1);
  q2 = q(2);
  q3 = q(3);
  q0 = sqrt(1.0 - q(1)^2 - q(2)^2 - q(3)^2);

  Tm = 2 / q0 * [[(q0^2 + q1^2),    (q1 * q2 - q0 * q3),  (q1 * q3 + q0 * q2)]
                 [(q1 * q2 + q0 * q3),  (q0^2 + q2^2),    (q2 * q3 - q0 * q1)]
                 [(q1 * q3 - q0 * q2),  (q2 * q3 + q0 * q1),  (q0^2 + q3^2)]];
end

function [Khi] = matrixKhi(R, m)
  % Eq. (21) of 10.1016/j.cma.2006.10.006

  Khi = 0.5 * [[((R(2, 3) - R(3, 2)) * m(1) + R(3, 1) * m(2) - R(2, 1) * m(3)),   (R(1, 1) * m(3) - R(1, 3) * m(1)),                            (R(1, 2) * m(1) - R(1, 1) * m(2))]
               [(R(2, 3) * m(2) - R(2, 2) * m(3)),                          ((R(3, 1) - R(1, 3)) * m(2) + R(1, 2) * m(3) - R(3, 2) * m(1)),     (R(2, 2) * m(1) - R(2, 1) * m(2))]
               [(R(3, 3) * m(2) - R(3, 2) * m(3)),                          (R(3, 1) * m(3) - R(3, 3) * m(1)),                            ((R(1, 2) - R(2, 1)) * m(3) + R(2, 3) * m(1) - R(1, 3) * m(2))]];

end

function [G1, G2, G3] = matrixGi(a1, a2, a3, r1, r2, r3)
  % Eq. (27) of 10.1016/j.cma.2006.10.006

  a21 = a2 - a1;
  a31 = a3 - a1;
  a32 = a3 - a2;
  a13 = a1 - a3;

  v = norm(cross(a21, a31));
  c = a1' * r1 + a2' * r2 + a3' * r3;

  G1 = zeros(6, 3);
  G1(1, 3) = -r1(2) / c;
  G1(2, 3) =   r1(1) / c;
  G1(3, 1) =   a32(1) / v;
  G1(3, 2) =   a32(2) / v;

  G2 = zeros(6, 3);
  G2(1, 3) = -r2(2) / c;
  G2(2, 3) =   r2(1) / c;
  G2(3, 1) =   a13(1) / v;
  G2(3, 2) =   a13(2) / v;

  G3 = zeros(6, 3);
  G3(1, 3) = -r3(2) / c;
  G3(2, 3) =   r3(1) / c;
  G3(3, 1) =   a21(1) / v;
  G3(3, 2) =   a21(2) / v;

end

function [P] = matrixP(a1, a2, a3, G1, G2, G3)
  % Eq. (26) of 10.1016/j.cma.2006.10.006

  Ai = zeros(5, 3);
  Ai(3, 1) = 1;
  Ai(4, 2) = 1;
  Ai(5, 3) = 1;

  I = zeros(5, 6);
  I(1, 1) = 1;
  I(2, 2) = 1;
  I(3, 4) = 1;
  I(4, 5) = 1;
  I(5, 6) = 1;

  P = zeros(15, 18);
  Ai(1, 3) = -a1(2);
  Ai(2, 3) =   a1(1);
  P1 = [I - Ai * G1', -Ai * G2', -Ai * G3'];

  Ai(1, 3) = -a2(2);
  Ai(2, 3) =   a2(1);
  P2 = [-Ai * G1', I - Ai * G2', -Ai * G3'];

  Ai(1, 3) = -a3(2);
  Ai(2, 3) =   a3(1);
  P3 = [-Ai * G1', -Ai * G2', I - Ai * G3'];

  P = [P1; P2; P3];
end

function [F1, F2] = matrixF(n)
  % Eq. (30) of 10.1016/j.cma.2006.10.006

  F1 = zeros(15, 3);
  F2 = zeros(18, 3);

  for i = 1:3
    j1 = 6 * (i - 1) + 1;
    j2 = j1 + 2;
    j3 = j2 + 1;
    j4 = 6 * i;

    sna = skew(n(j1:j2));
    snb = skew(n(j3:j4));

    i1 = 5*(i - 1) + 1;
    i2 = i1 + 1;
    F1(i1:i2, :) = sna(1:2, :);

    F2(j1:j2, :) = sna;
    F2(j3:j4, :) = snb;

  end
end

function [Kki] = matrixKki(q, m)
  % Eq. (41) of 10.1016/j.cma.2006.10.006

  q0s = 1.0 - q(1)^2 - q(2)^2 - q(3)^2;
  if q0s < 0
    q0s;
    q
    error('ojo');
  end
  q0 = sqrt(q0s);
  A = q(1) * m(1) + q(2) * m(2) + q(3) * m(3);

  H = zeros(3, 3);
  H(1, 1) = (q0s + q(1)^2) * A;
  H(2, 2) = (q0s + q(2)^2) * A;
  H(3, 3) = (q0s + q(3)^2) * A;
  H(1, 2) = q0s * (q(1) * m(2) - q(2) * m(1) - q0 * m(3)) + q(1) * q(2) * A;
  H(1, 3) = q0s * (q(1) * m(3) - q(3) * m(1) + q0 * m(2)) + q(1) * q(3) * A;
  H(2, 1) = q0s * (q(2) * m(1) - q(1) * m(2) + q0 * m(3)) + q(2) * q(1) * A;
  H(2, 3) = q0s * (q(2) * m(3) - q(3) * m(2) - q0 * m(1)) + q(2) * q(3) * A;
  H(3, 1) = q0s * (q(3) * m(1) - q(1) * m(3) - q0 * m(2)) + q(3) * q(1) * A;
  H(3, 2) = q0s * (q(3) * m(2) - q(2) * m(3) + q0 * m(1)) + q(3) * q(2) * A;

  Kki = (2 / q0^3) * H;

end
