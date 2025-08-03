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
% Implementation of a triangular finite element with 6 dfos (3 translations and 3 rotations) per node for the analysis of linear elastic isotropic shells with constant thickness.
% The element is formed by the superposition of a plate element (DKT) and a plane stress element (CST) with
% with addition to artificial drilling (rotation about the axis normal to the element plane) stiffness.
%
function [fs, ks, fintLocCoord, rotMat] = internalForcesShellTriangle(elemCoords, elemDisps, modelName, modelParams, thickness, rotMat)

  origin = 0;
  e1_parallel_side12 = 1;  % 0 if battini modification is used - 1 if not
  flag_first_mod  = 0; % 1 if battini modification is used - 0 if not / local rotations
  flag_second_mod = 0;  % 1 if battini modification is used - 0 if not / out of plane disps = 0
  flag_third_mod  = 0;  % 1 if battini modification is used - 0 if not / quaternions
  flag_OPT        = 1;  % 1 if battini modification is used - 0 if not / quaternions

  % material and geometric parameters
  young_modulus = modelParams(1);
  poisson_ratio = modelParams(2);
  h = thickness;

  % Nodes position vector in global reference frame
  elemCoords = elemCoords';
  r1_g = elemCoords(1:3);
  r2_g = elemCoords(4:6);
  r3_g = elemCoords(7:9);
  if origin == 0
    rc_g = (r1_g + r2_g + r3_g) / 3;
  elseif origin == 1
    rc_g = r1_g;
  end

  % Disps and spatial rotations in global reference frame
  Ug = switchToTypeIndexing(elemDisps);
  u1_g    = Ug(1:3);
  u2_g    = Ug(7:9);
  u3_g    = Ug(13:15);
  t1_g    = Ug(4:6);
  t2_g    = Ug(10:12);
  t3_g    = Ug(16:18);

  % Updated position vector in global reference frame
  p1_g = r1_g + u1_g;
  p2_g = r2_g + u2_g;
  p3_g = r3_g + u3_g;
  if origin == 0
    pog = (p1_g + p2_g + p3_g) / 3;
  elseif origin == 1
    pog = p1_g;
  end
  % fprintf('p_g \n')
  % [ p1_g p2_g p3_g pog ]

  % Global rotation matrix
  % eq. (35) of 10.1016/j.cma.2006.10.006
  % R1_g = rotMat{1};
  % R2_g = rotMat{2};
  % R3_g = rotMat{3};
  R1_g = expm(skew(t1_g));
  R2_g = expm(skew(t2_g));
  R3_g = expm(skew(t3_g));
  % R1_g = expon(t1_g);
  % R2_g = expon(t2_g);
  % R3_g = expon(t3_g);
  % stop

  % Transformation matrices from global reference frame
  [To, x02, x03, y03] = edgeLocalAxisShellTriangle(r1_g, r2_g, r3_g);
  [Tr, ~, ~, ~]       = edgeLocalAxisShellTriangle(p1_g, p2_g, p3_g);

  % Rotation matrix from global reference frame to local reference frame in initial configuration
  Ro = To';
  % Rotation matrix from global reference frame to local reference frame in deformed configuration
  Rr = Tr';

  % nodal displacements in local reference frame in deformed configuration
  % eq. (1) of 10.1016/j.cma.2006.10.006
  r1_o = To * (r1_g - rc_g);
  r2_o = To * (r2_g - rc_g);
  r3_o = To * (r3_g - rc_g);
  % fprintf('ri_o \n')
  % [ r1_o r2_o r3_o ]

  u1_def = Rr' * (p1_g - pog) - r1_o;
  u2_def = Rr' * (p2_g - pog) - r2_o;
  u3_def = Rr' * (p3_g - pog) - r3_o;
  % fprintf('u_def \n')
  % [ u1_def u2_def u3_def ]
  % stop

  % eq. (7) of 10.1016/j.cma.2006.10.006
  a1_def = u1_def + r1_o;
  a2_def = u2_def + r2_o;
  a3_def = u3_def + r3_o;
  % [a1_def a2_def a3_def]

  if e1_parallel_side12 == 0
    Num = 0;
    Den = 0;
    a   = [a1_def, a2_def, a3_def];
    ro  = [r1_o, r2_o, r3_o];
    for i = 1:3
      ai = a(:, i);
      rio = ro(:, i);
      auxNum = ai(2) * rio(1) - ai(1) * rio(2);
      Num = Num + auxNum;

      auxDen = ai(1) * rio(1) + ai(2) * rio(2);
      Den = Den + auxDen;
    end
    tan_theta = Num / Den;
    theta = rad2deg(atan(tan_theta));
    % num = a1_def(2)*r1_o(1) - a1_def(1)*r1_o(2) + a2_def(2)*r2_o(1) - a2_def(1)*r2_o(2) + a3_def(2)*r3_o(1) - a3_def(1)*r3_o(2)
    % den = a1_def(1)*r1_o(1) + a1_def(2)*r1_o(2) + a2_def(1)*r2_o(1) + a2_def(2)*r2_o(2) + a3_def(1)*r3_o(1) + a3_def(2)*r3_o(2)
    % stop
    % Falta rotar aun y escribir coordenadas de nodos actualizadas
  end

  % eq. (27) of 10.1016/j.cma.2006.10.006
  [G1, G2, G3] = matrixGi(a1_def, a2_def, a3_def, r1_o, r2_o, r3_o, e1_parallel_side12, origin);
  G = [G1; G2; G3];
  % [ sum(G(:,1)) sum(G(:,2)) sum(G(:,3)) ]
  % stop

  % Rotation matrix from local reference frame to nodal reference frame in deformed configuration
  % eq. (2) of 10.1016/j.cma.2006.10.006
  R1_def = Rr' * R1_g * Ro;
  R2_def = Rr' * R2_g * Ro;
  R3_def = Rr' * R3_g * Ro;
  % stop

  % Nodal rotations in local reference frame in deformed configuration
  % eq. (13) of 10.1016/j.cma.2006.10.006
  v1_def = rotationVector(R1_def, flag_first_mod);
  v2_def = rotationVector(R2_def, flag_first_mod);
  v3_def = rotationVector(R3_def, flag_first_mod);
  % [ v1_def v2_def v3_def ] ;
  % stop
  % Local displacement vector in local reference frame in deformed configuration
  % eq. (12) of 10.1016/j.cma.2006.10.006
  pl_full = zeros(18, 1);
  uz_dofs = [3, 9, 15];
  index_full = (1:18);
  % im = [1, 2, 7, 8, 13, 14];              % Membrane dofs (u, v)
  % ib = [3, 4, 5, 9, 10, 11, 15, 16, 17];  % bending dofs (w, rx, ry)
  % drill_dofs = [6, 12, 18];              % (rz)

  if flag_second_mod == 1
    pl = zeros(15, 1);
    pl(1:2)   = u1_def(1:2);
    pl(3:5)   = v1_def;
    pl(6:7)   = u2_def(1:2);
    pl(8:10)  = v2_def;
    pl(11:12) = u3_def(1:2);
    pl(13:15) = v3_def;
    index_full(uz_dofs) = [];
  else
    pl = zeros(18, 1);
    pl(1:3)   = u1_def;
    pl(4:6)   = v1_def;
    pl(7:9)   = u2_def;
    pl(10:12) = v2_def;
    pl(13:15) = u3_def;
    pl(16:18) = v3_def;
  end
  pl_full(index_full) = pl;

  % calculating the linear stiffness matrix and internal force vector of the shell element in local coordinates
  [Kl_full, fintLocCoord, Kb] = localShellTriangle(x02, x03, y03, young_modulus, poisson_ratio, h, pl_full, flag_OPT, r1_g, r2_g, r3_g, Tr);

  % Reduces Kl matrix to number of dofs considered
  Kl = Kl_full(index_full, index_full);
  % Kl =(Kl + Kl')/2;
  % local internal force vector
  fl = Kl * pl;

  % eq. (19) of 10.1016/j.cma.2006.10.006
  if flag_first_mod == 1 && flag_second_mod == 1
    % eq. (15) of 10.1016/j.cma.2006.10.006
    Ta1 = matrixTa(R1_def);
    Ta2 = matrixTa(R2_def);
    Ta3 = matrixTa(R3_def);
    Ba = eye(15);
    Ba(3:5, 3:5)      = Ta1;
    Ba(8:10, 8:10)    = Ta2;
    Ba(13:15, 13:15)  = Ta3;
  else
    Ba = eye(18);
    invTs_1 = invTsShell(v1_def);
    invTs_2 = invTsShell(v2_def);
    invTs_3 = invTsShell(v3_def);
    Ba(4:6, 4:6)      = invTs_1;
    Ba(10:12, 10:12)  = invTs_2;
    Ba(16:18, 16:18)  = invTs_3;
  end

  % eq. (18) of 10.1016/j.cma.2006.10.006
  fa = Ba' * fl;

  % eq. (21) of 10.1016/j.cma.2006.10.006
  if flag_first_mod == 1 && flag_second_mod == 1
    Kh = zeros(15, 15);
    Kh(3:5, 3:5)      = matrixKhi(R1_def, fl(3:5));
    Kh(8:10, 8:10)    = matrixKhi(R2_def, fl(8:10));
    Kh(13:15, 13:15)  = matrixKhi(R3_def, fl(13:15));
  else
    Kh = zeros(18, 18);
    Kh1 = dinvTsShell(v1_def, fl(4:6))   * invTs_1;
    Kh2 = dinvTsShell(v2_def, fl(10:12)) * invTs_2;
    Kh3 = dinvTsShell(v3_def, fl(16:18)) * invTs_3;
    Kh(4:6, 4:6)      = Kh1;
    Kh(10:12, 10:12)  = Kh2;
    Kh(16:18, 16:18)  = Kh3;
  end

  % eq. (20) of 10.1016/j.cma.2006.10.006
  Ka = Ba' * Kl * Ba + Kh;

  % eq. (26) of 10.1016/j.cma.2006.10.006
  [P, A] = matrixP(a1_def, a2_def, a3_def, G1, G2, G3, flag_second_mod);

  % eq. (25) of 10.1016/j.cma.2006.10.006
  E = blkdiag(Rr, Rr, Rr, Rr, Rr, Rr);

  % eq. (30) of 10.1016/j.cma.2006.10.006
  n = P' * fa; % eq. (31) of 10.1016/j.cma.2006.10.006

  [F1, F2] = matrixF(n, flag_second_mod);
  % [F1, F2] = matrixF(fa, flag_second_mod);
  % stop

  % eq. (29) of 10.1016/j.cma.2006.10.006
  fg = E * n;
  % [ switchToNodalIndexing(fg) fg n fa fl]
  % stop

  Kl = (P' * Ka * P - G * F1' * P - F2 * G');
  Kg = E * Kl * E';
  % Kg = (Kg+Kg')/2;
  % Rankin
  % F = 1/2 * (F1+F2);
  % Kg = E * (P' * Ka * P - G * F' * P - F * G') * E';

  % change of variables to rotation vector
  Ts1   = funTsShell(t1_g);
  Ts2   = funTsShell(t2_g);
  Ts3   = funTsShell(t3_g);
  Br = blkdiag(eye(3), Ts1, eye(3), Ts2, eye(3), Ts3);

  Kv1 = dTsShell(t1_g, fg(4:6));
  Kv2 = dTsShell(t2_g, fg(10:12));
  Kv3 = dTsShell(t3_g, fg(16:18));
  Kv = blkdiag(zeros(3, 3), Kv1, zeros(3, 3), Kv2, zeros(3, 3), Kv3);

  fr = Br' * fg;
  Kr = Br' * Kg * Br + Kv;

  % Kr= (Kr+Kr')/2;

  % fr = fg;
  % Kr = Kg;
  % Bm = eye(18);
  % Kk = zeros(18, 18);
  % if flag_third_mod == 1
  %   % eq. (39) of 10.1016/j.cma.2006.10.006
  %   Bm(4:6, 4:6)      = matrixTm(rot1_g);
  %   Bm(10:12, 10:12)  = matrixTm(rot2_g);
  %   Bm(16:18, 16:18)  = matrixTm(rot3_g);
  %   % eq. (40) of 10.1016/j.cma.2006.10.006
  %   Kk(4:6, 4:6)      = matrixKki(rot1_g, fg(4:6));
  %   Kk(10:12, 10:12)  = matrixKki(rot2_g, fg(10:12));
  %   Kk(16:18, 16:18)  = matrixKki(rot3_g, fg(16:18));
  % end

  % % eq.(38) of 10.1016/j.cma.2006.10.006
  % % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  % fm = Bm' * fg;
  % Km = Bm' * Kg * Bm + Kk;

  % shifting lines and columns to onsas convention of dofs ordering
  % ks = {switchToNodalIndexing(Km)};
  % fs = {switchToNodalIndexing(fm)};

  ks = {switchToNodalIndexing(Kr)};
  fs = {switchToNodalIndexing(fr)};

end
