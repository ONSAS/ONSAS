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
function [fs, ks, fintLocCoord,Kl_full] = internalForcesShellTriangle(elemCoords, elemDisps, modelName, modelParams, thickness)

  % material and geometric parameters
  young_modulus = modelParams(1);
  poisson_ratio = modelParams(2);
  h = thickness;

  % elem nodal coords in global reference frame
  elemCoords = elemCoords';

  % ==============================================================================

  % Nodes position vector in global reference frame 
  r1g = elemCoords(1:3);
  r2g = elemCoords(4:6);
  r3g = elemCoords(7:9);
  rcg = (r1g + r2g + r3g) / 3;

  % Nodal disps in global reference frame
  Ug = switchToTypeIndexing(elemDisps);

  % Global disps and rotations
  u1g = Ug(1:3);
  q1  = Ug(4:6);
  u2g = Ug(7:9);
  q2  = Ug(10:12);
  u3g = Ug(13:15);
  q3  = Ug(16:18);

  % Updated position vector in global reference frame
  p1g = r1g + u1g;
  p2g = r2g + u2g;
  p3g = r3g + u3g;
  pog = (p1g + p2g + p3g) / 3;
  % p1g
  % p2g
  % p3g
  % stop
  % ==============================================================================
  
  % Global rotation matrix
  % eq. (35) of 10.1016/j.cma.2006.10.006
  R1g = globalRotationMatrix(q1);
  R2g = globalRotationMatrix(q2);
  R3g = globalRotationMatrix(q3);

  % Transformation matrices from global reference frame
  [To, x02, x03, y03] = edgeLocalAxisShellTriangle(r1g, r2g, r3g);
  [Tr, ~, ~, ~]       = edgeLocalAxisShellTriangle(p1g, p2g, p3g);

  

  % Rotation matrix from global reference frame to local reference frame in initial configuration
  Ro = To';
  % Rotation matrix from global reference frame to local reference frame in deformed configuration
  Rr = Tr';
  % Ro
  % Rr

  % nodal displacements in local reference frame in deformed configuration
  % eq. (1) of 10.1016/j.cma.2006.10.006
  r1o = To * (r1g - rcg);
  r2o = To * (r2g - rcg);
  r3o = To * (r3g - rcg);

  u1def = Rr' * (p1g - pog) - r1o
  u2def = Rr' * (p2g - pog) - r2o
  u3def = Rr' * (p3g - pog) - r3o

  % eq. (7) of 10.1016/j.cma.2006.10.006
  a1 = u1def + r1o;
  a2 = u2def + r2o;
  a3 = u3def + r3o;

  (u3def(2)-u2def(2))
  % stop

  % % Haugen
  %   u = u1g
  %   c = pog-rcg
  %   c_aux = ucg

  %   I = eye(3) ;
  %   Ro = Rr ;
  %   xo = r1g ;
  %   a = rcg ; 

  %   ud = u-c+(I-Ro)*(xo-a)
  %   ud_def = Rr'*ud
  %   stop

  e1_parallel_side12 = 1 ;
  if e1_parallel_side12 == 0
  
    % ri0 = Xi
    % ai = xi
    Num = 0 ;
    Den = 0 ;
    a = [a1,a2,a3] ;
    ro = [r1o,r2o,r3o] ;
    for i = 1:3
      ai = a(:,i) ;
      rio = ro(:,i) ;
      auxNum = ai(2)*rio(1) - ai(1)*rio(2) ;
      Num = Num + auxNum ;

      auxDen = ai(1)*rio(1) + ai(2)*rio(2) ;
      Den = Den + auxDen ;
    end  
    tan_theta = Num/Den 
    theta=rad2deg(atan(tan_theta))
    % num = a1(2)*r1o(1) - a1(1)*r1o(2) + a2(2)*r2o(1) - a2(1)*r2o(2) + a3(2)*r3o(1) - a3(1)*r3o(2)
    % den = a1(1)*r1o(1) + a1(2)*r1o(2) + a2(1)*r2o(1) + a2(2)*r2o(2) + a3(1)*r3o(1) + a3(2)*r3o(2)
    % stop
    % Falta rotar aun y escribir coordenadas de nodos actualizadas
  end

  % eq. (27) of 10.1016/j.cma.2006.10.006
  [G1, G2, G3] = matrixGi(a1, a2, a3, r1o, r2o, r3o, e1_parallel_side12);
  G = [G1; G2; G3];
  % G
  % sum(G(:,1))
  % sum(G(:,2))
  % sum(G(:,3))
  % stop

  % ==============================================================================

  % Rotation matrix from local reference frame to nodal reference frame in deformed configuration 
  % eq. (2) of 10.1016/j.cma.2006.10.006
  R1def = Rr' * R1g * Ro;
  R2def = Rr' * R2g * Ro;
  R3def = Rr' * R3g * Ro;

  % R1g
  % R2g
  % R3g
  % stop

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Nodal rotations in local reference frame in deformed configuration 
  % eq. (13) of 10.1016/j.cma.2006.10.006
  v1def = rotationVector(R1def);
  v2def = rotationVector(R2def);
  v3def = rotationVector(R3def);

  v1def 
  v2def
  v3def
  % stop
  % ==============================================================================

  % Local displacement vector in local reference frame in deformed configuration
  % eq. (12) of 10.1016/j.cma.2006.10.006
  pl = zeros(15, 1);
  pl(1:2) = u1def(1:2);
  pl(3:5) = v1def;
  pl(6:7) = u2def(1:2);
  pl(8:10) = v2def;
  pl(11:12) = u3def(1:2);
  pl(13:15) = v3def;
  
  pl_full = zeros(18, 1) ;
  uz_dofs = [ 3, 9, 15 ] ;
  % index_full = [1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18];
  index_full = (1:18);
  index_full(uz_dofs) = [];

  pl_full(index_full) = pl;

  % ==============================================================================

  % calculating the stiffness matrix and internal force vector of the shell element in local coordinates
  [Kl_full, fintLocCoord] = localShellTriangle(x02, x03, y03, young_modulus, poisson_ratio, h, pl_full);

  % reducing to 15 dofs
  Kl = Kl_full(index_full, index_full);

  % local internal force vector
  fl = Kl * pl;
  % pl
  % fl
  % stop
  % ==============================================================================

  % eq. (15) of 10.1016/j.cma.2006.10.006
  Ta1 = matrixTa(R1def);
  Ta2 = matrixTa(R2def);
  Ta3 = matrixTa(R3def);
  % Ta1
  % Ta2
  % Ta3
  % stop
  
  % eq. (19) of 10.1016/j.cma.2006.10.006
  Ba = eye(15);
  Ba(3:5, 3:5) = Ta1;
  Ba(8:10, 8:10) = Ta2;
  Ba(13:15, 13:15) = Ta3;
  % Ba

  % eq. (18) of 10.1016/j.cma.2006.10.006
  fa = Ba'*fl;
  % fa = fl;
  % fa(3:5) = Ta1' * fl(3:5);
  % fa(8:10) = Ta2' * fl(8:10);
  % fa(13:15) = Ta3' * fl(13:15);
  
  % ==============================================================================

  % eq. (21) of 10.1016/j.cma.2006.10.006
  Kh = zeros(15, 15);
  Kh(3:5, 3:5) = matrixKhi(R1def, fl(3:5));
  Kh(8:10, 8:10) = matrixKhi(R2def, fl(8:10));
  Kh(13:15, 13:15) = matrixKhi(R3def, fl(13:15));
  % Kh
  % stop

  % eq. (20) of 10.1016/j.cma.2006.10.006
  % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  Ka = Ba' * Kl * Ba + Kh;
  % dif_K = Kl ./ Ka
  % stop

  % eq. (26) of 10.1016/j.cma.2006.10.006
  P = matrixP(a1, a2, a3, G1, G2, G3);
  % P_f = matrixP_full(a1, a2, a3, G1, G2, G3);
  % P

  % eq. (25) of 10.1016/j.cma.2006.10.006
  E = blkdiag(Rr, Rr, Rr, Rr, Rr, Rr);

  % eq. (30) of 10.1016/j.cma.2006.10.006
  n = P' * fa; % eq. (31) of 10.1016/j.cma.2006.10.006
  % fa_aux = zeros(18,1) ;
  % fa_aux(index_full) = fa ;
  % [fa_aux n]
  % stop

  [F1, F2] = matrixF(n);
  
  % eq. (29) of 10.1016/j.cma.2006.10.006
  % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
  fg = E * n;
  %
  Kl_aux = (P' * Ka * P  - G * F1' * P - F2 * G') ;
  %
  Kg = E * (P' * Ka * P  - G * F1' * P - F2 * G') * E';
  
  % Kl_full ./ Kl_aux
  % stop

  % im = [1, 2, 7, 8, 13, 14];              % Membrane dofs (u, v)
  % ib = [3, 4, 5, 9, 10, 11, 15, 16, 17];  % bending dofs (w, rx, ry)
  % drill_dofs = [6, 12, 18] ;              % (rz)
  
  % Membrane matrix
  % ======================================================
  % K_linear_m  = Kl_full(im, im) ;
  % K_NL_m      = Kl_aux(im, im) ;
  % dif_K_m     = K_linear_m ./ K_NL_m
  % ======================================================
  % (K_linear_m - K_NL_m)
  % stop
  % bending matrix
  % ======================================================
  % K_linear_b = Kl_full(ib,ib) ;
  % K_NL_b = Kl_aux(ib, ib) ;
  % dif_K_b = K_linear_b ./ K_NL_b
  % ======================================================
  % (K_linear_b - K_NL_b)
  % T_lin = blkdiag(To,To,To,To,To,To) ;
  % fL = T_lin'* Kl_full * T_lin * Ug
  % fNL = Kg * Ug
  % fL./fNL
  % stop

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
  % ===========================================================================================================================
  Bm = eye(18) ;
  Kk = zeros(18,18) ;
  % ===========================================================================================================================
  fm = Bm' * fg;
  Km = Bm' * Kg * Bm + Kk;

  % shifting lines and columns to onsas convention of dofs ordering
  ks = {switchToNodalIndexing(Km)};
  fs = {switchToNodalIndexing(fm)};

  % ks = {Km};
  % fs = {fm};

  % fl
  % fa
  % fg
  % fm
  
end

function [R] = globalRotationMatrix(q) % ok
  % Eq. (35) of 10.1016/j.cma.2006.10.006
  q1 = q(1);
  q2 = q(2);
  q3 = q(3);
  q0 = sqrt(1.0 - q1^2 - q2^2 - q3^2);

  R  = 2 * [[(q0^2 + q1^2 - 0.5),  (q1 * q2 - q0 * q3),    (q1 * q3 + q0 * q2)]
            [(q1 * q2 + q0 * q3),  (q0^2 + q2^2 - 0.5),    (q2 * q3 - q0 * q1)]
            [(q1 * q3 - q0 * q2),  (q2 * q3 + q0 * q1),    (q0^2 + q3^2  - 0.5)]];
end

function [v] = rotationVector(R) % ok
  % Eq. (13) of 10.1016/j.cma.2006.10.006
  v = zeros(3, 1);
  v(1) = .5 * (R(3, 2) - R(2, 3));
  v(2) = .5 * (R(1, 3) - R(3, 1));
  v(3) = .5 * (R(2, 1) - R(1, 2));
end

function [Ta] = matrixTa(R) % ok
  % Eq. (15) of 10.1016/j.cma.2006.10.006
  Ta = 0.5 * [[ R(2, 2) + R(3, 3) , -R(1, 2)          , -R(1, 3)          ]
              [-R(2, 1)           ,  R(1, 1) + R(3, 3), -R(2, 3)          ]
              [-R(3, 1)           , -R(3, 2)          ,  R(1, 1) + R(2, 2)]
             ];
end

function [Tmi] = matrixTm(qi) % ok
  % Eq. (37) of 10.1016/j.cma.2006.10.006
  q1 = qi(1);
  q2 = qi(2);
  q3 = qi(3);
  q0 = sqrt(1.0 - q1^2 - q2^2 - q3^2);

  Tmi = 2 / q0 * [[(q0^2 + q1^2)      ,  (q1 * q2 - q0 * q3)  ,  (q1 * q3 + q0 * q2)]
                 [(q1 * q2 + q0 * q3) ,  (q0^2 + q2^2)        ,  (q2 * q3 - q0 * q1)]
                 [(q1 * q3 - q0 * q2) ,  (q2 * q3 + q0 * q1)  ,  (q0^2 + q3^2)      ]
                 ];
end

function [Khi] = matrixKhi(R, m) % ok
  % Eq. (21) of 10.1016/j.cma.2006.10.006

  Khi = 0.5 * [[((R(2, 3) - R(3, 2)) * m(1) + R(3, 1) * m(2) - R(2, 1) * m(3)),   (R(1, 1) * m(3) - R(1, 3) * m(1)),                            (R(1, 2) * m(1) - R(1, 1) * m(2))]
               [(R(2, 3) * m(2) - R(2, 2) * m(3)),                          ((R(3, 1) - R(1, 3)) * m(2) + R(1, 2) * m(3) - R(3, 2) * m(1)),     (R(2, 2) * m(1) - R(2, 1) * m(2))]
               [(R(3, 3) * m(2) - R(3, 2) * m(3)),                          (R(3, 1) * m(3) - R(3, 3) * m(1)),                            ((R(1, 2) - R(2, 1)) * m(3) + R(2, 3) * m(1) - R(1, 3) * m(2))]];

end

function [G1, G2, G3] = matrixGi(a1, a2, a3, r1, r2, r3, e1_flag) % ok
  % Eq. (27) of 10.1016/j.cma.2006.10.006

  a21 = a2 - a1;
  a31 = a3 - a1;
  a32 = a3 - a2;
  a13 = a1 - a3;

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
  e1_flag
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
    G3(1, 3) = 0;
    G3(2, 3) = 0;    
  end

end

function [P] = matrixP(a1, a2, a3, G1, G2, G3) % ok
  % Eq. (26) of 10.1016/j.cma.2006.10.006

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
  A1 = Ai ;
  A1(1, 3) = -a1(2);
  A1(2, 3) =  a1(1);
  P1 = [I - A1 * G1', -A1 * G2', -A1 * G3'];

  A2=Ai ;
  A2(1, 3) = -a2(2);
  A2(2, 3) =   a2(1);
  P2 = [-A2 * G1', I - A2 * G2', -A2 * G3'];

  A3 = Ai ;
  A3(1, 3) = -a3(2);
  A3(2, 3) =   a3(1);
  P3 = [-A3 * G1', -A3 * G2', I - A3 * G3'];

  P = [P1; P2; P3];
end

function [P] = matrixP_full(a1, a2, a3, G1, G2, G3) % ok
  % Eq. (26) of 10.1016/j.cma.2006.10.006

  % Matrix A
  Ai = zeros(6, 3);
  Ai(3, 1) = 1;
  Ai(4, 2) = 1;
  Ai(5, 3) = 1;
  % Identity matrix
  I = eye(6) ;

  % Projector Matrix
  P = zeros(18, 18);
  % Node 1
  A1 = Ai ;
  A1(1, 3) = -a1(2);
  A1(2, 3) =  a1(1);
  %
  A1(3, 1) =  a1(2);
  A1(3, 2) = -a1(1);
  %
  P1 = [I - A1 * G1', -A1 * G2', -A1 * G3'];

  % Node 2
  A2 = Ai;
  A2(1, 3) = -a2(2);
  A2(2, 3) =   a2(1);
  %
  A2(3, 1) =  a2(2);
  A2(3, 2) = -a2(1);
  %
  P2 = [-A2 * G1', I - A2 * G2', -A2 * G3'];

  % node 3
  A3 = Ai ;
  A3(1, 3) = -a3(2);
  A3(2, 3) =   a3(1);
  %
  A3(3, 1) =  a3(2);
  A3(3, 2) = -a3(1);
  %
  P3 = [-A3 * G1', -A3 * G2', I - A3 * G3'];

  P = [P1; P2; P3];
end

function [F1, F2] = matrixF(n) % ok
  % Eq. (30) of 10.1016/j.cma.2006.10.006

  % F1 = zeros(15, 3);
  % F2 = zeros(18, 3);

  % for i = 1:3
  %   j1 = 6 * (i - 1) + 1;
  %   j2 = j1 + 2;
  %   j3 = j2 + 1;
  %   j4 = 6 * i;

  %   sna = skew(n(j1:j2));
  %   snb = skew(n(j3:j4));

  %   i1 = 5*(i - 1) + 1;
  %   i2 = i1 + 1;
  %   F1(i1:i2, :) = sna(1:2, :);

  %   F2(j1:j2, :) = sna;
  %   F2(j3:j4, :) = snb;

  % end
  
  n1 = n(1:3) ;
  n2 = n(4:6) ;

  n3 = n(7:9) ;
  n4 = n(10:12) ;

  n5 = n(13:15) ;
  n6 = n(16:18) ; 

  F1 = [ skew(n1)(1:2,:)' zeros(3,3) skew(n3)(1:2,:)' zeros(3,3) skew(n5)(1:2,:)' zeros(3,3) ]' ;
  F2 = [ skew(n1)' skew(n2)' skew(n3)' skew(n4)' skew(n5)' skew(n6)' ]' ;

end

function [Kki] = matrixKki(q, m) % ok
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
