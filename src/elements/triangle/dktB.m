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

function B = dktB(PSI, ETA, x02, x03, y03)
  % calculate the strain-displacement matrix for the triangular plate element (DKT)
  % psi and eta are the area coordinates of the triangular element
  % x02, x03 and y03 are the local coordinates of the nodes 2 and 3

  area02 = x02 * y03;

  % auxiliar variables - geometric
  X32 = x03 - x02;
  XX2 = x02 * x02;
  XX3 = x03 * x03;
  XX32 = X32 * X32;
  YY3 = y03 * y03;
  QQ4 = 1.0 / (XX32 + YY3);
  QQ5 = 1.0 / (XX3 + YY3);
  QQ6 = 1.0 / XX2;
  SU4 = 3.0 * QQ4 * YY3;
  SU5 = 3.0 * QQ5 * YY3;
  D4  = 3.0 * QQ4 * X32 * y03;
  D5  = 3.0 * QQ5 * x03 * y03;
  CL4 = 6.0 * QQ4 * X32;
  CL5 = -x03 * 6.0 * QQ5;
  CL6 = 6.0 * QQ6 * x02;
  SL4 = 6.0 * QQ4 * y03;
  SL5 = -y03 * 6.0 * QQ5;

  % auxiliar variables - coordinates
  PS2  = 1.0 - 2.0 * PSI;
  ET2  = 1.0 - 2.0 * ETA;
  SU5E = SU5 * ET2;

  % Derivada de betax com relacao a psi
  BXP1 = CL6 * PS2 + ETA * (CL5 - CL6);
  BXP2 = -ETA * D5;
  BXP3 = -4.0 + 6.0 * (PSI + ETA) - ETA * SU5;
  BXP4 = -CL6 * PS2 + ETA * (CL4 + CL6);
  BXP5 = ETA * D4;
  BXP6 = -2.0 + 6.0 * PSI + ETA * SU4;
  BXP7 = -ETA * (CL5 + CL4);
  BXP8 = ETA * (D4 - D5);
  BXP9 = -ETA * (SU5 - SU4);

  % Derivada de betay com relacao a psi
  BYP1 = ETA * SL5;
  BYP2 = 1.0 - ETA * SU5;
  BYP3 = -BXP2;
  BYP4 = ETA * SL4;
  BYP5 = -1.0 + ETA * SU4;
  BYP6 = -BXP5;
  BYP7 = -ETA * (SL5 + SL4);
  BYP8 = BXP9;
  BYP9 = -BXP8;

  % Derivada de betax com relacao a eta
  BXE1 = -CL5 * ET2 - PSI * (CL6 - CL5);
  BXE2 = D5 * (ET2 - PSI);
  BXE3 = -4.0 + 6.0 * (PSI + ETA) + SU5E - PSI * SU5;
  BXE4 = PSI * (CL4 + CL6);
  BXE5 = PSI * D4;
  BXE6 = PSI * SU4;
  BXE7 = CL5 * ET2 - PSI * (CL4 + CL5);
  BXE8 = D5 * ET2 + PSI * (D4 - D5);
  BXE9 = -2.0 + 6.0 * ETA + SU5E + PSI * (SU4 - SU5);

  % Derivada de betay com relacao a eta
  BYE1 = SL5 * (PSI - ET2);
  BYE2 = 1.0 + SU5E - PSI * SU5;
  BYE3 = -BXE2;
  BYE4 = PSI * SL4;
  BYE5 = BXE6;
  BYE6 = -BXE5;
  BYE7 = SL5 * ET2 - PSI * (SL4 + SL5);
  BYE8 = -1.0 + SU5E + PSI * (SU4 - SU5);
  BYE9 = -BXE8;

  % Determinacao de [bb]
  B(1, 1) = y03 * BXP1 / area02;
  B(2, 1) = (-x03 * BYP1 + x02 * BYE1) / area02;
  B(3, 1) = (-x03 * BXP1 + x02 * BXE1 + y03 * BYP1) / area02;
  B(1, 2) = y03 * BXP2 / area02;
  B(2, 2) = (-x03 * BYP2 + x02 * BYE2) / area02;
  B(3, 2) = (-x03 * BXP2 + x02 * BXE2 + y03 * BYP2) / area02;
  B(1, 3) = y03 * BXP3 / area02;
  B(2, 3) = (-x03 * BYP3 + x02 * BYE3) / area02;
  B(3, 3) = (-x03 * BXP3 + x02 * BXE3 + y03 * BYP3) / area02;
  B(1, 4) = y03 * BXP4 / area02;
  B(2, 4) = (-x03 * BYP4 + x02 * BYE4) / area02;
  B(3, 4) = (-x03 * BXP4 + x02 * BXE4 + y03 * BYP4) / area02;
  B(1, 5) = y03 * BXP5 / area02;
  B(2, 5) = (-x03 * BYP5 + x02 * BYE5) / area02;
  B(3, 5) = (-x03 * BXP5 + x02 * BXE5 + y03 * BYP5) / area02;
  B(1, 6) = y03 * BXP6 / area02;
  B(2, 6) = (-x03 * BYP6 + x02 * BYE6) / area02;
  B(3, 6) = (-x03 * BXP6 + x02 * BXE6 + y03 * BYP6) / area02;
  B(1, 7) = y03 * BXP7 / area02;
  B(2, 7) = (-x03 * BYP7 + x02 * BYE7) / area02;
  B(3, 7) = (-x03 * BXP7 + x02 * BXE7 + y03 * BYP7) / area02;
  B(1, 8) = y03 * BXP8 / area02;
  B(2, 8) = (-x03 * BYP8 + x02 * BYE8) / area02;
  B(3, 8) = (-x03 * BXP8 + x02 * BXE8 + y03 * BYP8) / area02;
  B(1, 9) = y03 * BXP9 / area02;
  B(2, 9) = (-x03 * BYP9 + x02 * BYE9) / area02;
  B(3, 9) = (-x03 * BXP9 + x02 * BXE9 + y03 * BYP9) / area02;

end
