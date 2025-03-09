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
% Integration Gauss Points based on https://keisan.casio.com/exec/system/1329114617
function [xIntPoints, wIntPoints] = gaussPointsAndWeights (numGaussPoints)
  if numGaussPoints == 1
    xIntPoints = 0;
    wIntPoints = 2;

  elseif numGaussPoints == 2

    xIntPoints = [-sqrt(1 / 3) sqrt(1 / 3)];

    wIntPoints = [1          1];

  elseif numGaussPoints == 3

    xIntPoints = [-sqrt(3 / 5)     0  sqrt(3 / 5)];

    wIntPoints = [5 / 9      8 / 9     5 / 9];

  elseif numGaussPoints == 4

    xIntPoints = [-sqrt(3 - 2 * sqrt(6 / 5)) / sqrt(7),  sqrt(3 - 2 * sqrt(6 / 5)) / sqrt(7) ...
                  -sqrt(3 + 2 * sqrt(6 / 5)) / sqrt(7),  sqrt(3 + 2 * sqrt(6 / 5)) / sqrt(7)];

    wIntPoints = [(18 + sqrt(30)) / 36                   (18 + sqrt(30)) / 36      ...
                  (18 - sqrt(30)) / 36                   (18 - sqrt(30)) / 36];

  elseif numGaussPoints == 5

    xIntPoints = [-0.9061798459386639927976, -0.5384693101056830910363, 0, ...
                  0.5384693101056830910363,  0.9061798459386639927976];

    wIntPoints = [0.2369268850561890875143,  0.4786286704993664680413, 0.5688888888888888888889, ...
                  0.4786286704993664680413,  0.2369268850561890875143];

  elseif numGaussPoints == 6

    xIntPoints = [-0.9324695142031520278123, -0.661209386466264513661, -0.2386191860831969086305, ...
                  0.238619186083196908631,  0.661209386466264513661, 0.9324695142031520278123];

    wIntPoints = [0.1713244923791703450403, 0.3607615730481386075698, 0.4679139345726910473899,  ...
                  0.46791393457269104739,   0.3607615730481386075698, 0.1713244923791703450403];

  elseif numGaussPoints == 7

    xIntPoints = [-0.9491079123427585245262, -0.7415311855993944398639, -0.4058451513773971669066, ...
                  0,  0.4058451513773971669066,  0.7415311855993944398639, ...
                  0.9491079123427585245262];

    wIntPoints = [0.1294849661688696932706,  0.2797053914892766679015, 0.38183005050511894495, ...
                  0.417959183673469387755,  0.38183005050511894495,  0.279705391489276667901, ...
                  0.129484966168869693271];

  elseif numGaussPoints == 8

    xIntPoints = [-0.9602898564975362316836, -0.7966664774136267395916, -0.5255324099163289858177, ...
                  -0.1834346424956498049395, 0.1834346424956498049395, 0.5255324099163289858177, ...
                  0.7966664774136267395916, 0.9602898564975362316836];

    wIntPoints = [0.1012285362903762591525, 0.2223810344533744705444, 0.313706645877887287338, ...
                  0.3626837833783619829652, 0.3626837833783619829652, 0.313706645877887287338, ...
                  0.222381034453374470544, 0.1012285362903762591525];

  elseif numGaussPoints == 9

    xIntPoints = [-0.9681602395076260898356, -0.8360311073266357942994, -0.6133714327005903973087, ...
                  -0.3242534234038089290385,  0,  0.3242534234038089290385, ...
                  0.6133714327005903973087,  0.8360311073266357942994,  0.9681602395076260898356];

    wIntPoints = [0.0812743883615744119719,  0.1806481606948574040585,  0.2606106964029354623187, ...
                  0.312347077040002840069,   0.330239355001259763165,  0.312347077040002840069, ...
                  0.260610696402935462319,   0.1806481606948574040585,  0.081274388361574411972];

  elseif numGaussPoints == 10

    xIntPoints = [-0.973906528517171720078,  -0.8650633666889845107321, -0.6794095682990244062343, ...
                  -0.4333953941292471907993, -0.1488743389816312108848,  0.1488743389816312108848, ...
                  0.4333953941292471907993,  0.6794095682990244062343,  0.8650633666889845107321, ...
                  0.973906528517171720078];

    wIntPoints = [0.0666713443086881375936,   0.149451349150580593146,   0.219086362515982043996, ...
                  0.2692667193099963550912,   0.2955242247147528701739,  0.295524224714752870174, ...
                  0.269266719309996355091,    0.2190863625159820439955,  0.1494513491505805931458, ...
                  0.0666713443086881375936];

  elseif numGaussPoints == 12

    xIntPoints = [-0.98156063424672E+00, ...
                  -0.90411725637047E+00, ...
                  -0.76990267419431E+00, ...
                  -0.58731795428662E+00, ...
                  -0.36783149899818E+00, ...
                  -0.12523340851147E+00, ...
                  0.12523340851147E+00, ...
                  0.36783149899818E+00, ...
                  0.58731795428662E+00, ...
                  0.76990267419431E+00, ...
                  0.90411725637047E+00, ...
                  0.98156063424672E+00];

    wIntPoints = [0.47175336386511E-01, ...
                  0.10693932599532E+00, ...
                  0.16007832854335E+00, ...
                  0.20316742672307E+00, ...
                  0.23349253653835E+00, ...
                  0.24914704581340E+00, ...
                  0.24914704581340E+00, ...
                  0.23349253653835E+00, ...
                  0.20316742672307E+00, ...
                  0.16007832854335E+00, ...
                  0.10693932599532E+00, ...
                  0.47175336386511E-01];

  elseif numGaussPoints == 14

    xIntPoints = [-9.86283808696812e-01, ...
                  -9.28434883663574e-01, ...
                  -8.27201315069765e-01, ...
                  -6.87292904811685e-01, ...
                  -5.15248636358154e-01, ...
                  -3.19112368927890e-01, ...
                  -1.08054948707344e-01, ...
                  1.08054948707344e-01, ...
                  3.19112368927890e-01, ...
                  5.15248636358154e-01, ...
                  6.87292904811685e-01, ...
                  8.27201315069765e-01, ...
                  9.28434883663574e-01, ...
                  9.86283808696812e-01];

    wIntPoints = [3.51194603317519e-02, ...
                  8.01580871597602e-02, ...
                  1.21518570687903e-01, ...
                  1.57203167158193e-01, ...
                  1.85538397477938e-01, ...
                  2.05198463721296e-01, ...
                  2.15263853463158e-01, ...
                  2.15263853463158e-01, ...
                  2.05198463721296e-01, ...
                  1.85538397477938e-01, ...
                  1.57203167158193e-01, ...
                  1.21518570687903e-01, ...
                  8.01580871597602e-02, ...
                  3.51194603317519e-02];

  elseif numGaussPoints == 16

    xIntPoints = [-9.89400934991650e-01, ...
                  -9.44575023073233e-01, ...
                  -8.65631202387832e-01, ...
                  -7.55404408355003e-01, ...
                  -6.17876244402644e-01, ...
                  -4.58016777657227e-01, ...
                  -2.81603550779259e-01, ...
                  -9.50125098376374e-02, ...
                  9.50125098376374e-02, ...
                  2.81603550779259e-01, ...
                  4.58016777657227e-01, ...
                  6.17876244402644e-01, ...
                  7.55404408355003e-01, ...
                  8.65631202387832e-01, ...
                  9.44575023073233e-01, ...
                  9.89400934991650e-01];

    wIntPoints = [2.71524594117541e-02, ...
                  6.22535239386479e-02, ...
                  9.51585116824928e-02, ...
                  1.24628971255534e-01, ...
                  1.49595988816577e-01, ...
                  1.69156519395003e-01, ...
                  1.82603415044924e-01, ...
                  1.89450610455069e-01, ...
                  1.89450610455069e-01, ...
                  1.82603415044924e-01, ...
                  1.69156519395003e-01, ...
                  1.49595988816577e-01, ...
                  1.24628971255534e-01, ...
                  9.51585116824928e-02, ...
                  6.22535239386479e-02, ...
                  2.71524594117541e-02];

  elseif numGaussPoints == 28

    xIntPoints = [-0.9964424975739544, ...
                  -0.9813031653708728, ...
                  -0.9542592806289382, ...
                  -0.9156330263921321, ...
                  -0.8658925225743951, ...
                  -0.8056413709171791, ...
                  -0.7356108780136318, ...
                  -0.6566510940388649, ...
                  -0.5697204718114017, ...
                  -0.4758742249551183, ...
                  -0.3762515160890787, ...
                  -0.2720616276351781, ...
                  -0.1645692821333808, ...
                  -0.05507928988403425, ...
                  0.05507928988403425, ...
                  0.1645692821333808, ...
                  0.2720616276351781, ...
                  0.3762515160890787, ...
                  0.4758742249551183, ...
                  0.5697204718114017, ...
                  0.6566510940388649, ...
                  0.7356108780136318, ...
                  0.8056413709171791, ...
                  0.8658925225743951, ...
                  0.9156330263921321, ...
                  0.9542592806289382, ...
                  0.9813031653708728, ...
                  0.9964424975739544];

    wIntPoints = [0.009124282593094496, ...
                  0.02113211259277136, ...
                  0.03290142778230417, ...
                  0.04427293475900421, ...
                  0.055107345675716686, ...
                  0.06527292396699941, ...
                  0.07464621423456871, ...
                  0.08311341722890121, ...
                  0.0905717443930328, ...
                  0.09693065799792995, ...
                  0.1021129675780607, ...
                  0.10605576592284632, ...
                  0.10871119225829412, ...
                  0.11004701301647522, ...
                  0.11004701301647522, ...
                  0.10871119225829412, ...
                  0.10605576592284632, ...
                  0.1021129675780607, ...
                  0.09693065799792995, ...
                  0.0905717443930328, ...
                  0.08311341722890121, ...
                  0.07464621423456871, ...
                  0.06527292396699941, ...
                  0.055107345675716686, ...
                  0.04427293475900421, ...
                  0.03290142778230417, ...
                  0.02113211259277136, ...
                  0.009124282593094496];

  else
    error('The number of gauss quadrature points introduced is not implemented yet. please open an issue!');
  end
end
