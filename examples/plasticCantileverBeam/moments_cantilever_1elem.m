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

% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan

% =========================================================================

% numerical example
% cantilever beam loaded with a vertical force at the free end

% =========================================================================

% displacements in time n + 1, dpn1 v1, v2, theta1, theta2, alpha, xd

% plastic curvature in time n / kappa_plas_n

function Ms = moments_cantilever_1elem( v1, v2, theta1, theta2 , xd, alpha, xin1, kappa_plas_n, Mc, My, kh1, kh2, l)

x = [0;l/2;l];

Bv1 = -6/l^2*(1-2*x/l) ;
Bv2 =  6/l^2*(1-2*x/l) ;

Bt1 = -2/l*(2-3*x/l) ;
Bt2 = -2/l*(1-3*x/l) ;

G_bar = -(1+3*(1-2*xd/l)*(1-2*x/l))/l ; 

% CREO QUE NO while % no convergio

kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alpha ;

kappa_plas_n1 = kappa_plas_n ;

Mnp1_test = E*I*(kappa_bar - kappa_plas_n1) ;

if xin1 <= (My - Mc)/kh1

    q = -kh1*xin1 ;

else

    q = -(My - Mc)*(1-kh2/kh1) - kh2*xin1 ;

end

phi_test = Mnp1_test- (Mc - q) ;

if phi_test <0

end

% end

disp(Mnp1_test) ;
Ms = zeros(3,1) ;