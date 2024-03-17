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

% algorithm for the update of the internal variables for elastoplasticity with hardening
% and for the computation of the moment in the bulk Mn1

% =========================================================================

% displacements in time n + 1, dpn1 v1, v2, theta1, theta2, alpha, xd
% plastic curvature in time n / kappa_plas_n

function [kappa_plas_n1, xin11, Mn1] = moments_plus_internal_variables( v1, v2, theta1, theta2 , xd, alpha, xin1, kappa_plas_n, Mc, My, kh1, kh2, E, Iy, l)

x = 0 ;
alpha = 0 ;

Bv1 = -6/l^2*(1-2*x/l) ;
Bv2 =  6/l^2*(1-2*x/l) ;

Bt1 = -2/l*(2-3*x/l) ;
Bt2 = -2/l*(1-3*x/l) ;

G_bar = -(1+3*(1-2*xd/l)*(1-2*x/l))/l ;

kappa_bar = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + G_bar*alpha ;

kappa_plas_test = kappa_plas_n ;

Mn1_test = E*Iy*(kappa_bar - kappa_plas_test) ;

if xin1 <= (My - Mc)/kh1

    q = -kh1*xin1 ;

else

    q = -(My - Mc)*(1-kh2/kh1) - kh2*xin1 ;

end

phi_test = abs(Mn1_test)- (Mc - q) ;

if phi_test <= 0

    kappa_plas_n1 = kappa_plas_n ;
    xin11 = xin1 ;
    Mn1 = Mn1_test ;

else

    % calculation of gamma_n1

    if xin1 + phi_test/(kh1 + E*Iy) <= (My - Mc)/kh1

        gamma_n1 = phi_test/(kh1 + E*Iy) ;

    else

        gamma_n1 = phi_test/(kh2 + E*Iy) ;

    end

    kappa_plas_n1 = kappa_plas_n + gamma_n1*sign(Mn1_test) ;
    xin11 = xin1 + gamma_n1 ;

    %{
    if gamma_n1 == 0

        Cep = E*Iy ;

    elseif gamma_n1 > 0 && xin11 <= (My - Mc)/kh1

        Cep = E*Iy*kh1/(E*Iy + kh1) ;

    elseif gamma_n1 > 0 && xin11 > (My - Mc)/kh1

        Cep = E*Iy*kh2/(E*Iy + kh2) ;

    end
    %}

    Cep = E*Iy ;

    Mn1 = Cep*(kappa_bar - kappa_plas_n1) ;

end