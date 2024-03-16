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

% deplazamientos en tiempo n + 1, dpn1 v1, v2, theta1, theta2, alpha, xd

% curvatura plástica en tiempo n / kappas plas n

function Ms = moments_cantilever_1elem( v1, v2, theta1, theta2 , xd, alpha, l, kappas_plas_n )

x = [0;l/2;l];

Bv1 = -6/l^2*(1-2*x/l) ;
Bv2 =  6/l^2*(1-2*x/l) ;

Bt1 = -2/l*(2-3*x/l) ;
Bt2 = -2/l*(1-3*x/l) ;

% Gtecho = -1 ;  xd 

% CREO QUE NO while % no convergio

kappas_techo = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2 + (xd*alpha); % corregir + Gtecho * alpha

kappas_plas_n1 = kappas_plas_n ; 

Mnp1_test = E*I*( kappas_techo - kappas_plas_n1 ) ;

% phi_test = 

if phi_test <0

end

% end

disp(Mnp1_test) ;
Ms = zeros(3,1) ;