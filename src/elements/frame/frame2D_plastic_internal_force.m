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

function [ fs , ks, params_plastic_2Dframe_np1] = frame2D_plastic_internal_force( elemNodesxyzRefCoords , ...
    elemCrossSecParams    , ...
    modelParams , ...
    elemDisps , params_plastic_2Dframe )

% \/
% called by the function assembler
% /\

% initial/deformed lengths
Bdif = [ -eye(3) eye(3) ] ;
l = sqrt(sum((Bdif*elemNodesxyzRefCoords').^2)) ;
A  = elemCrossSecParams{2}(1) ;
Iy = elemCrossSecParams{2}(3) ;

% --- elastoplastic params ---
E   = modelParams(1) ;
Mc  = modelParams(2) ;
My  = modelParams(3) ;
Mu  = modelParams(4) ;
kh1 = modelParams(5) ;
kh2 = modelParams(6) ;
Ks  = modelParams(7) ;

uvector     = elemDisps([1,7]) ;    % x
vvector     = elemDisps([3,9])  ;   % y
thetavector = elemDisps([6,12]) ;   % theta z

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;
xpi = [0 l/2 l] ;
wpi = [1/3 4/3 1/3]*l*0.5 ;


% ==========================================================
% variables de estado candidatas
% ==========================================================

% renombro a vars locales
kp_n      = params_plastic_2Dframe(1:3) ;
xi1_n     = params_plastic_2Dframe(4:6) ;
xi2_n     = params_plastic_2Dframe(7) ;
SH_bool_n = params_plastic_2Dframe(8) ; % true if in the n time is active the softening state
xd_n      = params_plastic_2Dframe(9) ;  % hinge coordinate
alfa_n    = params_plastic_2Dframe(10) ; % alpha in time n
tM_n      = params_plastic_2Dframe(11) ; % hinge moment
xdi_n     = params_plastic_2Dframe(12) ; % number of the integration point where is the hinge

% candidates for state var for time n+1
kp_np1      = kp_n ;
xi1_np1     = xi1_n ;
xi2_np1     = xi2_n ;
SH_bool_np1 = SH_bool_n ;
xd_np1      = xd_n ;
alfa_np1    = alfa_n ; % alpha in time n
tM_np1      = tM_n ; % hinge moment
xdi_np1     = xdi_n ; % number of the integration point where is the hinge


% ==========================================================
% calculo momentos
% ==========================================================

% integration (Gauss-Lobatto)
% and calculation of values of internal parameters at integration points

% initial values of bulk moments
[ Mnp1, tM_np1, Ghats] = frame_plastic_IPmoments( E, Iy, vvector, thetavector, xpi, xd_np1, l, alfa_np1, kp_np1, wpi)

max_abs_mom = max( abs ( Mnp1)) ;

% ==========================================================
% solve local equations
% ==========================================================

if ( SH_bool_n == false && max_abs_mom > Mu ) || SH_bool_np1 == true
  % solve softening step
  [alfa_np1, xi2_np1, xd_np1] = softening_step(xd_n, alfa_n, xi2_n, tM_np1, l, E, Iy, Mu, Ks) ;

  if SH_bool_n == false, SH_bool_np1 = true; end 

else
  % solve plastic bending step
   [ kp_np1, xi1_np1, Cep_np1] = plastic_hardening_step( E, Iy, vvector, thetavector, xpi, xi1_n, kp_n, My, Mc, kh1, kh2, Mnp1)
  
end

[ Mnp1, tM_np1, Ghats] = frame_plastic_IPmoments( E, Iy, vvector, thetavector, xpi, xd_np1, l, alfa_np1, kp_np1, wpi)


% ==========================================================
% solve global equations
% ==========================================================

[ Kfd, Kfalfa, Khd, Khalfa] = frame_plastic_matrices( )

if soft_hinge_boolean == true || soft_hinge_boolean_np1 == true
    Kelement = Kfd - Kfalfa*Khalfa^(-1)*Khd ;
else
    Kelement = Kfd ;
end






% \/ Softening hinge activated

for ii = 1:npi
    if abs(M1(ii)) >= Mu && soft_hinge_boolean_np1 == false
        soft_hinge_boolean_np1 = true ;
        xd = xpi(ii) ;
        xdi = ii ;
    end
end

if soft_hinge_boolean_np1 == true || soft_hinge_boolean == true

    [alfan1, xin21, xd] = soft_hinge(xd, alfan, xin2, tM, l, E, Iy, Mu, Ks) ;

    fprintf('\n | alpha = %8.8f | tM = %8.4f\n |', alfan1, tM) ;
    fprintf('\n | displacement y = %8.4f\n |', vvector(2)) ;

end

Fintout = zeros(12,1) ;
KTout = zeros(12,12) ;

dofsconv = [1 1+6 3 3+6 6 6+6] ;
Fintout(dofsconv) = Fint ;
KTout(dofsconv, dofsconv) = Kelement ;

fs = {Fintout} ;
ks = {KTout} ;

params_plastic_2Dframe_np1 = zeros(1,12) ;

params_plastic_2Dframe_np1(1:3) = kpn1 ;
params_plastic_2Dframe_np1(4:6) = xin11 ;
params_plastic_2Dframe_np1(7)   = xin21 ;
params_plastic_2Dframe_np1(8)   = soft_hinge_boolean_np1 ;
params_plastic_2Dframe_np1(9)   = xd ;
params_plastic_2Dframe_np1(10)  = alfan1 ;
params_plastic_2Dframe_np1(11)  = tM ;
params_plastic_2Dframe_np1(12)  = xdi ;

end