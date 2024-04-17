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
[ Mnp1, tM_np1] = frame_plastic_IPmoments( E, Iy, vvector, thetavector, xpi, xd_np1, l, alfa_np1, kp_np1, wpi)

max_abs_mom = max( abs ( Mnp1)) ;

% ==========================================================
% solve local equations
% ==========================================================

if ( SH_bool_n == false && max_abs_mom > Mu ) || SH_bool_np1 == true
  % solve softening step
  [alfa_np1, xi2_np1, xd_np1] = softening_step(xd_n, alfa_n, xi2_n, tM_np1, l, E, Iy, Mu, Ks) ;
else
  % solve plastic bending step
  [alfa_np1, xi2_np1, xd_np1] = hardening_step(xd_n, alfa_n, xi2_n, tM_np1, l, E, Iy, Mu, Ks) ;

end

# [soft_hinge_boolean, Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, M1xpi, xd, Fi] ...
# = integrand_plastic(soft_hinge_boolean, ii, xpi(ii), xd, l, A, ...
#   uvector, vvector, thetavector, alfan, xin1, kpn, E, Iy, My, Mc, kh1, kh2, Cep) ;


% params_plastic_2Dframe [kpn(1:3), xin1(4:6), xin2(7), soft_hinge_boolean(8), xd(9), alpha(10), tM(11), xdi(12)]



Cep    = 0 ;
Fint   = 0 ;





% ==========================================================
% calculo matrices
% ==========================================================

Kfd    = zeros(6, 6) ;
Kfalfa = zeros(6, 1) ;
Khd    = zeros(1, 6) ;
Khalfa = 0 ;



for ii = 1:npi

    [soft_hinge_boolean, Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, M1xpi, xd, Fi] ...
       = integrand_plastic(soft_hinge_boolean, ii, xpi(ii), xd, l, A, ...
         uvector, vvector, thetavector, alfan, xin1, kpn, E, Iy, My, Mc, kh1, kh2, Cep) ;

    % stiffness matrices / integration (Gauss-Lobatto)
    Kfd    = Kfd    + Kfdj    * wpi(ii) ;
    Kfalfa = Kfalfa + Kfalfaj * wpi(ii) ;
    Khd    = Khd    + Khdj    * wpi(ii) ;
    Khalfa = Khalfa + Khalfaj * wpi(ii) ;
    
    size(Kfd)
    size(Kfalfa) 
    size(Khd) 
    size(Khalfa)

    size(Kfdj)
    size(Kfalfaj) 
    size(Khdj) 
    size(Khalfaj)

    stop
    
    % values of internal parameters at integration points
    kpn1(ii)  = kpn1xpi ;
    xin11(ii) = xin11xpi ;
    
    % internal forces / integration (Gauss-Lobatto)
    Fint = Fint + Fi*wpi(ii) ;

    % moments at the integration points for candidate displacements (correspondings of time n + 1)
    M1(ii) = M1xpi ;

if soft_hinge_boolean_np1 == true

    tM = 0 ;

    for jj = 1:npi

        Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(jj)/l)) ;

        % integration (Gauss-Lobatto)
        % tM calculated with the moments M1 corresponding to time n + 1
        tM = tM - Ghatxpi*M1(jj)*wpi(jj) ;
    
    end

else
    
alfan1 = alfan ;

end

end

% integral + Ks
Khalfa = Khalfa + Ks ;

% element stiffness matrix
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