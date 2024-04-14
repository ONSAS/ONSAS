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
% nu  = modelParams(8) ;

uvector     = elemDisps([1,7]) ;    % x
vvector     = elemDisps([3,9]) ;    % y
thetavector = elemDisps([6,12]) ;   % theta z

Kfd    = zeros(6,6) ;
Kfalfa = zeros(6,6) ;
Khd    = zeros(6,6) ;
Khalfa = 0 ;

Cep    = 0 ;
Fint   = 0 ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;
xpi = [0 l/2 l] ;
wpi = [1/3 4/3 1/3]*l*0.5 ;

% integration (Gauss-Lobatto)
% and calculation of values of internal parameters at integration points

% initial values of bulk moments
M1 = zeros(npi,1) ;

% initial value at hinge --not present yet--
alfan1 = 0 ;

% params_plastic_2Dframe [kpn(1:3), xin1(4:6), xin2(7:9), soft_hinge_boolean(10), xd(11), alpha(12), tM(13)]
kpn  = params_plastic_2Dframe(1:3) ;
xin1 = params_plastic_2Dframe(4:6) ;
xin2 =  params_plastic_2Dframe(7) ;

soft_hinge_boolean = params_plastic_2Dframe(8) ;

xd      = params_plastic_2Dframe(9) ;
alfan   = params_plastic_2Dframe(10) ;
tM      = params_plastic_2Dframe(11) ;
xdi     = params_plastic_2Dframe(12) ;

% set initial values of the parameters for time n + 1
kpn1  = zeros(3,1) ;
xin11 = zeros(3,1) ;
xin21 = 0 ; 

for ii = 1:npi

    [soft_hinge_boolean, Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, M1xpi, xd, Fi] ...
       = integrand_plastic(soft_hinge_boolean, ii, xpi(ii), xd, l, A, ...
         uvector, vvector, thetavector, alfan, xin1, kpn, E, Iy, My, Mc, kh1, kh2, Cep) ;

    % stiffness matrices / integration (Gauss-Lobatto)
    Kfd    = Kfd    + Kfdj    * wpi(ii) ;
    Kfalfa = Kfalfa + Kfalfaj * wpi(ii) ;
    Khd    = Khd    + Khdj    * wpi(ii) ;
    Khalfa = Khalfa + Khalfaj * wpi(ii) ;
    
    % values of internal parameters at integration points
    kpn1(ii)  = kpn1xpi ;
    xin11(ii) = xin11xpi ;
    
    % internal forces / integration (Gauss-Lobatto)
    Fint = Fint + Fi*wpi(ii) ;

    M1(ii) = M1xpi ;

end

Khalfa = Khalfa + Ks ; % integral + Ks

% element stiffness matrix
if soft_hinge_boolean == true

    Kelement = Kfd - Kfalfa*Khalfa^(-1)*Khd ;

else

    Kelement = Kfd ;

end

for ii = 1:npi

if abs(M1(ii)) >= Mu && soft_hinge_boolean == false

    soft_hinge_boolean = true ;

    xd = xpi(ii) ;
    xdi = ii ;

end

end

if soft_hinge_boolean == true

    tM = 0 ;

    for ii = 1:npi

        Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(ii)/l)) ;

        % integration (Gauss-Lobatto)
        tM = tM - Ghatxpi*M1(ii)*wpi(ii) ;

    end

[soft_hinge_boolean, alfan1, xin21, xd] = soft_hinge(soft_hinge_boolean, xd, alfan, xin2, tM, l, E, Iy, Mu, Ks) ;

fprintf('\n | alpha = %8.4f | tM = %8.4f\n |' , alfan1, tM) ;

end

Fintout = zeros(12,1) ;
KTout = zeros(12,12) ;

dofsconv = [1 1+6 3 3+6 6 6+6] ;
Fintout(dofsconv) = Fint ;
KTout(dofsconv, dofsconv) = Kelement ;

if norm(elemDisps)>1e-8 && norm(Fint)<1e-8 && norm(KTout*elemDisps)>1e-8

Fintout = [Fintout KTout*elemDisps] ;

end

fs = {Fintout} ;
ks = {KTout} ;

params_plastic_2Dframe_np1 = zeros(1,12);

if soft_hinge_boolean == true

% once the hinge is formed, we assume that the plastic deformations in the bulk
% will not be changing any more

    params_plastic_2Dframe_np1(1:3) = kpn ;
    params_plastic_2Dframe_np1(4:6) = xin1 ;

else

    params_plastic_2Dframe_np1(1:3) = kpn1 ;
    params_plastic_2Dframe_np1(4:6) = xin11 ;

end

params_plastic_2Dframe_np1(7) = xin21 ;
params_plastic_2Dframe_np1(8) = soft_hinge_boolean ;
params_plastic_2Dframe_np1(9) = xd ;
params_plastic_2Dframe_np1(10) = alfan1 ;
params_plastic_2Dframe_np1(11) = tM ;
params_plastic_2Dframe_np1(12) = xdi ;

end