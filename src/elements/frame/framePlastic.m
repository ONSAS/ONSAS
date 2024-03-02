% Copyright 2023, ONSAS Authors (see documentation)
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

function [ soft_hinge_boolean, Fint, M1, Kelement, kpn1, xin11, xin21, alfan1, xd, tM] = framePlastic(soft_hinge_boolean, dnk, kpn, xin1, xin2, alfan, xd, tM, elemParams, elastoplasticParams)


% --- element params ---
l  = elemParams(1) ;
A  = elemParams(2) ;
Iy = elemParams(3) ;

% --- elastoplastic params ---
E   = elastoplasticParams(1) ;
Mc  = elastoplasticParams(2) ;
My  = elastoplasticParams(3) ;
Mu  = elastoplasticParams(4) ;
kh1 = elastoplasticParams(5) ;
kh2 = elastoplasticParams(6) ;
Ks  = elastoplasticParams(7) ;

uvector     = dnk(1:2) ;
vvector     = dnk(3:4) ;
thetavector = dnk(5:6) ;

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

% set initial values of internal parameters at integration points
kpn1  = kpn  ;
xin11 = xin1 ;
xin21 = xin2 ;

for ii = 1:npi

    [soft_hinge_boolean, Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, xin21xpi, M1xpi, xd, Fi, alfan1] = integrand_plastic(soft_hinge_boolean, ii, xpi(ii), xd, l, A, uvector, vvector, thetavector, alfan, xin1, xin2, kpn, E, Iy, My, Mc, Mu, kh1, kh2, Ks, Cep, tM) ;
    
    % stiffness matrices / integration (Gauss-Lobatto)
    Kfd    = Kfd    + Kfdj    * wpi(ii) ;
    Kfalfa = Kfalfa + Kfalfaj * wpi(ii) ;
    Khd    = Khd    + Khdj    * wpi(ii) ;
    Khalfa = Khalfa + Khalfaj * wpi(ii) ;
    
    % values of internal parameters at integration points
    kpn1(ii)  = kpn1xpi  ;
    xin11(ii) = xin11xpi ;
    xin21(ii) = xin21xpi ;
    
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

if M1(ii) >= Mu && soft_hinge_boolean == false

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

end

end