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

function [ Fint, Kelement, kpn1, xin11, xin21, alfan1, xd, tM] = framePlastic( dnk, kpn, xin1, xin2, alfan, xd, elemParams, elastoplasticParams )


% --- element params ---
l = elemParams(1) ;
A = elemParams(2) ;
Iy = elemParams(3) ;

% --- elastoplastic params ---
E = elastoplasticParams(1) ;
Mc = elastoplasticParams(2) ;
My = elastoplasticParams(3) ;
Mu = elastoplasticParams(4) ;
kh1 = elastoplasticParams(5) ;
kh2 = elastoplasticParams(6) ;
Ks = elastoplasticParams(7) ;

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

uvector     = dnk(1:2) ;
vvector     = dnk(3:4) ;
thetavector = dnk(5:6) ;

Kfd    = zeros(6,6) ;
Kfalfa = zeros(6,6) ;
Khd    = zeros(6,6) ;
Khalfa = 0 ;
Cep    = 0 ;
Fint   = 0 ;

% trial value of the moment at the discontinuity (tM at xd)
tM = 0 ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;
xpi = [0 l/2 l] ;
wpi = [1/3 4/3 1/3] * l * 0.5 ;

kpn1  = zeros(npi,1) ;
xin11 = zeros(npi,1) ;
xin21 = zeros(npi,1) ;

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
% integration (Gauss-Lobatto)
% and calculation of values of internal parameters at integration points

for ii = 1:npi

    [Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, xin21xpi, ~, tM, xd, Fi, alfan1] = integrand_plastic(ii, xpi(ii), xd, l, uvector, vvector, thetavector, alfan, xin1, kpn, kpn1, E, Iy, My, Mc, kh1, kh2, A, Ks, xin2, Mu, Cep,tM) ;
    
    % stiffness matrices / integration (Gauss-Lobatto)
    Kfd    = Kfd    + Kfdj    * wpi(ii) ;
    Kfalfa = Kfalfa + Kfalfaj * wpi(ii) ;
    Khd    = Khd    + Khdj    * wpi(ii) ;
    Khalfa = Khalfa + Khalfaj * wpi(ii) ;
    
    % values of internal parameters at integration points
    kpn1(ii)  = kpn1xpi ;
    xin11(ii) = xin11xpi ;
    xin21(ii) = xin21xpi ;
    
    % internal forces / integration (Gauss-Lobatto)
    Fint = Fint + Fi*wpi(ii) ;
end

% Khalfa = Khalfa + Ks ; % integral + Ks

% element stiffness matrix

Kelement = Kfd ;%- Kfalfa*Khalfa^(-1)*Khd ;

for ii = 1:npi

    Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(ii)/l)) ;
    
    N = bendingInterFuns (xpi(ii), l, 2) ;
    Bv = [N(1) N(3)] ;
    Btheta = [N(2) N(4)] ;
    
    % curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)
    
    khat = Bv*vvector + Btheta*thetavector + Ghatxpi*alfan ;
    
    % khat2 = dirac(xd)*alfan ;
    % kn = khat + khat2 ;
    kenxpi = khat - kpn(ii) ;
    
    % moment at integration points
    Mixpi = E*Iy*kenxpi ;
    
    % integration (Gauss-Lobatto)
    
    tM = tM - Ghatxpi*Mixpi*wpi(ii) ;

end

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\