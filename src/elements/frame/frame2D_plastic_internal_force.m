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

function [ fs , ks, params_plastic_2Dframe] = frame2D_plastic_internal_force( elemNodesxyzRefCoords , ...
    elemCrossSecParams    , ...
    modelParams , ...
    elemDisps , params_plastic_2Dframe )
    
%    (soft_hinge_boolean, dnk, kpn, xin1, xin2, alfan, xd, tM, elemParams, elastoplasticParams)


% initial/deformed lengths
Bdif = [ -eye(3) eye(3) ] ;
l = sqrt( sum( ( Bdif * elemNodesxyzRefCoords'    ).^2 ) ) ;
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
nu  = modelParams(8) ;

uvector     = elemDisps([1,7]) ;  % x
vvector     = elemDisps([3,9]) ;  % y
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

%
kpn = params_plastic_2Dframe(1:3);
xin1 = params_plastic_2Dframe(4:6);
xin2 =  params_plastic_2Dframe(7:9);
soft_hinge_boolean = params_plastic_2Dframe(10)
xd = params_plastic_2Dframe(11)
alfan = params_plastic_2Dframe(12)


% set initial values of internal parameters at integration points
kpn1  = kpn  ;
xin11 = xin1 ;
xin21 = xin2 ;

for ii = 1:npi

    [soft_hinge_boolean, Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, xin21xpi, M1xpi, xd, Fi, alfan1] = integrand_plastic(soft_hinge_boolean, ii, xpi(ii), xd, l, A, uvector, vvector, thetavector, alfan, xin1, xin2, kpn, E, Iy, My, Mc, Mu, kh1, kh2, Ks, Cep, 0) ;
    
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

    M1(ii) = M1xpi ;

end

Khalfa = Khalfa + Ks ; % integral + Ks

% element stiffness matrix
if soft_hinge_boolean == 1
disp('hola')
    Kelement = Kfd - Kfalfa*Khalfa^(-1)*Khd ;

else
    
    Kelement = Kfd ;

end

tM = 0 ;

for ii = 1:npi

    Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(ii)/l)) ;

    % integration (Gauss-Lobatto)
    tM = tM - Ghatxpi*M1(ii)*wpi(ii) ;

end

if tM >= Mu && soft_hinge_boolean == false

    soft_hinge_boolean = true ;

    xd = 0 ;

end


Fintout = zeros(12,1) ;
KTout = zeros(12,12) ;

dofsconv = [1 1+6 3 3+6 5 5+6] ;
Fintout(dofsconv) = Fint ;
KTout(dofsconv,dofsconv) = Kelement ;

fs = {Fintout} ;
ks = {KTout} ;