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
nu  = modelParams(8) ;

uvector     = elemDisps([1,7]) ;    % x
vvector     = elemDisps([3,9]) ;    % y
thetavector = elemDisps([6,12]) ;   % theta z

Kfd    = zeros(6,6) ;
Kfalfa = zeros(6,6) ;
Khd    = zeros(6,6) ;
Khalfa = 0 ;

Cep    = 0 ;
Fint   = 0 ;
alfan1 = 0 ;

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
xin2 =  params_plastic_2Dframe(7:9) ;

soft_hinge_boolean = params_plastic_2Dframe(10) ;

xd      = params_plastic_2Dframe(11) ;
alfan   = params_plastic_2Dframe(12) ;
tM      = params_plastic_2Dframe(13) ;
xdi     = params_plastic_2Dframe(14) ;

if soft_hinge_boolean == true

    tM = 0 ;

    for ii = 1:npi

        Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(ii)/l)) ;

        % integration (Gauss-Lobatto)
        tM = tM - Ghatxpi*M1(ii)*wpi(ii) ;

    end

end

% set candidate values of internal parameters for next time at integration points
kpn1  = kpn  ;
xin11 = xin1 ;
xin21 = xin2 ; 

for ii = 1:npi

        if soft_hinge_boolean == true

            thetavector(2) = thetavector(1) + alfan ;
            vvector(2) = vvector(1) + xd*thetavector(1) + (l-xd)*(alfan + thetavector(1)) ;

        end

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

    thetavector(2) = thetavector(1) + alfan ;
    vvector(2) = vvector(1) + xd*thetavector(1) + (l-xd)*(alfan + thetavector(1)) ;

end

% hinge_softening_module 

% plastic softening at the discontinuity
% the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity
% softening criterion (failure function) at integration points

if soft_hinge_boolean == true

qfailxpi = min(-Ks*xin2(xdi), Mu) ;

phifailxpi = abs(tM)-(Mu-qfailxpi) ;

if phifailxpi <= 0
   
    alfan1 = alfan ;
    xin21 = xin2 ;

else

    if  xin2(xdi)<=-Mu/Ks

        gamma2 = phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks) ;

    else

        gamma2 = abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)) ;
    
    end
    
    alfan1      = alfan     + gamma2*sign(tM) ;
    xin21       = xin2      + gamma2 ;

end

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

params_plastic_2Dframe_np1 = zeros(1,13);

if soft_hinge_boolean == true

    params_plastic_2Dframe_np1(1:3) = kpn ;
    params_plastic_2Dframe_np1(4:6) = xin1 ;

else

    params_plastic_2Dframe_np1(1:3) = kpn1 ;
    params_plastic_2Dframe_np1(4:6) = xin11 ;

end

params_plastic_2Dframe_np1(7:9) = xin21 ;
params_plastic_2Dframe_np1(10) = soft_hinge_boolean ;
params_plastic_2Dframe_np1(11) = xd ;
params_plastic_2Dframe_np1(12) = alfan1 ;
params_plastic_2Dframe_np1(13) = tM ;
params_plastic_2Dframe_np1(14) = xdi ;

params_plastic_2Dframe(1:3) = kpn  ;
params_plastic_2Dframe(4:6) = xin1;
params_plastic_2Dframe(7:9) = xin2;

params_plastic_2Dframe(10) = soft_hinge_boolean ;

params_plastic_2Dframe(11) = xd ;
params_plastic_2Dframe(12) = alfan ;
params_plastic_2Dframe(13) = tM ;
params_plastic_2Dframe(14) = xdi ;

end