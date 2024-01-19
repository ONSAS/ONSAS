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

% =========================================================================

function [ dnk1, kpn1, xin11, xin21, alfan1, xd, Fint, tM] = framePlastic( dnk, kpn, xin1, xin2, alfan, xd, Fint, tM, elemParams, elastoplasticParams, lambda)

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

  % uvector = [u1; u2] ;
  % vvector = [v1; v2] ;
  % thetavector = [theta1; theta2] ;

  uvector = dnk(1:2) ;
  vvector = dnk(3:4) ;
  thetavector = dnk(5:6) ;
  
  Bu = [-1/l 1/l] ;

  Kfd = zeros(6,6) ;
  Kfalfa = zeros(6,6) ;
  Khd = zeros(6,6) ;
  Khalfa = 0 ;

  % Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]

  npi = 3 ;
  xpi = [0 l/2 l] ;
  wpi = [1/3 4/3 1/3] ;

  kpn1 = zeros(npi,1) ;
  xin11 = zeros(npi,1) ;
  xin21 = zeros(npi,1) ;

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
  % integration (Gauss-Lobatto)
  % and calculation of values of internal parameters at integration points

  for ii = 1:npi

    [Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, xin21xpi, ~, tM, xd, Fi] = integrand(ii, xpi(ii), xd) ;


    % stiffness matrices / integration (Gauss-Lobatto)

    Kfd = Kfd + Kfdj*wpi(ii) ;
    Kfalfa = Kfalfa + Kfalfaj*wpi(ii) ;
    Khd = Khd + Khdj*wpi(ii) ;
    Khalfa = Khalfa + Khalfaj*wpi(ii) ;

    % values of internal parameters at integration points

    kpn1(ii) = kpn1xpi ;
    xin11(ii) = xin11xpi ;
    xin21(ii) = xin21xpi ;

    % internal forces / integration (Gauss-Lobatto)

    Fint = Fint + Fi*wpi(ii) ;

  end

  Khalfa = Khalfa + Ks ; % integral + Ks

  % element stiffness matrix

  Kelement = Kfd - Kfalfa*Khalfa^(-1)*Khd ;

  Kelement(1,:) = [] ;
  Kelement(2,:) = [] ;
  Kelement(3,:) = [] ;

  Kelement(:,1) = [] ;
  Kelement(:,2) = [] ;
  Kelement(:,3) = [] ;

  Fint(1) = [] ;
  Fint(2) = [] ;
  Fint(3) = [] ;
 
  % system of equilibrium equations

  deltad = Kelement\([0 lambda 0]' - Fint) ;

  deltad = [0; 0; 0; deltad] ;

  Fint = [0; 0; 0; Fint] ;

  dnk1 = dnk + deltad ;

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

    function [Kfd, Kfalfa, Khd, Khalfa, kpn1xpi, xin11xpi, xin21xpi, M1xpi, tM, xd, Fi] = integrand(jj, xpi, xd)

    % elastoplasticity with hardening
    % the usual trial-corrector (return mapping) algorithm
    % is used at each of the Gauss–Lobatto integration points.

    N = bendingInterFuns (xpi, l, 2) ;
    Bv = [N(1) N(3)] ;
    Btheta = [N(2) N(4)] ;

    Bd = [Bu 0 0 0 0; 0 0 Bv Btheta] ;

    Ghat = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi/l)) ;

    % curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)

    khat = Bv*vvector + Btheta*thetavector + Ghat*alfan ;
    % khat2 = dirac(xd)*alfan ;
    % kn = khat + khat2 ;
    ken = khat - kpn(jj) ;

    % moment

    Mxpi = E*Iy*ken ;

    % yield criterion

    if xin1(jj) <= (My-Mc)/kh1

        qxpi = -kh1*xin1(jj) ;

    else

        qxpi = -(My-Mc)*(1-kh2/kh1)-kh2*xin1(jj) ;

    end
    
    % qxpi = piecewise(xin1(jj) <= (My-Mc)/kh1, -kh1*xin1(jj), -(My-Mc)*(1-kh2/kh1)-kh2*xin1(jj)) ;

    phixpi = abs(Mxpi) - (Mc - qxpi) ;

    % test values

    % kpn1test = kpn ;
    % xin11test = xin1 ;
    phitest =  phixpi ;

    % gamma values calculations (gamma derivative is the plastic multiplier)
    % the new values of internal variables are computed

    if phitest <= 0
  
        gamma = 0 ;
        kpn1xpi = kpn(jj) ;
        xin11xpi = xin1(jj) ;
        M1xpi = Mxpi ;

    else

        if xin1(jj) + phitest/(kh1+E*Iy)<=(My-Mc)/kh1

            gamma = phitest/(kh1+E*Iy) ;

        else

            gamma = phitest/(kh2+E*Iy) ;
        
        end

        % gamma = piecewise(xin1(jj) + phitest/(kh1+E*Iy)<=(My-Mc)/kh1, phitest/(kh1+E*Iy), phitest/(kh2+E*Iy)) ;
        
        kpn1xpi = kpn(jj) + gamma*sign(Mxpi) ;
        xin11xpi = xin1(jj) + gamma ;
        M1xpi = E*Iy*(khat-kpn1(jj)) ;

    end

    % elastoplastic tangent bending modulus

    if gamma == 0

        Cep = E*Iy ;

    elseif gamma > 0 && xin11xpi <= (My-Mc)/kh1

        Cep = E*Iy*kh1/(E*Iy + kh1) ;

    elseif gamma > 0 && xin11xpi > (My-Mc)/kh1

        Cep = E*Iy*kh2/(E*Iy + kh2) ;

    end

    % Cep = piecewise(gamma == 0, E*Iy , gamma > 0 & xin11xpi <= (My-Mc)/kh1, E*Iy*kh1/(E*Iy + kh1), gamma > 0 & xin11xpi > (My-Mc)/kh1, E*Iy*kh2/(E*Iy + kh2)) ;
  
    % stiffness matrices

    Kfd = Bd'*[E*A 0; 0 Cep]*Bd ;

    Kfalfa = Bd'*[E*A 0; 0 Cep]*[0 Ghat]' ;

    Khd = [0 Ghat]*[E*A 0; 0 Cep]*Bd ;

    Khalfa = Ghat*Cep*Ghat ;

    epsilon = Bu*uvector ;

    Fi= Bd'*[E*A*epsilon; M1xpi] ;

    % plastic softening at the discontinuity
    % the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity

    % test values
  
    % alfan1test = alfan ;
    % xin2test = xin2 ;

    % trial value of the moment at the discontinuity (tM at xd)
    tM = 0 ;

    % softening criterion (failure function) at integration points

    qfailxpi = min(-Ks*xin2(jj), Mu) ;

    phifailxpi = abs(tM)-(Mu-qfailxpi) ;
    
        if phifailxpi <= 0

            alfan1 = alfan ;
            xin21xpi = xin2(jj) ;

            else
                gamma2 = piecewise(xin2(jj)<=-Mu/Ks, phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks),abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2))) ;
                alfan1=alfan + gamma2*sign(tM) ;
                xin21xpi = xin2(jj) + gamma2 ;

                xd = xpi ;       
        end
    end

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
end