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

function [ kpn1, xin11, xin21, alfan1, xd ] = FramePlastic( dn, kpn, xin1, xin2, alfan, xd, elemCoords, elemCrossSecParams, massMatType, density, hyperElasModel, hyperElasParams, Ut, Udotdotte, intBool, matFintBool, elem )
  
  ndofpnode = 6 ;
  
  % --- material constit params ---
  E   = hyperElasParams(1) ;
  nu  = hyperElasParams(2) ;
  G   = E/(2*(1+nu)) ;

  [A, J, Iy, Iz] = crossSectionProps ( elemCrossSecParams, density ) ;

  % numerical example
  % cantilever beam of rectangular cross-section loaded with a vertical force at the free end
  
  % Mc, My, Mu / from the moment-curvature diagram
  % kh1, kh2   / hardening modules
  % Ks         / from the moment-rotation jump diagram

  l = 2.5 ;         % m
  E = 300000000 ;   % K(N/m^2) KPa
  EI = 77650 ;      % KNm^2
  Iy = EI/E ;       % m^4
  Mc = 37.9 ;       % KNm
  My = 268 ;
  Mu = 374 ;
  kh1 = 29400 ;     % KNm^2
  kh2 = 272 ;
  Ks = -18000 ;     % KNm

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
	
  % --- elem lengths and rotation matrix
  [ local2globalMats, l ] = beamParameters( elemCoords ) ;
  R = RotationMatrix(ndofpnode, local2globalMats) ;

  % --- set the local degrees of freedom corresponding to each behavior
  LocAxialdofs  = [ 1 7 ] ;
  LocTorsndofs  = [ 2 8 ] ;
  LocBendXYdofs = [ 3 6 9 12 ] ;
  LocBendXZdofs = [ 5 4 11 10 ] ;

  KL = zeros ( 2*ndofpnode, 2*ndofpnode ) ;

  Kaxial = E*A/l * [ 1 -1 ; ...
                    -1  1 ] ;
  KL( LocAxialdofs , LocAxialdofs ) = Kaxial ;

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

  % uvector = [u1; u2] ;
  % vvector = [v1; v2] ;
  % thetavector = [theta1; theta2] ;

  uvector = dn(1,:)' ;
  vvector = dn(2,:)' ;
  thetavector = dn(3,:)' ;
  
  Bu = [-1/l 1/l] ;

  Kfd = zeros(6,6) ;
  Kfalfa = zeros(6,6) ;
  Khd = zeros(6,6) ;
  Khalfa = zeros(6,6) ;

  % Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]

  npi = 3 ;
  xpi = [a (a+b)/2 b] ;
  wpi = [1/3 4/3 1/3] ;

  kpn1 = zeros(npi,1) ;
  xin11 = zeros(npi,1) ;
  xin21 = zeros(npi,1) ;

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
  % integration (Gauss-Lobatto)
  % and calculation of values of internal parameters at integration points

  for j = 1:npi

    [Kfdj, Kfalfaj, Khdj, Khalfaj, kpn1xpi, xin11xpi, xin21xpi, M1xpi, tM, xd] = integrand(j, xpi(j), xd) ;

    % integration (Gauss-Lobatto)

    Kfd = Kfd + Kfdj*wpi(j) ;
    Kfalfa = Kfalfa + Kfalfaj*wpi(j) ;
    Khd = Khd + Khdj*wpi(j) ;
    Khalfa = khalfa + Khalfaj*wpi(j) ;

    % values of internal parameters at integration points

    kpn1(j) = kpn1xpi ;
    xin11(j) = xin11xpi ;
    xin21(j) = xin21xpi ;

  end

  Khalfa = Khalfa + Ks ; % integral + ks

  function [Kfd, Kfalfa, Khd, Khalfa, kpn1xpi, xin11xpi, xin21xpi, M1xpi, tM, xd] = integrand(j, xpi, xd)

    % elastoplasticity with hardening
    % the usual trial-corrector (return mapping) algorithm
    % is used at each of the Gauss–Lobatto integration points.

    N = bendingInterFuns (xpi, l, 2) ;
    Bv = [N(1) N(3)] ;
    Btheta = [N(2) N(4)] ;

    Bd = [Bu 0 0; 0 Bv Btheta] ;

    Ghat = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi/l)) ;

    % curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)

    khat = Bv*vvector + Btheta*thetavector + Ghat*alfan ;
    khat2 = dirac(xd)*alfan ;
    kn = khat + khat2 ;
    ken = khat - kpn(j) ;

    % moment

    Mxpi= E*Iy*ken ;

    % yield criterion

    qxpi = piecewise(xin1(j) <= (My-Mc)/kh1, -kh1*xin1(j), -(My-Mc)*(1-kh2/kh1)-kh2*xin1(j)) ;
    phixpi = abs(Mxpi) - (Mc - qxpi) ;

    % test values

    % kpn1test = kpn ;
    % xin11test = xin1 ;
    phitest =  phixpi ;

    % gamma values calculations (gamma derivative is the plastic multiplier)
    % the new values of internal variables are computed

    if phitest <= 0
  
        gamma = 0 ;
        kpn1xpi = kpn(j) ;
        xin11xpi = xin1(j) ;
        M1xpi = Mxpi ;

    else

        gamma = piecewise(xin1(j) + phitest/(kh1+E*I)<=(My-Mc)/kh1, phitest/(kh1+E*I), phitest/(kh2+E*I)) ;
        kpn1xpi = kpn(j) + gamma*sign(M) ;
        xin11xpi = xin1(j) + gamma ;
        M1xpi = E*Iy*(khat-kpn1(j)) ;

    end

    % elastoplastic tangent bending modulus

    Cep = piecewise(gamma == 0, E*Iy , gamma > 0 & xin11xpi <= (My-Mc)/kh1, E*Iy*kh1/(E*Iy + kh1), gamma > 0 & xin11xpi > (My-Mc)/kh1, E*Iy*kh2/(E*Iy + kh2)) ;
  
    % stiffness matrices

    Kfd = Bd'*[E*A 0; 0 Cep]*Bd ;

    Kfalfa = Bd'*[E*A 0; 0 Cep]*[0 Ghat]' ;

    Khd = [0 Ghat]*[E*A 0; 0 Cep]*Bd ;

    Khalfa = Ghat*Cep*Ghat ;

    % plastic softening at the discontinuity
    % the standard trial-corrector (return mapping) algorithm is used also for softening rigid plasticity

    % test values
  
    % alfan1test = alfan ;
    % xin2test = xin2 ;

    % trial value of the moment at the discontinuity (tM at xd)
    tM = 0 ;

    for j=1:npi

        Ghatxpi = -1/l*(1+3*(1-2*xd/l)*(1-2*xpi(j)/l)) ;

        % curvatures (time n) / k, ke, kp, khat (continuous part of the curvature), khat2 (localized part of the curvature)

        khat = Bv*vvector + Btheta*thetavector + Ghatxpi*alfan ;
        khat2 = dirac(xd)*alfan ;
        kn = khat + khat2 ;
        kenxpi = khat - kpn(j) ;
        
        % moment at integration points
        Mixpi = E*Iy*kenxpi ;

        % integration (Gauss-Lobatto)
    
        tM = tM - Ghatxpi*Mixpi*wpi(j) ;

    end

    % softening criterion (failure function) at integration points

    qfailxpi = min(-Ks*xin2(j), Mu) ;
    phifailxpi = abs(tM)-(Mu-qfail) ;
    
        if phifailxpi <= 0

            alfan1 = alfan ;
            xin21xpi = xin2(j) ;

            else
                gamma2 = piecewise(xin2(j)<=-Mu/Ks, phifailxpi/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2)+Ks),abs(tM)/((4*E*Iy)/l^3*(l^2-3*l*xd+3*xd^2))) ;
                alfan1=alfan + gamma2*sign(tM) ;
                xin21xpi = xin2(j) + gamma2 ;

                xd = xpi ;       
        end
  end
end