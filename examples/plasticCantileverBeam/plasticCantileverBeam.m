clear ;
close all ;

% numerical example
% cantilever beam of rectangular cross-section loaded with a vertical force at the free end
          
% Mc, My, Mu / from the moment-curvature diagram
% kh1, kh2   / hardening modules
% Ks         / from the moment-rotation jump diagram
        
l = 2.5 ;         % m
A = 0.4*0.3 ;     % m^2
E = 300000000 ;   % K(N/m^2) KPa
EI = 77650 ;      % KNm^2
Iy = EI/E ;       % m^4
Mc = 37.9 ;       % KNm
My = 268 ;
Mu = 374 ;
kh1 = 29400 ;     % KNm^2
kh2 = 272 ;
Ks = -18000 ;     % KNm

nu = 0.3 ;
tol1 = 0.01 ;
tol2 = 0.01 ;

% initial values
dn = [0 0 0 0 0 0] ;
dn1 = [ 1 1 1 1 1 1] ;
Fint = 0 ;
tM = 0 ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]

npi = 3 ;

kpn = zeros(npi,1) ;
xin1 = zeros(npi,1) ;
xin2 = zeros(npi,1) ;

alfan = 0 ;

xd = 0 ;

% --- element params ---
elemParams = [l A Iy] ;

% --- elastoplastic params ---
elastoplasticParams = [E Mc My Mu kh1 kh2 Ks] ;

for lambda = 1:1000

while norm(dn1-dn) > tol1 && abs(lambda - Fint) > tol2

[dn1, kpn1, xin11, xin21, alfan1, xd, Fint, tM] = framePlastic(dn, kpn, xin1, xin2, alfan, xd, Fint, tM, elemParams, elastoplasticParams, lambda) ;

end

end