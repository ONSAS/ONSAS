clear ;
close all ;

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

nu = 0.3 ;

% --- material constit params ---
hyperElasParams(1) = E ;
hyperElasParams(2) = nu ;
G   = E/(2*(1+nu)) ;

[ dn1, kpn1, xin11, xin21, alfan1, xd] = framePlastic( dn, kpn, xin1, xin2, alfan, xd, elemCrossSecParams, density, hyperElasModel, hyperElasParams) ;
