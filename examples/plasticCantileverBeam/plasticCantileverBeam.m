% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% =========================================================================

% numerical example
% cantilever beam loaded with a vertical force at the free end

% =========================================================================

close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src'] ) );
          
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

freedofs = [2 4 6 ];

nu = 0.3 ;
tol1 = 1e-8 ;
tol2 = 1e-4 ;
tolk = 10 ;

% initial values
dn = [0 0 0 0 0 0]' ;
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

matdes = dn ;

%Mc / l 
for n = 2:100

    fprintf('Fuerza  %d \n',n) ;

    k = 0 ; % set iterations zero

    Fext = [0 0 0 n-1 0 0]' ;

    dnk = matdes(:,n-1) ;

    nonconverge = true ;

    while nonconverge && k < tolk

        k = k + 1 ;

        fprintf('Iteración %d\n',k) ;

        [Fint, Kelement, kpn1, xin11, xin21, alfan1, xd, tM] = framePlastic(dnk, kpn, xin1, xin2, alfan, xd, Fint, tM, elemParams, elastoplasticParams) ;
        
        residualForce = Fint - Fext ;

        Krelement = Kelement( freedofs, freedofs) ;
 
        residualForceRed = residualForce(freedofs) ;
      
        % system of equilibrium equations

        deltadred = Krelement\ ( - residualForceRed ) ;
        
        deltad = zeros(6,1);      
        deltad(freedofs) = deltadred ;
            
        dnk1 = dnk + deltad ;

        dnk = dnk1 ;

        norm1 = norm( deltadred ) ;
        norm2 = norm( residualForceRed ) ;

        nonconverge = norm1 > tol1 &&  norm2 > tol2 ;

    end

    matdes = [matdes dnk1] ;

end