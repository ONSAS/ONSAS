% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan

% =========================================================================

% numerical example
% cantilever beam loaded with a vertical force at the free end

% =========================================================================

close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src'] ) );
          
% Mc, My, Mu / from the moment-curvature diagram
% kh1, kh2   / hardening modules
% Ks         / from the moment-rotation jump diagram
        
l = 2.5 ;           % m
A = 0.4*0.3 ;       % m^2
E = 30000000 ;      % KN/m^2 KPa
EI = 77650 ;        % KN.m^2
Iy = EI/E ;         % m^4
Mc = 37.9 ;         % KN.m
My = 268 ;
Mu = 374 ;
kh1 = 29400 ;       % KN.m^2
kh2 = 272 ;
Ks = -18000 ;       % KN.m

freedofs = [2 4 6]; % u2 v2 theta2

nu   = 0.3 ;
tol1 = 1e-8 ;
tol2 = 1e-4 ;
tolk = 3 ;

% initial values
dn   = [0 0 0 0 0 0]' ;
Fint = 0 ;
tM   = 0 ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;

kpn  = zeros(npi,1) ;
xin1 = zeros(npi,1) ; % xi en tiempo n
xin2 = zeros(npi,1) ;

xin11 = zeros(npi,1) ; % xi del instante de tiempo n+1
xin21 = zeros(npi,1) ; % xi cuando hay rotula en n+1

alfan = 0 ;

xd = 0 ;

%Final_force = 120;
Final_force = 30;

load_case = [0 0 0 -1 0 0]' ; % load applied in vertical direction (Y)  [ u1 u2 v1 v2 t1 t2 ]
load_factors = 0:Final_force ;

% --- element params ---
elemParams = [l A Iy] ;

% --- elastoplastic params ---
elastoplasticParams = [E Mc My Mu kh1 kh2 Ks] ;

matdes = zeros (6, Final_force+1) ;

matdes(:,1) = dn ;

gxin = zeros(Final_force,1) ;

% header
fprintf(' timeInd  iteration  deltaU  residForce xin1 \n');

for ind = 2:length(load_factors)

    curr_load_factor = load_factors(ind) ;

    Fext = load_case * curr_load_factor ;

    dnk = matdes(:,ind-1) ;

    % iteration vars
    converged_bool = false ;
    k = 0 ; % set iterations zero

    # xin1 = xin11 ; % setea el valor del n tomando el de tiempo n+1 
    gxin(ind-1,1) = xin1(1) ;

    # fprintf('xin1 \n') ;
    # disp(xin11) ;

    fprintf('-------------------------------------------------\n') ;


    xin11 = xin1 ;

    while converged_bool == false && k < tolk

        k = k + 1 ;

        # fprintf('= = = = Iteración %d = = = =\n', k) ;



        [Fint, Kelement, kpn1, xin11, xin21, alfan1, xd, tM] = framePlastic(dnk, kpn, xin11, xin2, alfan, xd, elemParams, elastoplasticParams) ;
        
        # fprintf('Fint \n') ;
        # disp(Fint)
        # fprintf('Fext \n') ;
        # disp(Fext) ;
        
        residualForce = Fext - Fint ;

        Krelement = Kelement( freedofs, freedofs) ;
 
        residualForceRed = residualForce(freedofs) ;
      
        % system of equilibrium equations
        # fprintf('Krelement \n') ;
        # disp(Krelement) ;
        # fprintf('residualForceRed \n') ;
        # disp(residualForceRed) ;
        deltadred = Krelement\ (  residualForceRed ) ;

        %fprintf('deltadred \n') ;
        %disp(deltadred) ;
        
        deltad = zeros(6,1) ;      
        deltad(freedofs) = deltadred ;
        
        %fprintf('deltad \n') ;
        %disp(deltad) ;

        dnk1 = dnk + deltad ;

        # fprintf('dnk1 \n') ;
        # disp(dnk1) ;

        dnk = dnk1 ;

        norm1 = norm( deltadred ) ;
        norm2 = norm( residualForceRed ) ;

        fprintf('  %4i  %3i  %12.4e  %12.4e  %12.4e   \n', curr_load_factor, k, norm(deltadred), norm(residualForceRed), xin11(1) ) ;

        converged_bool = norm1 < tol1 || norm2 < tol2 ;

    end
    xin1 = xin11 ; % setea el valor del n tomando el de tiempo n+1 

    matdes(:, ind) = dnk1 ;

end

lw = 2.5 ; ms = 0.5 ; plotfontsize = 16 ;

figure('Name','Cantilever Beam / Plasticity','NumberTitle','off');
hold on, grid on
plot(abs(matdes(6,:)), load_factors,'b-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(matdes(4,:)), load_factors, 'k-o' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;
labx = xlabel('Generalized displacements (m, rad)');   laby = ylabel('Load Factor \lambda (KN)') ;
legend('Degree of Freedom y','Degree of Freedom \theta','location','Southeast') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity') ;

figure('Name','Cantilever Beam / Plasticity','NumberTitle','off');
hold on, grid on
plot(gxin, load_factors(1:(length(load_factors)-1))*2.5, 'k-o' , 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
labx = xlabel('Curvature accumulated \xi');   laby = ylabel('Moment applied (KN.m)') ;
legend('Internal parameter of plasticity \xi','location','Southeast');
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity') ;