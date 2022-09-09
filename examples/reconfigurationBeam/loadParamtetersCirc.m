
% #md Reconfiguration problem validation (Drag reduction of flexible plates by reconfiguration, Gosselin, etAl 2010)
%----------------------------
function [L, w, EI, rhoA, nuA, NR, cy_vec, uy_vec  ] = loadParamtetersCirc()
    % values extracted from https://github.com/lm2-poly/Reconfiguration-Beam/blob/master/reconfiguration.m
    NR      =100;  %NUMBER OF CYCD POINTS
    cycdmin   =-2 ;  %SMALLEST VALUE OF CYCD = 10^cymin
    cycdmax   =5  ;  %LARGEST VALUE OF CYCD = 10^cymax with cy = rho * L^3 * U^2 / 16 EI 
    cycdvec_indexes   = linspace(cycdmin, cycdmax, NR) ;
    % Geometric params
    l  = 1         ;
    d  = l/100     ;
    J = pi * d ^ 4 / 64 ; Iyy = J / 2 ; Izz = J / 2 ;  
    % Material params
    E = 3e7 ;  nu = 0.3 ; rho = 700 ; G = E / (2 * (1+nu)) ; B = E*Izz;
    % Fluid properties
    rhoF = 1020 ; nuA = 1.6e-5 ; c_d = feval('dragCircular', 0, 0) ;
    % fill U vector for each cy 
    uy_vec = [] ;
    cycd_vec = [] ;
    for cy_index = cycdvec_indexes
        cycd = 10^cy_index          ;   
        cycd_vec = [cycd_vec, cycd]   ;
        uy_vec =[uy_vec, sqrt(cycd * 4 * B / (c_d * l^3 * d * rhoF))];
    end


end