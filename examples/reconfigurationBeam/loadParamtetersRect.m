
% #md Reconfiguration problem validation (Drag reduction of flexible plates by reconfiguration, Gosselin, etAl 2010)
%----------------------------
function [L, w, EI, rhoA, nuA, NR, cy_vec, uy_vec  ] = loadParamteters()
    % values extracted from https://github.com/lm2-poly/Reconfiguration-Beam/blob/master/reconfiguration.m
    NR      =100;  %NUMBER OF CYCD POINTS
    cymin   =-2 ;  %SMALLEST VALUE OF CYCD = 10^cymin
    cymax   =5  ;  %LARGEST VALUE OF CYCD = 10^cymax with cy = rho * L^3 * U^2 / 16 EI 
    cyvec_indexes   = linspace(cymin, cymax, NR) ;
    % values extracted from TABLE 2 Specimen 1 
    L  = 15.8*1e-2 ;
    w  = 3.8*1e-2  ;
    EI = 1824*1e-6 ;
    % fluid properties
    rhoA = 1.225 ; nuA = 1.6e-5 ; 
    % fill U vector for each cy 
    uy_vec = [] ;
    cy_vec = [] ;
    for cy_index = cyvec_indexes
        cy = 10^cy_index        ;   
        cy_vec = [cy_vec, cy]   ;
        uy_vec =[uy_vec, sqrt(cy * 2 * EI / (rhoA * (L/2)^3)) ];
    end

end