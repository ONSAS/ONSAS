% #md Reconfiguration problem validation (Drag reduction of flexible plates by reconfiguration, Gosselin, etAl 2010)
%----------------------------
function [l, d, Izz, E, nu, rhoS, rhoF, nuF, dragCoefFunction, NR, cycd_vec, uy_vec ] = loadParametersCirc()
    % values extracted from https://github.com/lm2-poly/Reconfiguration-Beam/blob/master/reconfiguration.m
    NR      =10;  %NUMBER OF CYCD POINTS
    cycdmin =-2 ;  %SMALLEST VALUE OF CYCD = 10^cymin
    cycdmax =5  ;  %LARGEST VALUE OF CYCD = 10^cymax with cy = rho * L^3 * U^2 / 16 EI 
    cycdvec_indexes   = linspace(cycdmin, cycdmax, NR) ;
    % Geometric params
    l  = 1         ;
    d  = l/100     ;
    Izz = pi * d ^ 4 / 64 ;  
    % Material params
    E = 3e7 ;  nu = 0.3 ; rhoS = 700 ; B = E*Izz;
    % Fluid properties
    dragCoefFunction = 'dragCircular';
    rhoF = 1.225 ; nuF = 1.6e-5 ; c_d = feval(dragCoefFunction, 0, 0) ;
    % fill U vector for each cy 
    uy_vec = [] ;
    cycd_vec = [] ;
    for cy_index = cycdvec_indexes
        cycd = 10^cy_index          ;   
        cycd_vec = [cycd_vec, cycd]   ;
        uy_vec =[uy_vec, sqrt(cycd * 2 * B / (c_d * l^3 * d * rhoF))];
    end
end