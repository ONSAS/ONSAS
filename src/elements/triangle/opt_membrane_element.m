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

% Implementation based on paper: A study of optimal membrane triangles with drilling freedoms - Carlos A. Felippa
% DOI: doi:10.1016/S0045-7825(03)00253-6

function [Km, Kb, Kh] = opt_membrane_element(E, nu, h, A, r1, r2, r3)

    % geometric properties
    %
    % coordinate differences
    xy12 = r1 - r2;
    xy21 = -xy12;

    xy23 = r2 - r3;
    xy32 = -xy23;

    xy31 = r3 - r1;
    xy13 = -xy31;

    r0 = ( r1 + r2 + r3 ) / 3;
    xy10 = r1 - r0;
    xy20 = r2 - r0;
    xy30 = r3 - r0;
    % side lengths
    l12 = norm(xy12);
    l23 = norm(xy23);
    l31 = norm(xy31);
    % triangle heights
    a12 = 2 * A / l12;
    a23 = 2 * A / l23;
    a31 = 2 * A / l31;
    % median lengths
    m1 = 3/2 * sqrt(norm(xy10));
    m2 = 3/2 * sqrt(norm(xy10));
    m3 = 3/2 * sqrt(norm(xy10));
    % side lengths projected normal to medians
    b1 = 2 * A / m1;
    b2 = 2 * A / m2;
    b3 = 3 * A / m3;

    % E mat
    Eaux = [ 1  nu  0   ;...
            nu  1   0   ;...
            0   0   (1 - nu) / 2];
    Emat = E / (1 - nu^2) * Eaux;
    % Natural strains
    %
    % Straingage rosette transformation
    y12 = xy12(2);
    y21 = xy21(2);

    y23 = xy23(2);
    y32 = xy32(2);

    y13 = xy13(2);
    y31 = xy31(2);

    x12 = xy12(1);
    x21 = xy21(1);

    x23 = xy23(1);
    x32 = xy32(1);

    x13 = xy13(1);
    x31 = xy31(1);

    Maux = [    y23*y13*l12^2       y31*y21*l23^2    y12*y32*l31^2 ;...
                x23*x13*l12^2       x31*x21*l23^2    x12*x32*l31^2 ;...
                (y23*x31 + x32*y13)*l12^2    (y31*x12 + x13*y21)*l23^2     (y12*x23 + x21*y32)*l31^2 ];
    Te =  1 / (4 * A^2) * Maux ;

    Enat = Te' * Emat * Te; % eq 13

    % Hierarchical rotations
    %
    Maux = [ x32   y32    4*A x13   y13  0       x21    y21 0 ;...
             x32   y32      0 x13   y13  4*A     x21    y21 0 ;...
             x32   y32      0 x13   y13  0       x21    y21 4*A];
    Tthetau = 1 / (4 * A) * Maux; 

    % Stiffness template
    %
    % Basic stiffness
    %
    V = A * h;
    alpha_b = 3/2; % OPT
    Laux = [ y23                            0                               x32 ;...
            0                               x32                             y23 ;...
            1/6*alpha_b*y23*(y13-y21)       1/6*alpha_b*x32*(x31-x12)       1/3*alpha_b*(x31*y13 - x12*y21) ;...
            y31                             0                               x13 ;...
            0                               x13                             y31 ;...
            1/6*alpha_b*y31*(y21-y32)       1/6*alpha_b*x13*(x12-x23)       1/3*alpha_b*(x12*y21 - x23*y32) ;...
            y12                             0                               x21 ;...
            0                               x21                             y12 ;...
            1/6*alpha_b*y12*(y32-y13)       1/6*alpha_b*x21*(x23-x31)       1/3*alpha_b*(x23*y32-x31*y13) ];
    
    L = 1/2 * h * Laux; % eq 22
    Kb = 1 / V * L * Emat * L'; % eq 21

    % Higher order stiffness
    %
    % OPT free parameters
    beta0 = max( 1/2 * (1 - 4 * nu^2), 0.01 );
    beta1 = 1;
    beta2 = 2;
    beta3 = 1;
    beta4 = 0;
    beta5 = 1;
    beta6 = -1;
    beta7 = -1;
    beta8 = -1;
    beta9 = -2; 

    % eq 29
    Q1 = 2 * A / 3 * [  beta1/l12^2 beta2/l12^2 beta3/l12^2 ;...
                        beta4/l23^2 beta5/l23^2 beta6/l23^2 ;...
                        beta7/l31^2 beta8/l31^2 beta9/l31^2 ];

    Q2 = 2 * A / 3 * [  beta9/l12^2 beta7/l12^2 beta8/l12^2 ;...
                        beta3/l23^2 beta1/l23^2 beta2/l23^2 ;...
                        beta6/l31^2 beta4/l31^2 beta5/l31^2 ];
    
    Q3 = 2 * A / 3 * [  beta5/l12^2 beta6/l12^2 beta4/l12^2 ;...
                        beta8/l23^2 beta9/l23^2 beta7/l23^2 ;...
                        beta2/l31^2 beta3/l31^2 beta1/l31^2 ]; % diferencia con paper, parece que la eq 29 tiene un error!!

    % eq 30
    Q4 = 0.5 * (Q1 + Q2);
    Q5 = 0.5 * (Q2 + Q3);
    Q6 = 0.5 * (Q3 + Q1);

    Ktheta = 3/4 * beta0 * h * A * ( Q4' * Enat * Q4 + Q5' * Enat * Q5 + Q6' * Enat * Q6 ); % eq 31
    Kh =  Tthetau' * Ktheta * Tthetau ;

    % Membrane stiffness
    Km = Kb + Kh ;

end