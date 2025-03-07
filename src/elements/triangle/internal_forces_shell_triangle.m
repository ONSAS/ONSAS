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
%
% Implementation of a triangular finite element with 6 dfos (3 translations and 3 rotations) per node for the analysis of linear elastic isotropic shells with constant thickness.
% The element is formed by the superposition of a plate element (DKT) and a plane stress element (CST) with
% with addition to artificial drilling (rotation about the axis normal to the element plane) stiffness.
%
function [ fs, ks, fintLocCoord ] = internal_forces_shell_triangle(elemCoords, elemDisps, modelName, modelParams, thickness )


    %material and geometric parameters
    young_modulus = modelParams(1);  
    poisson_ratio = modelParams(2) ;
    h = thickness;

    elemCoords = elemCoords';

    r1g = elemCoords(1:3);
    r2g = elemCoords(4:6);
    r3g = elemCoords(7:9);
    rog = (r1g + r2g + r3g) / 3;

    Ug = switchToTypeIndexing( elemDisps ) ;

    u1g = Ug( 1: 3);
    q1  = Ug( 4: 6);
    u2g = Ug( 7: 9);
    q2  = Ug(10:12);
    u3g = Ug(13:15);
    q3  = Ug(16:18);

    p1g = r1g + u1g;
    p2g = r2g + u2g;
    p3g = r3g + u3g;
    pog = (p1g + p2g + p3g) / 3;

    % eq. (35) of 10.1016/j.cma.2006.10.006
    R1g = global_rotation_matrix(q1);
    R2g = global_rotation_matrix(q2);
    R3g = global_rotation_matrix(q3);

    [To, x02, x03, y03] = edge_local_axis_shell_triangle(r1g,r2g,r3g);
    [Tr, x02, x03, y03] = edge_local_axis_shell_triangle(p1g,p2g,p3g);

    Ro = To';
    Rr = Tr';

    # To
    # r1g
    # rog
    r1o = To*(r1g - rog);
    r2o = To*(r2g - rog);
    r3o = To*(r3g - rog);

    % eq. (1) of 10.1016/j.cma.2006.10.006
    u1def = Tr*(p1g - pog) - r1o;
    u2def = Tr*(p2g - pog) - r2o;
    u3def = Tr*(p3g - pog) - r3o;

    % eq. (2) of 10.1016/j.cma.2006.10.006
    R1def = Tr*R1g*Ro;
    R2def = Tr*R2g*Ro;
    R3def = Tr*R3g*Ro;

    % eq. (13) of 10.1016/j.cma.2006.10.006
    v1def = rotation_vector(R1def);
    v2def = rotation_vector(R2def);
    v3def = rotation_vector(R3def);

    % eq. (12) of 10.1016/j.cma.2006.10.006
    pl = zeros(15,1);
    pl(1:2) = u1def(1:2);
    pl(3:5) = v1def;
    pl(6:7) = u2def(1:2);
    pl(8:10) = v2def;
    pl(11:12) = u3def(1:2);
    pl(13:15) = v3def;

    pl_full = zeros(18,1);
    index_full = [1,2, 4,5,6, 7,8, 10,11,12, 13,14, 16,17,18];
    pl_full(index_full) = pl;

    % calculating the stiffness matrix and internal force vector of the shell element in local coordinates
    [Kl_full, fintLocCoord] = local_shell_triangle(x02, x03, y03, young_modulus, poisson_ratio, h, pl_full);
    index_full = [1,2, 4,5,6, 7,8, 10,11,12, 13,14, 16,17,18];
    % reducing to 15 dofs
    Kl = Kl_full( index_full, index_full);

    % local internal force vector
    fl = Kl * pl;

    % eq. (15) of 10.1016/j.cma.2006.10.006
    Ta1 = matrix_Ta(R1def);
    Ta2 = matrix_Ta(R2def);
    Ta3 = matrix_Ta(R3def);

    % eq. (18) of 10.1016/j.cma.2006.10.006
    fa = fl ;
    fa(3:5) = Ta1 * fl(3:5);
    fa(8:10) = Ta2 * fl(8:10);
    fa(13:15) = Ta3 * fl(13:15);

    % eq. (19) of 10.1016/j.cma.2006.10.006
    Ba = eye(15);
    Ba(3:5,3:5) = Ta1;
    Ba(8:10,8:10) = Ta2;
    Ba(13:15,13:15) = Ta3;

    % eq. (21) of 10.1016/j.cma.2006.10.006
    Kh = zeros(15,15);
    Kh( 3: 5, 3: 5) = matrix_Khi(R1def, fl( 3: 5));
    Kh( 8:10, 8:10) = matrix_Khi(R2def, fl( 8:10));
    Kh(13:15,13:15) = matrix_Khi(R3def, fl(13:15));

    % eq. (20) of 10.1016/j.cma.2006.10.006
    % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
    Ka = Ba' * Kl * Ba + Kh;

    # r1o
    # u1def

    % eq. (7) of 10.1016/j.cma.2006.10.006
    a1 = u1def + r1o;
    a2 = u2def + r2o;
    a3 = u3def + r3o;

    % eq. (27) of 10.1016/j.cma.2006.10.006
    [G1,G2,G3] = matrix_Gi(a1, a2, a3, r1o, r2o, r3o);
    G = [G1; G2; G3];

    # a1
    # a2
    # G1
    # G2
    # G3

    % eq. (26) of 10.1016/j.cma.2006.10.006
    P = matrix_P(a1, a2, a3, G1, G2, G3);

    % eq. (25) of 10.1016/j.cma.2006.10.006
    E = blkdiag(Rr,Rr,Rr,Rr,Rr,Rr);

    % eq. (30) of 10.1016/j.cma.2006.10.006
    # P
    # fa
    n = P' *fa; % eq. (31) of 10.1016/j.cma.2006.10.006
    [F1, F2] = matrix_F(n);

    % eq. (29) of 10.1016/j.cma.2006.10.006
    % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
    # E
    # n
    fg = E * n;
    Kg = E * ( P' * Ka * P  - G*F1'*P - F2*G') * E';

    % eq. (39) of 10.1016/j.cma.2006.10.006
    Bm = eye(18);
    Bm( 4: 6, 4: 6) = matrix_Tm(q1);
    Bm(10:12,10:12) = matrix_Tm(q2);
    Bm(16:18,16:18) = matrix_Tm(q3);

    % eq. (40) of 10.1016/j.cma.2006.10.006
    Kk = zeros(18,18);
#     q1
# fg(4:6)
    Kk( 4: 6, 4: 6) = matrix_Kki(q1, fg( 4: 6));
    Kk(10:12,10:12) = matrix_Kki(q2, fg(10:12));
    Kk(16:18,16:18) = matrix_Kki(q3, fg(16:18));

    % eq.(38) of 10.1016/j.cma.2006.10.006
    % this could be done much more efficiently avoiding unnecessary multiplications by zero or 1
    fm = Bm' * fg;
    Km = Bm' * Kg * Bm + Kk;

    % shifting lines and columns to onsas convention of dofs ordering
    ks = {switchToNodalIndexing( Km )}; 
    fs = {switchToNodalIndexing( fm )};

end


function [Kel,fintLocCoord] = local_shell_triangle(x02, x03, y03, E, nu, h, Ul);

    % Calculate the area of the triangle
    area = x02 * y03 / 2;

    % membrane stiffness
    aux1 = h *  E / ( 1 - nu^2) ; 
    aux2 = nu*aux1;
    Dm = [ [aux1, aux2 , 0 ]; [aux2, aux1, 0] ; [0, 0, aux1*(1-nu)/2] ];
    Bm = CST_B(x02, x03, y03);
    Km = area * Bm' * Dm * Bm ;
    im = [1, 2, 7, 8, 13, 14];
    % membrane forces (constant)
    Ulm = Ul(im);
    N = Dm * Bm * Ulm;

    % bending stiffness
    aux1 = E * h^3 / (12 * ( 1- nu^2) ); 
    aux2 = nu*aux1;
    Db = [ [aux1, aux2 , 0 ]; [aux2, aux1, 0] ; [0, 0, aux1*(1-nu)/2] ];
    ib = [3,4,5, 9,10,11, 15,16,17];
    int_point = [ [1./6. 1./6.]; [2./3., 1./6.]; [1./6., 2./3.]];
    Kb = zeros(9,9);
    wgt = area / 3.0;
    for ipt = 1:3;
        psi = int_point(ipt,1);
        eta = int_point(ipt,2);
        Bb = DKT_B(psi, eta, x02, x03, y03);
        Kb = Kb + wgt * Bb' * Db * Bb;
    end
    % plate moments (calculated at the center of the element)
    Ulb = Ul(ib);
    Bb = DKT_B(1./3., 1./3., x02, x03, y03);
    curv = Bb * Ulb;
    M = Db * curv;    

    % assembling the stiffness matrix of the shell element in local coordinates
    Kel = zeros(18,18);

    Kel(im,im) = Km;
    Kel(ib,ib) = Kb;
    
    k_dr = min( min( abs( Kb ) ) ) * 1.e-4;
    Kel(6 , 6) = k_dr;
    Kel(12,12) = k_dr;
    Kel(18,18) = k_dr;

    % returning moments in local element coordiantes 
    fintLocCoord = zeros(1,3);
    fintLocCoord = fintLocCoord + M;  % Store the moments in fintLocCoord
end

function [T, x02, x03, y03] = edge_local_axis_shell_triangle(p1,p2,p3);
%Calculates the matrix for transformation of basis between global and local axis;
%p1, p2 and p3 are the position vector for the nodes in global coordinates;
%the local x axis is paralel to the side connecting nodes 1 and 2
%the local z axis is normal to the element plane;
%the origin of local axis is located at node 1
    
    p12 = p2 - p1;
    p13 = p3 - p1;
    
    au_zl =  cross(p12,p13);
    # stop
    u_zl = au_zl / norm(au_zl);
    
    x02 = norm(p12);
    u_xl = p12 / x02;

    u_yl = cross(u_zl, u_xl);

    T = [ u_xl u_yl u_zl]';

    x03 = dot(u_xl, p13);
    y03 = dot(u_yl, p13);
end


function [ B ] = CST_B(x02, x03, y03);
    %calculate the strain-displacement matrix for the constant stress triangular element (CST)
    %x02, x03 and y03 are the local coordinates of the nodes 2 and 3 

    area02 = x02*y03;

    %( with y01 = y02 = x01 = 0)
    %bi = (yj - yk )/2A
    b1 = - y03/area02;
    b2 =   y03/area02;
    b3 =  0;

    %ci = (xk - xj)/2A
    c1 = (x03 - x02)/area02;
    c2 = - x03/area02;
    c3 = x02/area02;

    B = [   [b1,  0, b2,  0, b3,  0];
            [ 0, c1,  0, c2,  0, c3];
            [c1, b1, c2, b2, c3, b3]; ];

end

function [ B ] = DKT_B(PSI, ETA, x02, x03, y03);
    %calculate the strain-displacement matrix for the triangular plate element (DKT)
    %psi and eta are the area coordinates of the triangular element
    %x02, x03 and y03 are the local coordinates of the nodes 2 and 3 

    area02 = x02*y03;

    %auxiliar variables - geometric
    X32 = x03 - x02;
    XX2 = x02 * x02;
    XX3 = x03 * x03;
    XX32= X32 * X32;
    YY3 = y03 * y03;
    QQ4 = 1.d0 / (XX32 + YY3);
    QQ5 = 1.d0 / (XX3 + YY3);
    QQ6 = 1.d0 / XX2;
    SU4 = 3.d0 * QQ4 * YY3;
    SU5 = 3.d0 * QQ5 * YY3;
    D4  = 3.d0 * QQ4 * X32 * y03;
    D5  = 3.d0 * QQ5 * x03 * y03;
    CL4 = 6.d0 * QQ4 * X32;
    CL5 = -x03 * 6.d0 * QQ5;
    CL6 = 6.d0 * QQ6 * x02 ;
    SL4 = 6.d0 * QQ4 * y03;
    SL5 = -y03 * 6.d0 * QQ5;

    %auxiliar variables - coordinates
    PS2  = 1.d0 -2.d0 * PSI;
    ET2  = 1.d0 -2.d0 * ETA;
    SU5E = SU5 * ET2;

    %Derivada de betax com relacao a psi
    BXP1 = CL6 * PS2 + ETA*(CL5-CL6);
    BXP2 = -ETA * D5;                                                    
    BXP3 = -4.d0 + 6.d0 * (PSI+ETA) - ETA * SU5;                               
    BXP4 = -CL6 * PS2 + ETA * (CL4+CL6);                                     
    BXP5 = ETA * D4 ;                                                    
    BXP6 = -2.d0 + 6.d0 * PSI + ETA * SU4    ;                                 
    BXP7 = -ETA * (CL5+CL4)               ;                              
    BXP8 = ETA * (D4-D5)               ;                                
    BXP9 = -ETA * (SU5-SU4);

    % Derivada de betay com relacao a psi                                             
    BYP1 = ETA * SL5;
    BYP2 = 1.d0 - ETA * SU5;
    BYP3 = -BXP2;                                                    
    BYP4 = ETA * SL4;
    BYP5 = -1.d0 + ETA * SU4;                                              
    BYP6 = -BXP5   ;                                                 
    BYP7 = -ETA * (SL5+SL4);                                            
    BYP8 = BXP9 ;                                                    
    BYP9 = -BXP8   ;                                                 

    % Derivada de betax com relacao a eta
    BXE1 = -CL5 * ET2 - PSI * (CL6-CL5)  ;                                  
    BXE2 = D5 * (ET2-PSI)   ;                                            
    BXE3 = -4.d0 + 6.d0 * (PSI+ETA) + SU5E - PSI * SU5    ;                   
    BXE4 = PSI * (CL4+CL6)  ;                                            
    BXE5 = PSI * D4    ;                                                 
    BXE6 = PSI * SU4   ;                                                 
    BXE7 = CL5 * ET2 - PSI * (CL4+CL5);                                      
    BXE8 = D5 * ET2 + PSI * (D4-D5)  ;                                       
    BXE9 = -2.d0 +6.d0 * ETA + SU5E + PSI * (SU4-SU5);

    % Derivada de betay com relacao a eta
    BYE1 = SL5 * (PSI-ET2)   ;                                           
    BYE2 = 1.d0 + SU5E - PSI * SU5   ;                                     
    BYE3 = -BXE2 ;                                                   
    BYE4 = PSI * SL4 ;                                                  
    BYE5 = BXE6   ;                                                  
    BYE6 = -BXE5   ;                                                
    BYE7 = SL5 * ET2 - PSI * (SL4+SL5) ;                                    
    BYE8 = -1.d0 + SU5E + PSI * (SU4-SU5) ;                                    
    BYE9 = -BXE8;

    % Determinacao de [bb]
    B(1,1) = y03 * BXP1 / area02   ;                                        
    B(2,1) = (-x03 * BYP1 + x02 * BYE1) / area02   ;                             
    B(3,1) = (-x03 * BXP1 + x02 * BXE1 + y03 * BYP1) / area02  ;                    
    B(1,2) = y03 * BXP2 / area02 ;                                          
    B(2,2) = (-x03 * BYP2 + x02 * BYE2) / area02   ;                             
    B(3,2) = (-x03 * BXP2 + x02 * BXE2 + y03 * BYP2) / area02  ;                    
    B(1,3) = y03 * BXP3 / area02    ;                                       
    B(2,3) = (-x03 * BYP3 + x02 * BYE3) / area02    ;                            
    B(3,3) = (-x03 * BXP3 + x02 * BXE3 + y03 * BYP3) / area02  ;                    
    B(1,4) = y03 * BXP4 / area02   ;                                        
    B(2,4) = (-x03 * BYP4 + x02 * BYE4) / area02    ;                            
    B(3,4) = (-x03 * BXP4 + x02 * BXE4 + y03 * BYP4) / area02  ;                    
    B(1,5) = y03 * BXP5 / area02   ;                                        
    B(2,5) = (-x03 * BYP5 + x02 * BYE5) / area02     ;                           
    B(3,5) = (-x03 * BXP5 + x02 * BXE5 + y03 * BYP5) / area02  ;                    
    B(1,6) = y03 * BXP6 / area02    ;                                       
    B(2,6) = (-x03 * BYP6 + x02 * BYE6) / area02     ;                           
    B(3,6) = (-x03 * BXP6 + x02 * BXE6 + y03 * BYP6) / area02  ;                    
    B(1,7) = y03 * BXP7 / area02    ;                                       
    B(2,7) = (-x03 * BYP7 + x02 * BYE7) / area02     ;                           
    B(3,7) = (-x03 * BXP7 + x02 * BXE7 + y03 * BYP7) / area02  ;                    
    B(1,8) = y03 * BXP8 / area02    ;                                       
    B(2,8) = (-x03 * BYP8 + x02 * BYE8) / area02    ;                            
    B(3,8) = (-x03 * BXP8 + x02 * BXE8 + y03 * BYP8) / area02  ;                    
    B(1,9) = y03 * BXP9 / area02    ;                                       
    B(2,9) = (-x03 * BYP9 + x02 * BYE9) / area02    ;                            
    B(3,9) = (-x03 * BXP9 + x02 * BXE9 + y03 * BYP9) / area02  ;        

end

function [R] = global_rotation_matrix(q);
    % Eq. (35) of 10.1016/j.cma.2006.10.006
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q0 = sqrt( 1.0 - q(1)^2 - q(2)^2 - q(3)^2 ) ;

    R  = 2*[    [ (q0^2 + q1^2 - 0.5),  (q1*q2 - q0*q3),        (q1*q3 + q0*q2)     ];
                [ (q1*q2 + q0*q3),      (q0^2 + q2^2 - 0.5),    (q2*q3 - q0*q1)     ];
                [ (q1*q3 - q0*q2),      (q2*q3 + q0*q1),        (q0^2 + q3^2  -0.5) ] ];
end

function [v] = rotation_vector(R);
    % Eq. (13) of 10.1016/j.cma.2006.10.006
    v = zeros(3,1);
    v(1) = .5 * ( R(3,2) - R(2,3));
    v(2) = .5 * ( R(1,3) - R(3,1));
    v(3) = .5 * ( R(2,1) - R(1,2));
end

function [Ta] = matrix_Ta(R);
    % Eq. (15) of 10.1016/j.cma.2006.10.006
    Ta = [  [  R(2,2) + R(3,3), -R(1,2),        -R(1,3)];
            [ -R(2,1),           R(1,1)+R(3,3), -R(2,3)];
            [ -R(3,1),          -R(3,2),         R(1,1) + R(2,2) ];
        ];
end

function [Tm] = matrix_Tm(q)
    % Eq. (37) of 10.1016/j.cma.2006.10.006
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q0 = sqrt( 1.0 - q(1)^2 - q(2)^2 - q(3)^2 ) ;

    Tm = 2/q0*[ [ (q0^2 + q1^2),    (q1*q2 - q0*q3),  (q1*q3 + q0*q2)];
                [ (q1*q2 + q0*q3),  (q0^2 + q2^2),    (q2*q3 - q0*q1)];
                [ (q1*q3 - q0*q2),  (q2*q3 + q0*q1),  (q0^2 + q3^2) ] ];
end


function [Khi] = matrix_Khi(R, m);
    % Eq. (21) of 10.1016/j.cma.2006.10.006

    Khi = 0.5 * [   [ ((R(2,3)-R(3,2))*m(1) + R(3,1)*m(2) - R(2,1)*m(3)),   (R(1,1)*m(3) - R(1,3)*m(1)),                            (R(1,2)*m(1) - R(1,1)*m(2)) ];
                    [ (R(2,3)*m(2) - R(2,2)*m(3)),                          ((R(3,1)-R(1,3))*m(2) + R(1,2)*m(3) - R(3,2)*m(1)),     (R(2,2)*m(1) - R(2,1)*m(2)) ];
                    [ (R(3,3)*m(2) - R(3,2)*m(3)),                          (R(3,1)*m(3) - R(3,3)*m(1)),                            ((R(1,2)-R(2,1))*m(3) + R(2,3)*m(1) - R(1,3)*m(2)) ] ];

end


function [G1,G2,G3] = matrix_Gi(a1, a2, a3, r1, r2, r3);
    % Eq. (27) of 10.1016/j.cma.2006.10.006

    a21 = a2-a1;
    a31 = a3-a1;
    a32 = a3-a2;
    a13 = a1-a3;

    v = norm(cross(a21,a31)); 
    c = a1'*r1 + a2'*r2 + a3'*r3;

    G1 = zeros(6,3);
    G1(1,3) = - r1(2)/c;
    G1(2,3) =   r1(1)/c;
    G1(3,1) =   a32(1)/v;
    G1(3,2) =   a32(2)/v;

    G2 = zeros(6,3);
    G2(1,3) = - r2(2)/c;
    G2(2,3) =   r2(1)/c;
    G2(3,1) =   a13(1)/v;
    G2(3,2) =   a13(2)/v;
    
    G3 = zeros(6,3);
    G3(1,3) = - r3(2)/c;
    G3(2,3) =   r3(1)/c;
    G3(3,1) =   a21(1)/v;
    G3(3,2) =   a21(2)/v;
    
end

function [P] = matrix_P(a1, a2, a3, G1, G2, G3);
    % Eq. (26) of 10.1016/j.cma.2006.10.006

    Ai = zeros(5,3);
    Ai(3,1) = 1;
    Ai(4,2) = 1;
    Ai(5,3) = 1;
    
    I = zeros(5,6) ;
    I(1,1) = 1;
    I(2,2) = 1;
    I(3,4) = 1;
    I(4,5) = 1;
    I(5,6) = 1;
    
    P = zeros(15,18) ;
    Ai(1,3) = - a1(2) ;
    Ai(2,3) =   a1(1) ;
    P1 = [ I - Ai*G1', -Ai*G2' , -Ai*G3' ];

    Ai(1,3) = - a2(2);
    Ai(2,3) =   a2(1); 
    P2 = [ -Ai*G1', I - Ai*G2', -Ai*G3' ];

    Ai(1,3) = - a3(2);
    Ai(2,3) =   a3(1);
    P3 = [ -Ai*G1', -Ai*G2', I - Ai*G3' ];

    P = [P1; P2; P3];
end


function [F1, F2] = matrix_F(n);
    % Eq. (30) of 10.1016/j.cma.2006.10.006

    F1 = zeros(15,3);
    F2 = zeros(18,3);

    for i = 1:3
        j1 = 6*(i-1)+1;
        j2 = j1+2;
        j3 = j2+1;
        j4 = 6*i;

        sna = skew(n(j1:j2));
        snb = skew(n(j3:j4));
        
        F1(j1:j2-1,:) = sna(1:2,:);  

        F2(j1:j2,:) = sna;
        F2(j3:j4,:) = snb;

    end
end



function [Kki] = matrix_Kki(q, m);
    % Eq. (41) of 10.1016/j.cma.2006.10.006

    q0s = 1.0 - q(1)^2 - q(2)^2 - q(3)^2;
    if q0s<0
        q0s
        error('ojo')
    end
    q0 = sqrt(q0s);
    A = q(1)*m(1) + q(2)*m(2) + q(3)*m(3);
    
    H = zeros(3, 3);
    H(1, 1) = (q0s + q(1)^2) * A;
    H(2, 2) = (q0s + q(2)^2) * A;
    H(3, 3) = (q0s + q(3)^2) * A;
    H(1, 2) = q0s * (q(1)*m(2) - q(2)*m(1) - q0*m(3)) + q(1)*q(2)*A;
    H(1, 3) = q0s * (q(1)*m(3) - q(3)*m(1) + q0*m(2)) + q(1)*q(3)*A;
    H(2, 1) = q0s * (q(2)*m(1) - q(1)*m(2) + q0*m(3)) + q(2)*q(1)*A;
    H(2, 3) = q0s * (q(2)*m(3) - q(3)*m(2) - q0*m(1)) + q(2)*q(3)*A;
    H(3, 1) = q0s * (q(3)*m(1) - q(1)*m(3) - q0*m(2)) + q(3)*q(1)*A;
    H(3, 2) = q0s * (q(3)*m(2) - q(2)*m(3) + q0*m(1)) + q(3)*q(2)*A;

    Kki = (2/q0^3) * H;
    
end
    