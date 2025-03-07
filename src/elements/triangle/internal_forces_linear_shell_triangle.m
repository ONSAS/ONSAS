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
function [ fs, ks, fintLocCoord ] = internal_forces_linear_shell_triangle(elemCoords, elemDisps, modelName, modelParams, thickness)

    %material and geometric parameters
    E = modelParams(1);  
    nu = modelParams(2) ;
    h = thickness;

    p1 = elemCoords(1:3);
    p2 = elemCoords(4:6);
    p3 = elemCoords(7:9);

    elemDisps_sortT = switchToTypeIndexing( elemDisps ) ;

    [T, x02, x03, y03] = edge_local_axis_shell_triangle(p1,p2,p3);
    Te = blkdiag(T,T,T,T,T,T) ;

    area = x02*y03 / 2;

    % membrane stiffness
    aux1 = h *  E / ( 1 - nu^2) ; 
    aux2 = nu*aux1;
    Dm = [ [aux1, aux2 , 0 ]; [aux2, aux1, 0] ; [0, 0, aux1*(1-nu)/2] ];

    Bm = CST_B(x02, x03, y03);
    Km = area * Bm' * Dm * Bm ;

    % bending stiffness
    aux1 = E * h^3 / (12 * ( 1- nu^2) ); 
    aux2 = nu*aux1;
    Db = [ [aux1, aux2 , 0 ]; [aux2, aux1, 0] ; [0, 0, aux1*(1-nu)/2] ];

    int_point = [ [1./6. 1./6.]; [2./3., 1./6.]; [1./6., 2./3.]];

    im = [1, 2, 7, 8, 13, 14];
    ib = [3,4,5, 9,10,11, 15,16,17];

    fintLocCoord = zeros(1,3);

    dispTe = Te*elemDisps_sortT ;
 
    
    % TODO: temporal FIXX generalize!
    TM = T(1:2,1:2);
    % ------------------------------

    Kb = zeros(9,9);
    wgt = area/3.d0;
    for ipt = 1:3;
        psi = int_point(ipt,1);
        eta = int_point(ipt,2);
        Bb = DKT_B(psi, eta, x02, x03, y03);

        Kb = Kb + wgt * Bb' * Db * Bb;

        fint_ip = Db * Bb * dispTe(ib) ;
        Mmat    = (TM') * [ fint_ip(1) fint_ip(3) ; fint_ip(3) fint_ip(2) ] * TM ;
        fint_ip = [ Mmat(1,1) Mmat(2,2) Mmat(1,2) ] ;
        fintLocCoord = fintLocCoord + fint_ip * wgt/area ;
    end
    
    %assembling the stiffness matrix of the shell element in local coordinates
    Ke = zeros(18,18);

    Ke(im,im) = Km;
    Ke(ib,ib) = Kb;
    
    k_dr = min( min( abs( Kb ) ) ) * 1.e-4;
    Ke(6 , 6) = k_dr;
    Ke(12,12) = k_dr;
    Ke(18,18) = k_dr;

    [Kl_funcao , ~ ] = local_shell_triangle(x02, x03, y03, E, nu, h, zeros(18,1))

    Ke
    Kl_funcao
    norm( Ke - Kl_funcao)
stop


    % calculating the stiffness matrix of the shell element in global coordinates
    Ke = Te' * Ke *Te ;

    %shifting lines and coluns to onsas convention of dofs order
    aux_r = [1,4,2,5,3,6];
    aux_onsas = [aux_r, aux_r+6 , aux_r+12];

    K = Ke(aux_onsas, aux_onsas);

    f = K*elemDisps;

    ks = {K} ; fs = {f};

end

function [T, x02, x03, y03] = edge_local_axis_shell_triangle(p1,p2,p3)
%Calculates the matrix for transformation of basis between global and local axis;
%p1, p2 and p3 are the position vector for the nodes in global coordinates;
%the local x axis is paralel to the side connecting nodes 1 and 2
%the local z axis is normal to the element plane;
%the origin of local axis is located at node 1
    
    p12 = p2 - p1;
    p13 = p3 - p1;
    
    au_zl =  cross(p12,p13);
    u_zl = au_zl / norm(au_zl);
    
    x02 = norm(p12);
    u_xl = p12 / x02;

    u_yl = cross(u_zl, u_xl);

    T = [ u_xl; u_yl; u_zl];

    x03 = dot(u_xl, p13);
    y03 = dot(u_yl, p13);
end


