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

function [ fs, ks, fintLocCoord ] = internal_forces_shell_triangle(elemCoords, elemDisps, modelName, modelParams, thickness)
    
    %calculating local coordinates and coordinates transformation matrix [T]
    p1 = elemCoords(1:3);
    p2 = elemCoords(4:6);
    p3 = elemCoords(7:9);
    
    p12 = p2 - p1;
    p13 = p3 - p1;
    %au_zl =  cross(p12,p13);
    
    %u_zl = au_zl / norm(au_zl);
    %u_xl = cross( [0,1,0] , u_zl)
    %u_yl = cross( u_zl , u_xl)

    %u_xl = p12 / norm(p12);
    %u_zl = au_zl / norm(au_zl);
    %u_yl = cross(u_zl, u_xl);

    %T = [ u_xl; u_yl; u_zl];
    
    T = local_axis_shell_triangle(p1,p2,p3)

    elemCoords_l = [ 0,0,0, (T*p12')' , (T*p13')' ];
    elemCoords_l([3, 6, 9] ) = 0;
    
    %membrane: interal forces and stiffness matrix
    paramOut = [] ;%not used
    planeStateFlag = 1; %plane stress
    dotdotdispsElem = 0;
    density = 0;
    previous_state = cell( 1, 3) ;
    previous_state(:,1) = {zeros( 1, 3 )} ;
    previous_state(:,2) = {zeros( 1, 3 )} ;
    previous_state(:,3) = {0} ;

    elemDisps_m = elemDisps(1:2:end);
    [ fsm, ksm ] = elementTriangSolid( ...
    elemCoords_l, elemDisps_m, modelName, [0,modelParams], paramOut, thickness, planeStateFlag, dotdotdispsElem, density, previous_state );

    %plate: calculation of internal forces and stiffness matrix
    % aux_plate = [2,4,5, 8,10,11, 14, 16,17];
    aux_plate = [5, 2,4, 11,8,10, 17, 14, 16];
    elemDisps_p = elemDisps(aux_plate) ;

    [ fsp, ksp , fintLocCoord_p] = internal_forces_plate_triangle( elemCoords_l, elemDisps_p, modelName, modelParams, thickness );

    %assembling the shell element
    ks = zeros(18,18);
    fint = zeros(18,1);
    
    %assembling the contribution of the plate element
    il = [3,4,5, 9,10,11, 15,16,17];
    fint(il) = fsp{1};
    ks(il,il) = ksp{1};

%   ksp = ksp{1};
%    for ii = 1:9 %adding the plate stiffness coefficients
%        ks(il(ii),il) = ksp(ii,:);
%    end

    %assembling the contribution of the membrane element
    il = [1,2, 7,8, 13,14];
    fsm = fsm{1};
    aux_m = [1,2, 4,5, 7,8];
    fint(il) = fsm(aux_m);
    ksm = ksm{1};
    ks(il,il) = ksm(aux_m, aux_m);

    %artificial stiffness for the drilling dofs
    il = [6, 12, 18];
    ktheta = min( min( abs(ksp{1}) ))*1.e-4;
    ks(il,il) = ktheta * eye(3);

    %transforming to globla coordinates
    Te = blkdiag(T,T,T,T,T,T);

    Ke = Te' * ks *Te;

    %shifting lines and coluns to onsas convention of dofs order
    aux_r = [1,4,2,5,3,6];
    aux_onsas = [aux_r, aux_r+6 , aux_r+12];

    K = Ke(aux_onsas, aux_onsas);

    f = K*elemDisps;

    ks = {K} ; fs = {f};

    fintLocCoord = [ 0 0 0];

function [T] = local_axis_shell_triangle(p1,p2,p3);

    p12 = p2 - p1;
    p13 = p3 - p1;
    au_zl =  cross(p12,p13);
    u_zl = au_zl / norm(au_zl);

    if abs(u_zl[2]) < 1 - 1.e-3 ;
        u_xl = cross( [0,1,0] , u_zl)
    else
        u_xl = [1,0,0]
    else;
        u_yl = cross( u_zl , u_xl)

    T = [ u_xl; u_yl; u_zl];