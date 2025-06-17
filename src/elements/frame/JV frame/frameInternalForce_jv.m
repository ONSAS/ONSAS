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


function  [fs, ks, stress, rotData, localInternalForces] = frameInternalForce_jv(elemCoords, elemCrossSecParams, elemConstitutiveParams, Ue, rotMat)

    % element coordinates
    xs = elemCoords(:);
    r1_g = xs(1:3);
    r2_g = xs(4:6);


    % ----- material and geometric params ------
    E   = elemConstitutiveParams(2);
    nu  = elemConstitutiveParams(3);
    G   = E / (2 * (1 + nu));
    % ----- extract cross section properties ---
    [Area, J, Iyy, Izz, ~] = crossSectionProps (elemCrossSecParams, 0); % select a ficticious elemrho
    % ------------------------------------------

    % Disps and spatial rotations in global reference frame
    Ug = switchToTypeIndexing(Ue);
    u1_g    = Ug(1:3);
    t1_g    = Ug(4:6);
    u2_g    = Ug(7:9);
    t2_g    = Ug(10:12);
    % -------------------------------

    % Updated position vector in global reference frame
    p1_g = r1_g + u1_g;
    p2_g = r2_g + u2_g;

    % Global rotation matrix
    % Rg1 = expm(skew(t1_g));
    % Rg2 = expm(skew(t2_g));
    Rg1 = rotMat{1};
    Rg2 = rotMat{2};

    % Element undeformed and deformed lengths - R0 rotation matrix
    
    [x21, d21, l, l0] = corotLenCoords(xs, Ug);
    Ro = beamRefConfRotMat(x21);
    To = Ro' ;

    % Rotation matrix Rr
    Rr = RotationMatrixRr(x21, d21, l, Rg1, Rg2, Ro) ; 
    Tr = Rr' ; 

    % % Co rot mats
    % [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices_jv(Ug, elemCoords, rotMat);

    % Choice of element origin
    rc_g = r1_g ;
    pog  = p1_g ;
    
    % nodal displacements in local reference frame in deformed configuration
    r1_o = To * (r1_g - rc_g);
    r2_o = To * (r2_g - rc_g);

    u1_def = Rr' * (p1_g - pog) - r1_o ;
    u2_def = Rr' * (p2_g - pog) - r2_o ;

    a1_def = u1_def + r1_o;
    a2_def = u2_def + r2_o;

    % Rotation matrix from local reference frame to nodal reference frame in deformed configuration
    R1_def = Rr' * Rg1 * Ro;
    R2_def = Rr' * Rg2 * Ro;

    v1_def = rotationVector(R1_def, 0);
    v2_def = rotationVector(R2_def, 0);
    
    
    u_def = u2_def(1) - u1_def(1) ; 
    % --- local force vector and tangent stiffness matrix ---
    [fl_aux, kl_aux, strain, stress] = beamLocalStaticForces (u_def, v1_def, v2_def, l0, E, G, Area, Iyy, Izz, J);

    % dofs_aux = [ 1, 4, 5, 6, 7, ];

    % fl_aux
    % fl = zeros(12,1) ;
    % fl(1)       = -fl_aux(1) ;
    % fl(7)       =  fl_aux(1) ;
    % fl(4:6)     =  fl_aux(2:4) ;
    % fl(10:12)   =  fl_aux(5:7) ;
    fl = fl_aux ;
    Kl = kl_aux ;


    % auxiliary q base
    q1g = Rg1 * Ro * [0 1 0]';
    q2g = Rg2 * Ro * [0 1 0]';
    qg  = (q1g + q2g) / 2;

    
    
    % Change of variables    
    %   -- term to equal virtual work of moments in local an global --
    % transformation to the new local coordinates
    invTs_1 = invTs_jv(v1_def);
    invTs_2 = invTs_jv(v2_def);

    % Ba = zeros(12);
    % Ba(1:3, 1:3)      = eye(3) ;
    % Ba(4:6, 4:6)      = invTs_1;
    % Ba(10:12, 10:12)  = invTs_2;
    
    Ba = zeros(7,7) ;
    Ba = blkdiag(1, invTs_1, invTs_2) ;

    % Kl = zeros(12,12) ;

    Kh1 = dinvTs_jv(v1_def, fl(2:4))    * invTs_1 ;
    Kh2 = dinvTs_jv(v2_def, fl(5:7))  * invTs_2 ;
    Kh = blkdiag(0, Kh1, Kh2) ;

    fa = Ba'*fl ;
    Ka = Ba'*Kl*Ba + Kh ;

    Gt = matrixGt_jv(q1g, q2g, Rr, l) ;
    P = matrixP_jv(Gt) ;

    EE = blkdiag(Rr, Rr, Rr, Rr) ;

    r = [-Rr(:,1)' zeros(1,3)  Rr(:,1)' zeros(1,3)];
    Bg = [ r ; P*EE'] ;

    fg_jv = Bg'*fa ;


    % stop
    
    


    % mauricio
    dg = switchToTypeIndexing(Ue);
    % global thetas
    tg1 = dg(4:6);
    tg2 = dg(10:12);
    % length and coords of the element
    [x21, d21, l, l0] = corotLenCoords(xs, dg);
    % --- auxiliary vector and matrices  ---
    % aux zero matrices
    [I3, O3, O1, II] = corotZeros();


    [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices(Ue, elemCoords);
    % ===========================================================================
    % axial displacement
    u   = l - l0;
    tl1 = logar(Rroof1);
    tl2 = logar(Rroof2);

    locDisp = [u tl1' tl2'];
    % -------------------------------

    [nu, nu11, nu12, nu21, nu22, e1, e2, e3, r, Gaux, P, EE] = corotVecMatAuxStatic( ...
                                                                                    R0, Rr, Rg1, Rg2, l, II, O3, O1);
    % -------------------------------
    [fl, kl, strain, stress] = beamLocalStaticForces (u, tl1, tl2, l0, E, G, Area, Iyy, Izz, J);
    % -- term to equal virtual work of moments in local an global --
    % transformation to the new local coordinates
    De1 = invTs(tl1);
    De2 = invTs(tl2);

    % matrix for transformation between global and relative rotations/moments
    H  = [1   O1   O1; ...
            O1' De1   O3; ...
            O1'  O3  De2];

    localInternalForces = [fl(1) fl(3) fl(4)];
    fe = H' * fl;

    Dh1 = dinvTs(tl1, fl(2:4)) * De1;
    Dh2 = dinvTs(tl2, fl(5:7)) * De2;

    Kh = [0   O1   O1
            O1' Dh1   O3
            O1'  O3  Dh2];

    ke = H' * kl * H + Kh;

    %  -------------transformation to the global coordinates-------
    B = [r'
        -nu / l * e3' (1 - nu12 / 2) * e1' + nu11 / 2 * e2'  nu / l * e3' 1 / 2 * (-nu22 * e1' + nu21 * e2')
        -e3' / l e2' e3' / l 0 0 0
        e2' / l e3' -e2' / l 0 0 0
        -nu / l * e3' 1 / 2 * (-nu12 * e1' + nu11 * e2')  nu / l * e3' (1 - nu22 / 2) * e1' + nu21 / 2 * e2'
        -e3' / l 0 0 0 e3' / l e2'
        e2' / l 0 0 0 -e2' / l e3'];

    fg = B' * fe;

    A  = (I3 - e1 * e1') / l;

    Dr = [A  O3 -A  O3
        O3 O3  O3 O3
        -A  O3  A  O3
        O3 O3  O3 O3];

    F = P' * fe(2:7);

    sF = [skew(F(1:3))
        skew(F(4:6))
        skew(F(7:9))
        skew(F(10:12))];

    nab = [0
            (nu * (fe(2) + fe(5)) + fe(3) + fe(6)) / l
            (fe(4) + fe(7)) / l];

    Kg = B' * ke * B + Dr * fe(1) - EE * sF * Gaux' * EE' + EE * Gaux * nab * r';

    Dg1 = funTs(tg1);
    Dg2 = funTs(tg2);

    q = [fg(1:3)
        Dg1' * fg(4:6)
        fg(7:9)
        Dg2' * fg(10:12)];

    Dk1 = dTs(tg1, fg(4:6));
    Dk2 = dTs(tg2, fg(10:12));

    H = [I3 O3  O3 O3
        O3 Dg1 O3 O3
        O3 O3  I3 O3
        O3 O3  O3 Dg2];

    Kt = H' * Kg * H;

    Kt(4:6, 4:6) = Kt(4:6, 4:6) + Dk1;
    Kt(10:12, 10:12) = Kt(10:12, 10:12) + Dk2;

    % make the internal tangent matrix symmetric to improve convergence order
    Kt = (Kt + Kt') / 2;

    % --- Write it back to ONSAS nomenclature [fx mx,....]  ---
    Finte = switchToNodalIndexing(q);

    KTe = zeros(size(Kt));
    KTe = switchToNodalIndexing(Kt);

    % -------------------------------------------------------


    % u
    % [u1_def u2_def]
    % [ v1_def tl1 v2_def tl2 ]
    % stop

    % fprintf('Ro \n')
    % R0_m ./Ro
    % fprintf('Rr \n')
    % Rr_m ./Rr
    % fprintf('Rg \n')
    % Rg1_m
    % Rg1
    % fprintf('------ \n')
    % Rg2_m 
    % Rg2
    % fprintf('Rroof \n')
    % R1_def
    % Rroof1_m
    % fprintf('------ \n')
    % R2_def
    % Rroof2_m

    % [fa fe] % pequeña diferencia en fuerza interna axial
    % fa./fe
    % Ka./ke

    [ fg_jv fg ]




function [Rr] = RotationMatrixRr(x21, d21, l, Rg1, Rg2, R0)

    % deformed x axis
    e1 = (x21 + d21) / l;
    % auxiliary vectors q
    q1 = Rg1 * R0 * [0 1 0]';
    q2 = Rg2 * R0 * [0 1 0]';
    q  = (q1 + q2) / 2;

    % deformed z local axis
    e3 = cross(e1, q);
    e3 = e3 / norm(e3); % normalization

    % deformed y local axis
    e2 = cross (e3, e1);

    % rotation matrix
    Rr = [e1 e2 e3];
    



