% Copyright 2023, ONSAS Authors (see documentation)
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
% This function returns Drag, lift and pitch moment forces of an inner point in a frame element.
% cross section considered in global coordinates.
% the computation is done according to (M.C. Vanzilli, J.M. Perez Zerpa, 2022)

function integFluidForce = integFluidForce( x, ddotg, udotFlowElem                                    ,...
                                            l0, tl1, tl2, Rr                                          ,...
                                            vecChordUndef, dimCharacteristic, I3, O3, P, G, EE, L2, L3,...
                                            aeroCoefs, densityFluid, viscosityFluid                   ,...
                                            VIVBool, q, p, constantLiftDir, uniformUdot, tlift1, tlift2, ...
                                            fluidFlowBool, ILVIVBool, uBEMbool )
  
%% -----------------------------------------------------------------------
% Bernoulli weight function
[N1, N2, N3, N4, N5, N6, N7, N8] = bernoulliInterpolWeights(x, l0) ;
% Auxiliary matrices
[P1, P2, N, N1, N2] = corotVecMatAuxDyn( N1, N2, N3, N4, N5, N6, N7, N8, tl1, tl2, G, I3, O3, P ) ;
% -----------------------------------------------------------------------

%% -----------------------------------------------------------------------
% Local displacements of a generic cross section
ul = P1 * [ tl1; tl2 ]                   ; % Eq.(38)  T-N Le J.-M. Battini et al 2014

% Auxiliary matrices H
H1  = N + P1 * P - 1 * skew( ul ) * G' ;
H2  = P2 * P + G'; %Ec 72 se puede usar para comprobar con ec A.10

% angular local rotation
thethaRoof  = P2 * [tl1 ; tl2] ; % Eq. 39 Le, Battini 2014

% local Rroof rotation matrix is
Rroofx      = expon( thethaRoof ) ;
%---------------------------------------------------------

%% ----------------------------------------------------------------
% Kinematic velocities for the generic cross section
% cross section centroid rigid velocity in global coordinates:
% if uniform
% -----------------------------------------------------------------
if ~isempty( VIVBool ) && ~isempty( constantLiftDir ) && ~isempty( uniformUdot )
    if uniformUdot
      udotG = (ddotg(1:3) + ddotg(7:9))/2; % nodal velocities averaged
    else
      udotG = Rr * H1 * EE' * ddotg ; % Eq.(61)  T-N Le J.-M. Battini et al 2014
    end
else
    udotG = Rr * H1 * EE' * ddotg ; % Eq.(61)  T-N Le J.-M. Battini et al 2014
end

%% -----------------------------------------------------------------
% Compute relative velocity of interpolated cross section
% cross section absolute fluid flow velocity in global coordinates interpolated with linear shape functions:
udotFlowG = udotFlowElem(1:3) * N1 + udotFlowElem(4:6) * N2 ;
% Relative, perpendicular and projected  flow velocity of the cross section to compute drag lift and moment:
[VpiRelGflow, VpiRelGperpflow, VrelGflow] = computeVpiRels( udotFlowG, [0 0 0]', Rroofx, Rr, L2, L3 ) ;  
% -----------------------------------------------------------------
% uBEM compute relative local velocity
if ~isempty( uBEMbool ) && uBEMbool

    % load induced velocity node 1 and node 2
    global inducedVelNod1    ;   global inducedVelNod2;

    % Geom parameteres
    global elemTwist; 

    % Interpolated chord vector in evaluated gaussian section interpolated with linear shape functions
    vecChordUndef     = vecChordUndef(1:3) * N1 + vecChordUndef(4:6) * N2 ;    
    
    % Interpolated twist vector in evaluated gaussian section interpolated with linear shape functions
    twistVector       = elemTwist(1:3)' * N1 + elemTwist(4:6)' * N2 ;
    
    % Interpolated twist vector in evaluated gaussian section interpolated with linear shape functions
    dimCharacteristic = norm(vecChordUndef) ;    
    
    % cross section induced velocity flow velocity in global coordinates interpolated with linear shape functions
    inducedVelGn1     = inducedVelNod1 * N1  + inducedVelNod2 * N2 ;     % induced wake velocity of previous time step    
    
    % Kinematic velocities for the generic cross section
    % cross section centroid rigid velocity in global coordinates
    udotG = Rr * H1 * EE' * ddotg ; % Eq.(61)  T-N Le J.-M. Battini et al 2014
    
    [VpiRelG, VpiRelGperp, VrelG] = uBEMcomputeVpiRels( udotFlowG, udotG,...
                                       inducedVelGn1, Rroofx, Rr, L2, L3 ) ;

elseif ~isempty( fluidFlowBool ) && fluidFlowBool % Leclercq validation
    % VpiRelG along y
    [VpiRelG, VpiRelGperp, VrelG] = computeVpiRels( udotFlowG, [0 udotG(2) 0]', Rroofx, Rr, L2, L3 )  ;
else
    [VpiRelG, VpiRelGperp, VrelG] = computeVpiRels( udotFlowG, udotG, Rroofx, Rr, L2, L3 )  ;
end
  
%-----------------------------------------------------------------

% ------------ Compute relative incidence angle  ------------
% the chord vector orientation in the deformed coordinates to compute incidence flow angle is:
if ~isempty( uBEMbool ) && uBEMbool
    %foilTwist   = dot( (Rroofx'*Rr'*deg2rad( twistVector ) )', [1,0,0] ) ;
    %chordVec    = expon( [foilTwist 0 0] )*( L2*(Rroofx'*Rr')*vecChordUndef' ) ;
    vecChordUndef    = expon( deg2rad(twistVector) )*( (Rroofx'*Rr')*vecChordUndef ) ;
    tch = ( vecChordUndef / norm( vecChordUndef )) ;
else
    tch = (vecChordUndef / norm( vecChordUndef )) ;
end

%disp('vpi')
%norm( VpiRelG)

% Calculate relative incidence angle in the deformed configuration
if( norm( VpiRelG ) == 0 )
    td = tch ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
elseif ~isempty( uniformUdot ) && uniformUdot % Verification for small disp
    td = VpiRelGflow/ norm( VpiRelGflow ) ; % constant along x
    tlconst = VpiRelGperpflow/ norm( VpiRelGperpflow) ; % constant along y
else % the drag direction at a generic cross section in deformed coordinates is:
    td = VpiRelG / norm( VpiRelG ) ;
end

if isnan(  norm( VpiRelG)  ),  stop, end

if ~isempty( uBEMbool ) && uBEMbool
    cosBeta  = dot( tch, td ) / ( norm(td) * norm(tch) ) ;
    sinBeta  = dot( cross(tch, td), [1 0 0] ) / ( norm( td ) * norm( tch ) ) ;
    betaRelG = sign( sinBeta ) * acos( cosBeta ) ;
else
    cosBeta  = dot( tch, td ) / ( norm(td) * norm(tch) ) ;
    sinBeta  = dot( cross(td,tch), [1 0 0] ) / ( norm( td ) * norm( tch ) ) ;
    betaRelG = sign( sinBeta ) * acos( cosBeta ) ;
end
% ------------------------------------------------------------------

%% -----------------------------------------------------------------
% Compute the lift, drag and pitch coefs corresponding to the uBEM model
if ~isempty( uBEMbool ) && uBEMbool    
    % Interpolated Cl, Cd and Cm in evaluated gaussian section interpolated with linear shape functions   
    global cdstat1; global clstat1; global cmstat1;
    global cdstat2; global clstat2; global cmstat2;
    
    c_d  = cdstat1 * N1 + cdstat2 * N2 ;
    c_l  = clstat1 * N1 + clstat2 * N2 ;
    c_m  = cmstat1 * N1 + cmstat2 * N2 ;
else
  %-----------------------------------------------------------------
    
    userDragCoef   = aeroCoefs{1} ;
    userLiftCoef   = aeroCoefs{2} ;
    userMomentCoef = aeroCoefs{3} ;
    
    % -----------------------------------------------------------------
    % Computation of Renynolds number
    
    % Computation of Reynolds number
    
    Re = norm(udotFlowG) * dimCharacteristic / viscosityFluid ;
    
    % -----------------------------------------
    % ------------ Read Cd, Cl, Cm  ------------
    
    
    % Check fluid coefficients existence and the load it values if not set 0:
    if ~isempty( userDragCoef )
        c_d = feval( userDragCoef, betaRelG, Re  ) ;
        c_d_il = 0.1; % IL VIV
    else
        c_d = 0 ;
    end
    if ~isempty( userLiftCoef )
        c_l = feval( userLiftCoef, betaRelG, Re  ) ;
    else
        ~isempty( VIVBool ) && VIVBool && error('The lift CL0 coef function must be defined for VIVBool problems ') ;
        c_l = 0 ;
    end
    if ~isempty( userMomentCoef )
        c_m = feval( userMomentCoef, betaRelG, Re ) ;
    else
        c_m = 0 ;
    end
end
%--------------------------------------------------------------------

%% -------------------------------------------------------------------
% ------------ Compute drag, lift and pitch moment forces  ------------
  % The cross section fluid forces in deformed coordinates is:
  % drag cross section force vector in deformed coordinates

if ~isempty( uniformUdot ) && uniformUdot
    
    fdl = 1/2 * densityFluid * c_d * dimCharacteristic * norm( VpiRelG)^2 * td    ;

else
    
    fdl = 1/2 * densityFluid * c_d * dimCharacteristic * norm( VpiRelG) * VpiRelG     ;

    if isnan(norm(fdl)),stop,end
end
if ~isempty(  ILVIVBool ) && ILVIVBool
    fdl_il =  1/2 * densityFluid * c_d_il * p/2 * dimCharacteristic * norm( VpiRelGflow )^2 * td   ;  % U^2 along VpiRelG
else
    fdl_il =  [0 0 0]' ;
end
% lift cross section force vector in deformed coordinates
if ~isempty( VIVBool ) && ~isempty( constantLiftDir ) && ~isempty( uniformUdot )
    if uniformUdot 
        fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelG )^2 * tlconst;
    elseif constantLiftDir % lift direction is constant
        %prom the lift direction in global coordinates
        tlift = (tlift1 + tlift2) / 2 ;
        % transform the lift direction into deformed coordinates to re use the Eq in line 330
        tlift_defCoords = Rroofx' * Rr' * tlift / norm( Rroofx' * Rr' * tlift  ) ;
        % compute the lift force in deformed coordinates
        if ~isempty( fluidFlowBool ) && fluidFlowBool
            fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelGflow )^2 * tlift_defCoords ;
        else
            fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelG )^2 * tlift_defCoords ;
        end
    else % lift direction is variable
        if ~isempty( fluidFlowBool ) && fluidFlowBool 
            fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelGflow ) * VpiRelGperp ;
%        elseif ~isempty( ILVIVBool ) && ILVIVBool % Trim validation
%          fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelG ) * VpiRelGperp ;
        else % Trim validation
            fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelG ) * VpiRelGperp ; %note that if there is VIV effect q is 2
        end
    end

else % no WOM and a variable lift direction

    fll =  1/2 * densityFluid * c_l * q / 2 * dimCharacteristic * norm( VpiRelG ) * VpiRelGperp ; %note that if there is VIV effect q is 2

end


% drag + lift cross section force vector in deformed coordinates
fal =  fdl + fll + fdl_il;
% torsional moment fluid load in deformed coordinates
ma =  1/2 * densityFluid * c_m * (VpiRelG)' * VpiRelG * dimCharacteristic * ( [1 0 0]' ) ;
% --------------------------------------------------------------------

%% --------------------------------------------------------------------
% Compute the element fluid load forces vector in global coordinates
% compute the integral term of the current cross section in rigid coordinates
integralTermAeroForceRigid  =   H1' * Rroofx * fal + H2' * Rroofx * ma ;
% rotate to global coordinates with EE matrix for rigid configuration formulation
integFluidForce  =  EE *( integralTermAeroForceRigid ) ; %Rotate from rigid to global coordinates
%-----------------------------------------------------------------
end
