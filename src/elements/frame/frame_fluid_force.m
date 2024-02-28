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
% This function computes fluid forces as proposed in https://arxiv.org/abs/2204.10545
function [fHydroElem, tMatHydroElemU] = frame_fluid_force( elemCoords           , ...
                                     elemCrossSecParams                         , ...
                                     Ue, Udote, Udotdote                        , ...
                                     aeroCoefs, chordVector, aeroNumericalParams ,...
                                    analysisSettings, nextTime, currElem         ,...
                                    computeAeroStiffnessMatrix)

  % Check all required parameters are defined
  assert( ~isempty( analysisSettings.fluidProps), ' empty analysisSettings.fluidProps.' )


  %% Declare booleans for VIV model

  % Declare booleans for VIV model

  global VIVBool
  global ILVIVBool
  global constantLiftDir
  global uniformUdot
  global fluidFlowBool
  global uBEMbool; global DWMbool; 

  AMBool = analysisSettings.addedMassBool ;


  % Implementation Booleans for internal test, baseBool changes the local angles computation
  baseBool = false ;
  
  %% ------------------------------------------------------
  % Define chord and aero coef parameters of the element
  % -------------------------------------------------------
  if uBEMbool ;
      % Boolean to run uBEM theory for HAWT model
      global uBEMdataCoords ;
      global polars ; aeroData   = polars ; 
      % blade data including, radius, twist, chord and thick
      [ aoast, liftCoef, dragCoef, momCoef ] = uBEMAeroProps( aeroData ) ;
      aeroCoefs{1} = aoast    ;
      aeroCoefs{2} = liftCoef ;
      aeroCoefs{3} = dragCoef ;
      aeroCoefs{4} = momCoef  ;
  end

  % extract fluid properties
  densityFluid   = analysisSettings.fluidProps{1,1} ;
  viscosityFluid = analysisSettings.fluidProps{2,1} ;
  userFlowVel    = analysisSettings.fluidProps{3,1} ;

  % check user Flow Vel is not empty
  assert( ~isempty( userFlowVel ), 'empty user windvel' )

  % extract nonLinearity in aero force boolean
  numGaussPoints         = aeroNumericalParams{1};
  geometricNonLinearAero = aeroNumericalParams{3} ;

  % Boolean to compute the fluid force with ut = 0, this should be used if the fluid loads are computed in the reference
  % configuration and nonlinear effects are not considered.
  if ~geometricNonLinearAero
    Ue = zeros(12,1) ;
  end

  % fluid velocity at the nodes of the element evaluated in deformed configuration (spatial points):
  udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6) , nextTime ) ;
  udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime ) ;
  % compact them into a single vector for the element
  udotFlowElem  = [udotFlowNode1; udotFlowNode2] ;

  % Elem reference coordinates:
  xs = elemCoords(:) ;
  
  % Load element properties to fluid loads
  % length of the chord vector
  if ~isempty( uBEMbool ) && uBEMbool
      global elemIDsection; global elemRadius; global elemBeta;
      clear chordVector
      % Frame geom properties
      [elemIDsection, elemRadius, chordVector, elemBeta, nodIdx1, nodIdx2] = uBEMframeProps(uBEMdataCoords, xs);
      dimCharacteristic = norm( ( chordVector(:,1) + chordVector(:,2) )/2 ) ;
  else
      dimCharacteristic = norm( chordVector ) ;
  end

  % compute corotational rotation matrices
  [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices( Ue, elemCoords ) ;

  % --- global kinematics ---
  % permut indexes according to Battini's nomenclature
  dg = switchToBattiniNom( Ue ) ;
  ddotg    = switchToBattiniNom( Udote ) ;
  ddotdotg = switchToBattiniNom( Udotdote ) ;
  % -------------------------------

  % length and coords of the element
  [x21, d21, l, l0] = corotLenCoords(xs ,dg) ;

  % --- auxiliary vector and matrices  ---
  % aux zero matrices
  [I3, O3, O1, II] = corotZeros() ;

  % auxiliary q base
  q1g = Rg1 * R0 * [0 1 0]' ;
  q2g = Rg2 * R0 * [0 1 0]' ;
  qg  = ( q1g + q2g ) / 2   ;

  [nu, nu11, nu12, nu21, nu22, e1, e2, e3, r, G, P, EE ] = corotVecMatAuxStatic(...
                                                                  R0, Rr, Rg1, Rg2, l, II, O3, O1);
  % -------------------------------
  
  % local rotations
  if ~baseBool ;
    Rroof1 = Rr' * Rg1 * R0 ;% Eq(30) J.-M. Battini 2002
    Rroof2 = Rr' * Rg2 * R0 ;% Eq(30) J.-M. Battini 2002
  else
    Rroof1 = Rr' * R0 * Rg1 ;
    Rroof2 = Rr' * R0 * Rg2 ;
  end

  tl1 = logar( Rroof1 ) ; % Eq(31) J.-M. Battini 2002
  tl2 = logar( Rroof2 ) ;% Eq(31) J.-M. Battini 2002
  
  % auxiliary matrix created for uBEM method to compute twist angle
  L1 = [ 1 0 0   ;
         0 0 0   ;
         0 0 0]  ;

  % auxiliary matrix created to project transversal velocity
  L2 = [ 0 0 0   ;
         0 1 0   ;
         0 0 1 ] ;

  % 90 degrees rotation matrix
  L3 = expon( [pi/2 0 0] ) ;

  % --------------------------------------------------------

  % Extract points and weights for numGausspoints selected
  [xIntPoints, wIntPoints] = gaussPointsAndWeights( numGaussPoints ) ;

  % Fluid velocity definition
  % fluid velocity at the nodes of the element evaluated in deformed configuration (spatial points):
  % check user Flow Vel is not empty
  assert( ~isempty( userFlowVel ), 'empty user windvel' )
  if ~isempty( uBEMbool ) && uBEMbool
      % induced velocity vector of node 1 and node 2 in global corotational system
      % Init static wake velocity in HAWT global coordiantes
      global wWake; global iterCounter; global nonWakebool;
      
      % compute section wind velocity from global HAWT system to
      % global corotational system
      global timeIdx; timeIdx = round(nextTime/analysisSettings.deltaT);
    
      % Init induced vel components at first enter to iteration loop
      global inducedVel; global inducedVelNod1; global inducedVelNod2;
      if nonWakebool && iterCounter == 0
          inducedVel  = wWake{1} ;
          iterCounter = 1 ;
      elseif ~nonWakebool && iterCounter == 0
          inducedVel  = wWake{timeIdx} ;
          iterCounter = 1 ;
      end
      
      inducedVelNod1 = inducedVel(nodIdx1,:)';
      inducedVelNod2 = inducedVel(nodIdx2,:)';

      %% -------------------------------------------------------------------------------- 
      % Node 1
      % Velocity params
      
      % Node 1 frame velocity
      udotFrame1 = Udote( 1:2:6 ) ;
      % Node 1 flow velocity
      udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6), nextTime )  ;
      % Node 1 induced velocity of previous time step

      % Compute relative velocity of node 1 
      [VpiRelNode1, VpiRelperpNode1, VrelGnode1] = uBEMcomputeVpiRels( udotFlowNode1, ...
                                    udotFrame1, inducedVelNod1, Rroof1, Rr, L2, L3 ) ;
      
      % Compute angle of attack of node 1
      chRef1  = chordVector(:,1)  ; 
      tch1    = chRef1 / norm( chRef1 ) ;
      if( norm( VpiRelNode1 ) == 0 )
          td1 = tchNod1 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
      else % the drag direction at a generic cross section in deformed coordinates is:
          td1 = VpiRelNode1 / norm( VpiRelNode1 ) ;
      end
      cosBetanod1 = dot( td1, tch1 ) / ( norm(td1) * norm(tch1) ) ;
      sinBetanod1 = dot( cross(td1, tch1), [1 0 0] ) / ( norm( td1 ) * norm( tch1 ) ) ;
      betaRelnod1 = sign( sinBetanod1 ) * acos( cosBetanod1 ); 
      
      % Compute node 1 aero params
      global clstat1; global cdstat1; global cmstat1;
      [clstat1, cdstat1, cmstat1] = uBEMinterpAeroParams(aeroCoefs, elemIDsection(1), betaRelnod1);
      % ----------------------------------------------------------------------------------
      % Node 2
      % Velocity params
      
      % Node 2 frame velocity
      udotFrame2 = Udote( 7:2:12 ) ;
      % Node 2 flow velocity
      udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime )  ;
      % Node 2 induced velocity of previous time step

      % Compute relative velocity of node 2 
      [VpiRelNode2, VpiRelperpNode2, VrelGnode2] = uBEMcomputeVpiRels( udotFlowNode2, ...
                                    udotFrame2, inducedVelNod2, Rroof2, Rr, L2, L3 );
      
      % Compute angle of attack of node 2
      chRef2  = chordVector(:,2); 
      tch2    = chRef2 / norm( chRef2 ) ;
      if( norm( VpiRelNode2 ) == 0 )
          td2 = tchNod2 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
      else % the drag direction at a generic cross section in deformed coordinates is:
          td2 = VpiRelNode2 / norm( VpiRelNode2 ) ;
      end
      cosBetanod2 = dot( td2, tch2 ) / ( norm(td2) * norm(tch2) ) ;
      sinBetanod2 = dot( cross(td2, tch2), [1 0 0] ) / ( norm( td2 ) * norm( tch2 ) ) ;
      betaRelnod2 = sign( sinBetanod2 ) * acos( cosBetanod2 );  
    
      % Compute node 1 aero params
      global clstat2; global cdstat2; global cmstat2;
      [clstat2, cdstat2, cmstat2] = uBEMinterpAeroParams(aeroCoefs, elemIDsection(2), betaRelnod2);  
%       check  = [ udotFrame1'; udotFrame2'; VpiRelNode1'; VpiRelNode2' ]
%       param1 = [ nodIdx1, clstat1, betaRelnod1 ; nodIdx2, clstat2, betaRelnod2]
  else
      udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6), nextTime ) ;
      udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime ) ;
  end
  % compact them into a single vector for the element

  %% ----------------------------------------------------------
  % WOM or uBEM computation call for cases with VIVbool or uBEMbool equal to true
  % -----------------------------------------------------------
  if ~isempty( VIVBool )
    %% ----------------------------------------------------------------
    % WOM computational
    if VIVBool && ( norm(udotFlowNode1)*norm(udotFlowNode2) ) > 0 % VIV with non null flow
      % extract accelerations and velocities (global coordinates) of the nodes
      % node 1
      udotFrame1 = Udote( 1:2:6 )   ;   udotdotFrame1 = Udotdote( 1:2:6 )   ;
      % node 2
      udotFrame2 = Udote( 7:2:end ) ;   udotdotFrame2 = Udotdote( 7:2:end ) ;

      % projected velocities at nodes 1 and 2 in deformed coordinates
      % node 1
      [VpiRel1_defCords, VpiRelPerp1_defCords, Vrel1_glob] = computeVpiRels( udotFlowNode1, udotFrame1,...
                                                                             Rroof1, Rr, L2, L3 ) ;
      % node 2
      [VpiRel2_defCords, VpiRelPerp2_defCords, Vrel2_glob] = computeVpiRels( udotFlowNode2, udotFrame2,...
                                                                             Rroof2, Rr, L2, L3 ) ;

      % transform them into global coordinates
      %node 1
      VpiRel1     =  Rr * Rroof1 * VpiRel1_defCords      ;
      VpiRelPerp1 =  Rr * Rroof1 * VpiRelPerp1_defCords  ;
      %node 2
      VpiRel2     =  Rr * Rroof2 * VpiRel2_defCords      ;
      VpiRelPerp2 =  Rr * Rroof2 * VpiRelPerp2_defCords  ;

      % directions of lift forces

      if  ~isempty( constantLiftDir ) && constantLiftDir 
        % Compute fluid mean velocity at the nodes on the initial configuration
        % find the first veloctiy direction unitll is not null
        t0 = 0; timeStepNotNullVel = 0;
        udotFlowNode10 = feval( userFlowVel, elemCoords(1:3)', t0 ) ;
        udotFlowNode20 = feval( userFlowVel, elemCoords(4:6)', t0 ) ;
        while norm( udotFlowNode10 ) == 0 && norm( udotFlowNode20 ) == 0
          timeStepNotNullVel = timeStepNotNullVel + 1;
          t0 = timeStepNotNullVel*analysisSettings.deltaT ;
          udotFlowNode10 = feval( userFlowVel, elemCoords(1:3)', t0 ) ;
          udotFlowNode20 = feval( userFlowVel, elemCoords(4:6)', t0 ) ;
        end
        % Compute the direction of the axial vector on the initial configuration in global coordinates
        e1 = R0 * [1 0 0]';
        % Transform axial vector in the initial configuration to global cooridantes
        tlift1 = cross(e1, udotFlowNode10 ) / norm( cross(e1,udotFlowNode10) ) ;
        tlift2 = cross(e1, udotFlowNode20 ) / norm( cross(e1,udotFlowNode20) ) ;
      else 
        tlift1 = VpiRelPerp1 / norm( VpiRelPerp1 ) ;
        tlift2 = VpiRelPerp2 / norm( VpiRelPerp2 ) ;
      end


      % computes van der pol solution for current element
      % node 1
      [VpiRel1_defCords, VpiRelPerp1_defCords, Vrel1_glob] = computeVpiRels( udotFlowNode1, [0 0 0]',...
                                                                             Rroof1, Rr, L2, L3 ) ;
      % node 2
      [VpiRel2_defCords, VpiRelPerp2_defCords, Vrel2_glob] = computeVpiRels( udotFlowNode2, [0 0 0]',...
                                                                             Rroof2, Rr, L2, L3 ) ;
      VprojRel1     =  Rr * Rroof1 * VpiRel1_defCords      ; % Do NOT depend on section velocity
      VprojRel2     =  Rr * Rroof2 * VpiRel2_defCords      ;% Do NOT depend on section velocity
      tlflow1 = VpiRelPerp1_defCords/norm(VpiRelPerp1_defCords);
      tlflow2 = VpiRelPerp2_defCords/norm(VpiRelPerp2_defCords);
      q = WOMV4( VprojRel1, VprojRel2,  udotdotFrame1,udotdotFrame2,...
                 tlflow1, tlflow2, dimCharacteristic, nextTime, analysisSettings.deltaT, currElem, ILVIVBool) ;
      if ~isempty( ILVIVBool ) && ILVIVBool % In line VIV
          p = WOM_IL( VprojRel1, VprojRel2, udotdotFrame1,udotdotFrame2,...
                      dimCharacteristic, nextTime, analysisSettings.deltaT, currElem ) ;
      else
          p=0;
      end

    else
      q = 0 ; p = 0;% No lift with circular cross section!
      % declare lift constant directions which are not taken into account (in this case the lift direction is updated)
      tlift1 = [] ; tlift2 = [] ;
    end
  else
    q = 2 ; p = 0;
    % declare lift constant directions which are not taken into account (in this case the lift direction is updated)
    tlift1 = [] ; tlift2 = [] ;
  end

  %% ----------------------------------------------------------------------- 
  % Compute the element fluid force by the equivalent virtual work theory
  % ------------------------------------------------------------------------
  fDragLiftPitchElem = zeros(12,1) ;

  for ind = 1 : length( xIntPoints )

    %The Gauss integration coordinate is:
    xGauss = l0/2 * ( xIntPoints( ind ) + 1 ) ;
    %Integrate for different cross section inner to the element
    fDragLiftPitchElem =  fDragLiftPitchElem ...
               +l0/2 * wIntPoints( ind ) ...
                 * integFluidForce( xGauss, ddotg, udotFlowElem,...
                                    l0, tl1, tl2, Rr,...
                                    chordVector', dimCharacteristic,...
                                    I3, O3, P, G, EE, L2, L3,...
                                    aeroCoefs, densityFluid, viscosityFluid,...
                                    VIVBool, q, p, constantLiftDir, uniformUdot, tlift1, tlift2, fluidFlowBool, ILVIVBool, uBEMbool ) ;

    if isnan( norm(fDragLiftPitchElem)), error(' drag force is NaN'), end

  end
  % ------------------------------------------------------------------------


  %% -----------------------------------------------------------------------

  % express aerodynamic force in ONSAS nomenclature  [force1 moment1 force2 moment2  ...];
  % ------------------------------------------------------------------------
  fDragLiftPitchElem = swtichToONSASBase( fDragLiftPitchElem ) ;
  % ------------------------------------------------------------------------

  %% -----------------------------------------------------------------------
  % Compute the element fluid force by the equivalent virtual work theory
  % -----------------------------------------------------------------------
  fAddedMassElem = addedMassForce( AMBool                              ,...
                                l0, elemCoords, elemCrossSecParams     ,...
                                analysisSettings.deltaT, nextTime      ,...
                                userFlowVel, densityFluid ) ;


  % -------------------------------

  fHydroElem =  fDragLiftPitchElem + fAddedMassElem ;
  % -----------------------------------------------------------------------

  %% -------------------------------------------------------------------- 
  % compute tangent matrix (dFagElem/du) using Central Difference
  % ---------------------------------------------------------------------
  % fHydroElem(udotdot, udot, u + iu) - fHydroElem
  if computeAeroStiffnessMatrix
    tMatHydroElemU = dispTangMatElem( fHydroElem                     ,...
                                    elemCoords, elemCrossSecParams   ,...
                                    Ue, Udote, Udotdote              ,...
                                    aeroCoefs, chordVector, aeroNumericalParams ,...
                                    analysisSettings, nextTime, currElem, uBEMbool ) ;
  else
    tMatHydroElemU = [] ;
  end
  % ---------------------------------------------------------------------

end


%% -----------------------------------------------------------------------
% This function returns the tangent matrix of the hydrodinamic force vector 
% with respect to u employing a simple central difference alg.
function dispTangMatElem = dispTangMatElem( fHydroElem                                ,...
                                            elemCoords, elemCrossSecParams            ,...
                                            Ue, Udote, Udotdote                       ,...
                                            aeroCoefs, chordVector, aeroNumericalParams ,... 
                                            analysisSettings , nextTime, currElem, uBEMbool )
  % disp("entre")
  % initialize aerodynamic tangent matrix
  dispTangMatElem = zeros(12,12) ;
  % numerical step to compute the tangets
  h = 1e-10           ;
  elemTypeAero(5) = 0 ; % set compute tangents to false
  for indexIncrementU = 1:12
    e_i = zeros(12,1) ;  e_i(indexIncrementU) = 1 ;
    % increment displacement
    UplusDeltaU = Ue + h * e_i   ;
    % compute forces with u + h*ei at the index indexIncrementU}
    fhydro_incU = frame_fluid_force( elemCoords                                ,...
                                      elemCrossSecParams                        ,...
                                      UplusDeltaU, Udote, Udotdote              ,...
                                      aeroCoefs, chordVector, aeroNumericalParams ,...
                                      analysisSettings,nextTime, currElem,...
                                      false, uBEMbool ) ;
    
    % central difference
    dispTangMatElem(:, indexIncrementU ) = ( fhydro_incU - fHydroElem ) / h ;
  end % endfor
end % end function
