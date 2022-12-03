% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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

% This function computes the fluid loads within the quasi-steady theory for co-rotational dynamic frame elements proposed by Lee, Battini 2014
function [fHydroElem, tMatHydroElemU] = frame_fluid_force( elemCoords            ,...
                                     elemCrossSecParams                       ,...
                                     Ue, Udote, Udotdote                      ,...
                                     aeroCoefs, elemTypeAero, analysisSettings,...
                                     nextTime, currElem, hydroTangBoolU )

  % Check all required parameters are defined
  assert( ~isempty( analysisSettings.fluidProps), ' empty analysisSettings.fluidProps.' )
  assert( ~isempty( elemTypeAero), ' empty elements.elemTypeAero.' )
  assert( ~isempty( aeroCoefs )  , ' empty elements.aeroCoefs '    )

  % Declare booleans for VIV phenomenon
  % set boolean to set constant lift direction in VIV problems
  % set boolean to set unfirom u dot
  global VIVBool
  global constantLiftDir
  global uniformUdot
  global AMBool
  % Implementation Booleans for internal test, baseBool changes the local angles computation
  baseBool = false ;

  % extract fluid properties
  densityFluid   = analysisSettings.fluidProps{1,1} ;
  viscosityFluid = analysisSettings.fluidProps{2,1} ;
  userFlowVel    = analysisSettings.fluidProps{3,1} ;

  % check user Flow Vel is not empty
  assert( ~isempty( userFlowVel ), 'empty user windvel' )

  % extract nonLinearity in aero force boolean
  geometricNonLinearAero = analysisSettings.geometricNonLinearAero ;

  % Boolean to compute the fluid force with ut = 0, this should be used if the fluid loads are computed in the reference
  % configuration and nonlinear effects are not considered.
  if ~geometricNonLinearAero
    Ue = zeros(12,1) ;
  end

  % fluid velocity at the nodes of the element evaluated in deformed configuration (spatial points):
  udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6), nextTime ) ;
  udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime ) ;
  % compact them into a single vector for the element
  udotFlowElem  = [udotFlowNode1; udotFlowNode2] ;

  % Elem reference coordinates:
  xs = elemCoords(:) ;

  % Load element properties to fluid loads
  % chord vector
  vecChordUndef    = elemTypeAero( 1:3 )'  ;
  % length of the chord vector
  dimCharacteristic = norm( vecChordUndef ) ;

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


  % auxiliary matrix created to project transversal velocity
  L2 = [ 0 0 0
         0 1 0
         0 0 1 ] ;
  % auxillary matrix created to rotate 90 degrees (this will define the lift force once the drag is computed)
  L3 = expon( [pi/2 0 0] ) ;

  % Extract points and weights for numGausspoints selected
  numGaussPoints = elemTypeAero(4);
  [xIntPoints, wIntPoints] = gaussPointsAndWeights( numGaussPoints ) ;

  % WOM computation call for cases with VIVbool equal to true
  if ~isempty( VIVBool ) && ~isempty( constantLiftDir ) && ~isempty( uniformUdot )

    if VIVBool && ( norm(udotFlowNode1)*norm(udotFlowNode2) ) > 0
      % extract the accelerations and the velocities of nodes in global coordinates
      % node 1
      udotFrame1    = Udote( 1:2:6 )      ;
      udotdotFrame1 = Udotdote( 1:2:6 )   ;
      % node 2
      udotFrame2    = Udote( 7:2:end )    ;
      udotdotFrame2 = Udotdote( 7:2:end ) ;

      % projected velocities at nodes 1 and 2 in deformed coordinates
      %node 1
      [VpiRel1_defCords, VpiRelPerp1_defCords, Vrel1_glob] = computeVpiRels( udotFlowNode1, udotFrame1,...
                                                                             Rroof1, Rr, L2, L3 ) ;
      %node 2
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
      if ~constantLiftDir
        tlift1 = VpiRelPerp1 / norm( VpiRelPerp1 ) ;
        tlift2 = VpiRelPerp2 / norm( VpiRelPerp2 ) ;
      else
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
      end
      % compute van der pol solution for current element
      q = WOMV4( VpiRel1, VpiRel2, udotdotFrame1, udotdotFrame2,...
                 tlift1, tlift2, dimCharacteristic, nextTime, analysisSettings.deltaT, currElem ) ;
    else
      q = 0 ; % No lift with circular cross section!
      % declare lift constant directions which are not taken into account (in this case the lift direction is updated)
      tlift1 = [] ; tlift2 = [] ;
    end
  else
    q = 2 ;
    % declare lift constant directions which are not taken into account (in this case the lift direction is updated)
    tlift1 = [] ; tlift2 = [] ;
  end

  % Compute the element fluid force by the equivalent virtual work theory
  fDragLiftPitchElem = zeros(12,1) ;
  for ind = 1 : length( xIntPoints )
    %The Gauss integration coordinate is:
    xGauss = l0/2 * ( xIntPoints( ind ) + 1 ) ;
    %Integrate for different cross section inner to the element
    fDragLiftPitchElem =  fDragLiftPitchElem ...
               +l0/2 * wIntPoints( ind ) * integFluidForce( xGauss, ddotg, udotFlowElem,...
                                                           l0, tl1, tl2, Rr,...
                                                           vecChordUndef, dimCharacteristic,...
                                                           I3, O3, P, G, EE, L2, L3,...
                                                           aeroCoefs, densityFluid, viscosityFluid,...
                                                           VIVBool, q,  constantLiftDir, uniformUdot, tlift1, tlift2 ) ;
  end

  % express aerodynamic force in ONSAS nomenclature  [force1 moment1 force2 moment2  ...];
  fDragLiftPitchElem = swtichToONSASBase( fDragLiftPitchElem ) ;
  % -------------------------------

  % Compute the element fluid force by the equivalent virtual work theory
  fAddedMassElem = addedMassForce( AMBool                              ,...
                                l0, elemCoords, elemCrossSecParams     ,...
                                analysisSettings.deltaT, nextTime      ,...
                                userFlowVel, densityFluid ) ;
  % -------------------------------



  fHydroElem =  fDragLiftPitchElem + fAddedMassElem;

  % --- compute tangent matrix (dFagElem/du) using Central Difference  ---
  % fHydroElem(udotdot, udot, u + iu) - fHydroElem
  tMatHydroElemU = [] ;
  if hydroTangBoolU
    tMatHydroElemU = dispTangMatElem( fHydroElem                     ,...
                                    elemCoords, elemCrossSecParams   ,...
                                    Ue, Udote, Udotdote              ,...
                                    aeroCoefs, elemTypeAero          ,...
                                    analysisSettings, nextTime, currElem ) ;
  end

  % -------------------------------


end





% This function returns the tangent matrix of the hydrodinamic force vector with respect to u
% employing a simple central difference alg.
function dispTangMatElem = dispTangMatElem( fHydroElem                                ,...
                                            elemCoords, elemCrossSecParams            ,...
                                            Ue, Udote, Udotdote                       ,...
                                            aeroCoefs, elemTypeAero, analysisSettings ,...
                                            nextTime, currElem )
  % initialize aerodynamic tangent matrix
  dispTangMatElem = zeros(12,12) ;
  % numerical step to compute the tangets
  h = 1e-10                  ;
  for indexIncrementU = 1:12
    e_i = zeros(12,1)        ;
    e_i(indexIncrementU) = 1 ;
    % increment displacement
    UplusDeltaU = Ue + h * e_i   ;
    % compute forces with u + hu at the index indexIncrementU
    fhydro_incU = frame_fluid_force( elemCoords                                ,...
                                      elemCrossSecParams                        ,...
                                      UplusDeltaU, Udote, Udotdote              ,...
                                      aeroCoefs, elemTypeAero, analysisSettings ,...
                                      nextTime, currElem, false ) ;
    % central difference
    dispTangMatElem(:,indexIncrementU) = ( fhydro_incU - fHydroElem ) / h ;
  end % endfor
end % end function
