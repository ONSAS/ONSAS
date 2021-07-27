% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro  
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


% This script declares several matrices and vectors required for the analysis. In this script, the value of important magnitudes, such as internal and external forces, displacements, and velocities are computed for step/time 0.


function [ U, Udot, Udotdot ] = initialCondsProcessing( nNodes )

% create velocity and displacements vectors
U       = zeros( 6*nNodes,   1 ) ;
Udot    = zeros( 6*nNodes,   1 ) ;
Udotdot = zeros( 6*nNodes,   1 ) ;

% adds non homogeneous initial conditions
%~ if length( nonHomogeneousInitialCondU0 ) > 0
  %~ for i = 1 : size( nonHomogeneousInitialCondU0, 1 ) % loop over rows of matrix
    %~ dofs = nodes2dofs(nonHomogeneousInitialCondU0(i, 1 ), 6 ) ;
    %~ U( dofs ( nonHomogeneousInitialCondU0 (i, 2 ) ) ) = ...
      %~ nonHomogeneousInitialCondU0 ( i, 3 ) ;
  %~ end
%~ end % if nonHomIniCond

%~ if length( nonHomogeneousInitialCondUdot0 ) > 0
  %~ if numericalMethodParams(1) >= 3
    %~ for i=1:size(nonHomogeneousInitialCondUdot0, 1)
      %~ dofs = nodes2dofs( nonHomogeneousInitialCondUdot0(i, 1), 6 ) ;
      %~ Udot( dofs( nonHomogeneousInitialCondUdot0(i, 2 ))) = ...
        %~ nonHomogeneousInitialCondUdot0(i, 3 );
    %~ end
  %~ else
    %~ warning(' velocity initial conditions set for a static analysis method' ) ;
  %~ end
%~ end


% computation of initial acceleration for some cases
% ---------------------------------------------------

% --- initial tangent matrices ---
%~ [ mats ] = assembler ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, dynamicAnalysisBoolean, 2, Udotdott, booleanConsistentMassMat ) ;

%~ systemDeltauMatrix = mats{1}

%~ stop
%~ if dynamicAnalysisBoolean == 1,

  %~ massMat    = mats{2} ;

  %~ % --- computation of initial Udotdott for truss elements only!!!
  %~ Fext = computeFext( constantFext, variableFext, loadFactors(1), userLoadsFilename ) ;

  %~ Udotdott (neumdofs) = massMat( neumdofs, neumdofs ) \ ( Fext(neumdofs) -Fintt( neumdofs ) ) ;
