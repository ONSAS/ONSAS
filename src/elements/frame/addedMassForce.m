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

% This function returns the hydrodinamic mass force of the element in global coordinates.

function fam = addedMassForce( AMBool                                   ,...
                               l0, elemCoords, elemCrossSecParams       ,...
                               deltaT, nextTime, userFlowVel, densityFluid ) ;

  Aelem  = crossSectionProps( elemCrossSecParams, 0.0 ) ; % the 0.0 density does not affect A value

  if ~isempty( AMBool ) && AMBool     % linear(node1_x)  angular(node1_x)   % linear(node2_z)  angular(node2_z)
    % fill fluid acceleration vector [udotdot_f_x_1, wdotdot_f_x_1 .... udotdot_f_z_2, wdotdot_f_x_12 ]
    Udotdotflow = zeros(12, 1);
    % compute fluid element acceleration [udotdot_f_x_1, udotdot_f_y_1, udotdot_f_z_1, udotdot_f_x_2 ....]
    ddUf = computeddUf(nextTime, deltaT, userFlowVel,  elemCoords);
    Udotdotflow(1:2:12) = ddUf(1:6); % Irotationnal flow

    % lumped add mass formlation:
    % circular cross section implementation
    assert(elemCrossSecParams{1}(1:end) == 'circle')
    Ca            = 1          ;
    elementVolume = l0 * Aelem ;
    massNodeAdded = (1 + Ca) * densityFluid * elementVolume / 2 ;
    fam = massNodeAdded * Udotdotflow(1:12);

  else
    fam = zeros(12, 1);
  end
end
