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

% This function returns the fluid acceleration in global coordinates
function ddUf = computeddUf(nextTime, dt, userFlowVel,  elemCoords)
   t0 = (nextTime-dt);
   t1 = nextTime;
   udotdotFlowNode1 = (feval(userFlowVel, elemCoords(1:3)', t1) - feval(userFlowVel,  elemCoords(1:3)', t0))/dt ;
   udotdotFlowNode2 = (feval(userFlowVel, elemCoords(4:6)', t1) - feval(userFlowVel, elemCoords(4:6)', t0))/dt ;
   ddUf = [udotdotFlowNode1' udotdotFlowNode2'];
end
