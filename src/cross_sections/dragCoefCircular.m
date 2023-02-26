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

% This function retrevies drag for circular cross section  
% Reference for the drag formulation.
% https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0000722
% The function is suggetes by Frederick Gosselin in #issue 617 

function C_d = dragCoefCircular( betaRel, Re )

    C_d = 11 * Re^(-0.75) + 0.9 * (1- exp(- 1000 / Re) ) + 1.2 * (1 - exp(-(Re/4500)^(0.7)));

    if Re >= 2 * (10)^5
        warning("The drag coefficient relation breaks down as it does not account for the drag crisis.");
    end

end

