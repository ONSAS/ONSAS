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

function [P1, P2, N, N1, N2 ] = corotVecMatAuxDyn(...
N1, N2, N3, N4, N5, N6, N7, N8, tl1, tl2, Gaux, I3, O3, P )
  

  % Auxiliary shape function matrices variables for the cross section:
  P1 = [  0   0   0   0   0    0  ; ...
          0   0   N3  0   0    N4 ; ...
          0  -N3  0   0   -N4  0  ] ; % Eq.(38)  T-N Le J.-M. Battini et al 2014

  P2 = [  N1  0   0   N2  0    0  ; ...
           0  N5  0   0   N6   0  ; ...
           0  0   N5  0   0    N6 ] ; % Eq.(39)  T-N Le J.-M. Battini et al 2014

  N  = [ N1 * I3   O3   N2 * I3   O3 ] ; % Eq.(51)  T-N Le J.-M. Battini et al 2014

  ul = P1 * [ tl1; tl2 ]                   ; % Eq.(38)  T-N Le J.-M. Battini et al 2014

  H1 = N + P1 * P - 1 * skew( ul ) * Gaux' ; % Eq.(59)  T-N Le J.-M. Battini et al 2014
  H2 = P2 * P + Gaux'                      ; % Eq.(72)  T-N Le J.-M. Battini et al 2014

end