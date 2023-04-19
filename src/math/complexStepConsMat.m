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
 
function ConsMat = complexStepConsMat( stressFun, consParams, Egreen )

  compStep = 1e-10 ;

  for i=1:6

    direction = zeros( 3,3) ;

    if i<4
      direction(i,i) = 1 ;
    elseif i==4
      direction(2,3) = 1 ;
      direction(3,2) = 1 ;
    elseif i==5
      direction(1,3) = 1 ;
      direction(3,1) = 1 ;
    elseif i==6
      direction(1,2) = 1 ;
      direction(2,1) = 1 ;
    end

    EgreenComp = Egreen  + compStep * direction * j ;

    Scomp = feval( stressFun, consParams, EgreenComp, 0 ) ;

    dsde = mat2voigt( ( imag( Scomp ) / compStep ) , 1 ) ;

    ConsMat(1:3,i) = dsde(1:3) ;
    ConsMat(4:6,i) = dsde(4:6)*0.5 ;
    %~ ConsMat(:,i) = dsde ;
  end
