% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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

function [ Fext, vecLoadFactors ] = computeFext( factorLoadsFextCell, loadFactorsFuncCell, analysisSettings, evalTime, lengthFext, userLoadsFilename )

Fext = zeros( lengthFext, 1 ) ;

vecLoadFactors = [] ;

for i=1:length( factorLoadsFextCell )

  if ~isempty( factorLoadsFextCell{i} )
    
    vecLoadFactors(i) = loadFactorsFuncCell{i}(evalTime) ;

    Fext  = Fext + vecLoadFactors(i) * factorLoadsFextCell{i}  ;
  end
end

if ~isempty( userLoadsFilename )
  Fext = Fext + feval( userLoadsFilename, evalTime )  ;
end
