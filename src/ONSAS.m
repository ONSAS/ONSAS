% ==============================================================================
% --------     ONSAS: an Open Non-linear Structural Analysis System     --------
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


function ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

% convert input data to model structures
% --------------------------------------
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

modelCurrSol
modelProperties
BCsData

stop
% performe the time analysis
% --------------------------
ONSAS_solve( modelCurrSol, modelProperties, BCsData )


if otherParams.screenOutputBool
  fprintf([ '|-------------------------------------------------|\n'])
  fprintf(  '|  ONSAS finished in: %7.1e seconds /%5.2f mins |\n', totalTime, totalTime/60 )
  fprintf([ '|=================================================|\n\n\n'])
end

