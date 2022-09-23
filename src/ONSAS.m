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

function [ matUs, loadFactorsMat, cellFint ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

%mdFirst the input structs are converted to structs with the model information
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%mdAfter that the structs are used to perform the numerical time analysis
[ matUs, loadFactorsMat, cellFint ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

%md Finally the report is generated
outputReport( modelProperties.outputDir, modelProperties.problemName )
