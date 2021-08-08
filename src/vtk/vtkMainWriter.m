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


%md function for writing vtk files of deformed configurations of structures.
%md Creates the file filename with the nodes coordinates given in nodes,
%md the conectivity given in conect and with the point and element data
%md given in cellPointData and cellCellData, respectively.

function vtkMainWriter( modelCurrSol, modelProperties )

plotInd = find( modelProperties.timesPlotsVec == modelCurrSol.timeIndex ) ;

%md if the current time index is not in the plts indexes vector then ends execution
if length( plotInd ) == 0, return, end

%md filname counter starts in zero
filename = [ modelProperties.outputDir modelProperties.problemName '_' sprintf('%04i', plotInd) '.vtk']

%md data conversion
[ vtkNodes, vtkConec, cellPointData, cellCellData ] = vtkDataFormater( modelCurrSol, modelProperties ) ;

%md the function __vtkWriter__ writes the vtk file. it has no outputs and recieves vtk formatted nodes, conectivity and cell and point data.
vtkFileWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
