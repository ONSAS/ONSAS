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
 
%md function for writing vtk files of deformed configurations of structures.
%md Creates the file filename with the nodes coordinates given in nodes,
%md the conectivity given in conect and with the point and element data
%md given in cellPointData and cellCellData, respectively.

function vtkMainWriter( modelCurrSol, modelProperties )

plot_ind_float = modelCurrSol.currTime ...
            / ( modelProperties.plots_deltaTs_separation*modelProperties.analysisSettings.deltaT ) ;

plot_ind_round = round( plot_ind_float ) ;

if abs( plot_ind_float-plot_ind_round) < 1e-10,

  plotInd = plot_ind_round ;

  %md filname counter starts in zero
  filename = [ modelProperties.outputDir modelProperties.problemName '_' sprintf('%04i', plotInd) '.vtk'] ;

  %fprintf( [ '  writing vtk file ' modelProperties.problemName '_' sprintf('%04i', plotInd) '.vtk\n'] ) ;
  %md nodes and data conversion
  [ vtkNodes, vtkConec , vtkPointDataCell, vtkCellDataCell ] = vtkDataConversion( modelCurrSol, modelProperties ) ;
  %md the function __vtkWriter__ writes the vtk file. it has no outputs and recieves vtk formatted nodes, conectivity and cell and point data.
  vtkFileWriter( filename, vtkNodes, vtkConec , vtkPointDataCell, vtkCellDataCell ) ;

end
