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

%script for extracting variables from struct data sets.

arrfie    = cell(6,1);
arrfie{1} = modelCurrSol;
arrfie{2} = BCsData ;
arrfie{3} = modelProperties  ;
arrfie{4} = 'modelCurrSol';
arrfie{5} = 'BCsData' ;
arrfie{6} = 'modelProperties'  ;

for k=1:3
  fieldNames = fieldnames( arrfie{k} ) ;
  for i=1:length(fieldNames)
    eval( [ fieldNames{i} '= getfield(' arrfie{k+3} ',''' fieldNames{i} ''');' ] )
  end
end
