% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

arrfie    = cell(8,1);
arrfie{1} = modelCurrState;
arrfie{2} = BCsCurrState ;
arrfie{3} = auxIO  ;
arrfie{4} = modelProperties  ;
arrfie{5} = 'modelCurrState';
arrfie{6} = 'BCsCurrState' ;
arrfie{7} = 'auxIO'  ;
arrfie{8} = 'modelProperties'  ;

for k=1:4
  fieldNames = fieldnames( arrfie{k} ) ;
  for i=1:length(fieldNames)
    eval( [ fieldNames{i} '= getfield(' arrfie{k+4} ',''' fieldNames{i} ''');' ] )
  end
end
