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

% run all test examples from examples folder.

keyword = 'test.m' ;
%~ keyword = 'example' ;

fileslist = readdir('./');
keyfiles  = {} ;

totalRuns = 0 ;

for i=1:length(fileslist)
  if length( strfind( fileslist{i}, keyword ) ) > 0
    totalRuns = totalRuns +1 ; 
    keyfiles{totalRuns} = fileslist{i} ;
  end
end

keyfiles

current       = 1 ;
totalRuns

while current <= totalRuns

  save('-mat', 'exData.mat','current','totalRuns', 'keyfiles' );

  run( [ './' keyfiles{current} ] ) ;
  pause(0.5)
  fprintf([' === test ' problemName ' problem executed === \n\n']);
  
  load('exData.mat');
  current = current + 1 ; 
end

delete('exData.mat')
