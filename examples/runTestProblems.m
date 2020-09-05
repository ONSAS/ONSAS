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

addpath('../sources/');

if isThisOctave
  fileslist = readdir('./');
else
  fileslist = {} ;
  auxMatlab = dir('*.*');
  for k=1:length( auxMatlab )
    fileslist{k} = auxMatlab(k).name ;
  end
end

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

verifBoolean = 1 ;

while current <= totalRuns && verifBoolean == 1

  save('-mat', 'exData.mat','current','totalRuns', 'keyfiles' );

  run( [ './' keyfiles{current} ] ) ;

  pause(0.5)

  if verifBoolean
    fprintf([' === test ' problemName ' problem:  PASSED === \n\n']);
  else
    fprintf([' === test ' problemName ' problem FAILED   === \n\n']);
  end
  
  load('exData.mat');
  current = current + 1 ; 
end

delete('exData.mat')

fidVB = fopen('auxVerifBoolean.dat','w') ;

fprintf( fidVB, sprintf('%1i',verifBoolean ) );
fclose( fidVB );

