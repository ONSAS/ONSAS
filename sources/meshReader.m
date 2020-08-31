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


function [ nodesMat, conecMat ] = meshReader( mshFilename )

fid = fopen( mshFilename ,'r') ;

maxLengthLine = 200 ;

% ---- header reading --------------------------
X = fgets(fid); X = fgets(fid);
X = fgets(fid); X = fgets(fid);
X = fgets(fid); X = fgets(fid);
X = fgets(fid);

nnodes = str2num(fgets(fid)) ;

nodesMat = fscanf(fid,'%g %g %g %g\n' ,[4 nnodes])' ;
nodesMat(:,4) = [] ;

X = fgets(fid);

nElems = str2num( fgets(fid) ) ;

conecMat = fscanf(fid,'%g %g %g %g\n' ,[4 nElems ])' ;

fclose(fid);
