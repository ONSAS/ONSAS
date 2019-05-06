%Script for reading of the input file to be used in the analysis. During the execution of this script the user must choose an input file from a list of available files which is shown in the terminal.

%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

fprintf('==============================================\n');
fprintf( [ 'Welcome to ONSAS v' ONSASversion '.\n' ] )
fprintf('==============================================\n');
fprintf( [ 'Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli \n' ] ) ;
fprintf( [ 'This program comes with ABSOLUTELY NO WARRANTY. \n' ] ) ;
fprintf( [ 'This is free software, and you are welcome to redistribute it under certain conditions; read COPYNG file for more details. \n' ] ) ;

close all, clear -x previouslyDefinedSelectedFileVar ONSASversion
fileslist = readdir('./input');


if exist( 'previouslyDefinedSelectedFileVar' ) == 0
  
  fprintf('The content of the input folder is: \n')
  fprintf('--------------\n');
  
  for i=1:length(fileslist)
    string = fileslist{i};
    fprintf('  file %3i: %s\n',i, string); 
  end

  fprintf('--------------\n');
  
  selectedfile = input('  - Write the number of the selected input file:') ;

else
  selectedfile = previouslyDefinedSelectedFileVar ;

end

run( fileslist{ selectedfile }(1:end-2) ) ;

if plotParamsVector(1)>0
  fprintf(['  - Reading variables variables from input file: \n      ' fileslist{ selectedfile }  '...' ] );
end


if plotParamsVector(1)>0
  fprintf([' done.\n'] );
end
