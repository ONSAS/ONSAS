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

function [ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

%md sets the current version
ONSASversion = '0.1.10'  ;

%md set defaults
analysisSettings.Utp10       = checkOrSetDefault ( analysisSettings, 'iniMatUs'        , [] ) ;
otherParams.screenOutputBool = checkOrSetDefault ( otherParams     , 'screenOutputBool', 1  ) ;
otherParams.plotsFormat      = checkOrSetDefault ( otherParams     , 'plotsFormat', []  ) ;
otherParams.nodalDispDamping = checkOrSetDefault ( otherParams     , 'nodalDispDamping', 0 ) ;

%md welcome message function
welcomeMessage(ONSASversion, otherParams );

% creates outputdir in current location
% -------------------------------------
outputDir = [ './output/' otherParams.problemName '/' ] ;
createOutputDir( outputDir, otherParams ) ;
otherParams.outputDir = outputDir ;

%md process boundary conds information and elements
[ Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, ...
  diriDofs, neumDofs, KS, userLoadsFilename ] = boundaryCondsProcessing( mesh, ...
                           materials, elements, boundaryConds, initialConds ) ;

% process initial conditions
% --------------------------
[ U, Udot, Udotdot ] = initialCondsProcessing( size(Nodes,1) ) ;

currTime         = 0 ; timeIndex        = 1 ; convDeltau      = zeros( size(U) ) ;
timeStepIters    = 0 ; timeStepStopCrit = 0 ;

%md call assembler
[~, Stress ] = assembler ( Conec, elements, Nodes, materials, KS, U, Udot, Udotdot, analysisSettings, [ 0 1 0 ], otherParams.nodalDispDamping ) ;

systemDeltauMatrix = computeMatrix( Conec, elements, Nodes, materials, KS, analysisSettings, U, Udot, Udotdot, neumDofs, otherParams.nodalDispDamping ) ;

if isfield(otherParams,'spitMatrices') && otherParams.spitMatrices == true
  cd output
  save('-mat', 'matrices.mat', 'systemDeltauMatrix' );
  cd ..
  stop
end

[ Fext, vecLoadFactors ] = computeFext( factorLoadsFextCell, loadFactorsFuncCell, analysisSettings, 0, length(U), userLoadsFilename ) ;

%md prints headers for solver output file
printSolverOutput( outputDir, otherParams.problemName, 0                  ) ;
printSolverOutput( outputDir, otherParams.problemName, [ 2 timeIndex currTime 0 0 ] ) ;


nTimes = round( analysisSettings.finalTime / analysisSettings.deltaT )

% if length( otherParams.plotParamsVector ) > 1
%   nplots = min( [ nTimes otherParams.plotParamsVector(2) ] ) ;
% else
%   % default value: all
   nplots = nTimes ;
% end

timesPlotsVec = round( linspace( 1, nTimes, nplots )' ) ;

%md compress model structs
[ modelCurrSol, modelProperties, BCsData ] = modelCompress( ...
  timeIndex, currTime, U, Udot, Udotdot, Stress, convDeltau, systemDeltauMatrix, ...
  timeStepStopCrit, timeStepIters, factorLoadsFextCell, loadFactorsFuncCell, neumDofs, ...
  KS, userLoadsFilename, Nodes, Conec, materials, elements, analysisSettings, ...
  outputDir, vecLoadFactors, otherParams.problemName, otherParams.plotsFormat, ...
  timesPlotsVec, otherParams.nodalDispDamping );



%md writes vtk file
if strcmp( modelProperties.plotsFormat, 'vtk' )
  vtkMainWriter ( modelCurrSol, modelProperties )
end


if exist( 'controlDofs') ==0
  controlDofs = [] ;
  controlDofsAndFactors = [] ;
end

if length( controlDofs ) > 0
  controlDofsAndFactors = zeros( size( controlDofs,1 ) , 2 ) ;

  % control dof info
  for i=1:size(controlDofs,1)
    aux                = nodes2dofs( controlDofs(i,1), 6 ) ;
    controlDofsAndFactors(i,:) = [ aux( controlDofs(i, 2) ) controlDofs(i,3) ] ;
  end
end



% =========================================
% function for creation of output directory
% -----------------------------------------
function createOutputDir( outputDir, otherParams )

if exist( './output/' ) ~= 7
  if otherParams.screenOutputBool
    fprintf( '  - Creating directory ./output/ ...' );
  end

  mkdir('./', './output/' );

  if otherParams.screenOutputBool
    fprintf( ' done. \n' );
  end
end

% -----------------
if exist( outputDir ) == 7 % problemName is a directory
  % the content is erased
  if otherParams.screenOutputBool
    fprintf( ['|  - Cleaning output directory ...'] ) ;
  end
  if isThisOctave
    confirm_recursive_rmdir(0)
  end

  % delete
  [aux, msg] = rmdir( outputDir ,'s') ;

  % create empty
  mkdir( outputDir );

elseif exist( ['./' otherParams.problemName '/' ] ) ~= 7 % problemName is not a directory
  % it is created
  if otherParams.screenOutputBool
    fprintf( ['|  - Creating output directory ...'] ) ;
  end
  mkdir( outputDir );
end
if otherParams.screenOutputBool
  fprintf( ' done. \n' );
end


% =========================================
% function for welcome message
% -----------------------------------------
function welcomeMessage( ONSASversion, otherParams )

if otherParams.screenOutputBool
  fprintf([ '\n' ...
            '|=================================================|\n' ...
            '|         _ _             _ _     _ _     _ _     |\n' ...
            '|       /    /  /|   /  /       /    /  /         |\n' ...
            '|      /    /  / |  /  /_ _    /_ _ /  /_ _       |\n' ...
            '|     /    /  /  | /       /  /    /       /      |\n' ...
            '|    /_ _ /  /   |/   _ _ /  /    /   _ _ /       |\n' ...
            '|                                                 |\n' ...
            '|-------------------------------------------------|\n' ] );
  fprintf([ '| Welcome to ONSAS v' ONSASversion '.                       |\n' ...
            '| This program comes with ABSOLUTELY NO WARRANTY. |\n' ...
            '| Please read the COPYING.txt and README.md files |\n' ...
            '|-------------------------------------------------|\n'] ) ;
  fprintf( ['| Solving problem:  ' otherParams.problemName '\n' ] ) ;
end
