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

% Script for verification of the input variables definition. Default
% values ar assigned.

if exist('booleanScreenOutput') == 0 || booleanScreenOutput

  % default value
  booleanScreenOutput = 1 ;

  fprintf([ '\n\n' ...
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
            '| Read files COPYING.txt and README.md for more   |\n' ...
            '| information.                                    |\n' ...
            '|-------------------------------------------------|\n'] ) ;
end

if booleanScreenOutput
  fprintf( ['| Solving problem:  ' problemName '\n' ] ) ;
  fprintf(  '|  - input variables verification ... ') ;
end
if exist( 'loadsParams' ) == 0
  loadsParams = {} ; 
end



% --- verification of relevant variables ---
checkVarNamesList = { 'problemName', 'Nodes', 'Conec', 'dirOnsas', ...
                      'materialsParams', ...
                      'elementsParams', ...
                      'loadsParams', ...
                      'springsParams' } ;

for j = 1:length(checkVarNamesList)
  varName = checkVarNamesList{j} ;
  if exist( varName, 'var' ) == 0,
    error([ varName ' variable was not defined.'] );
  end
end
% ------------------------------------------

if exist( 'crossSecsParams' ) == 0
  crossSecsParams = {} ; 
end


if exist( 'numericalMethodParams' ) == 0
  numericalMethodParams = [ 0 ] ;
end
                      

if exist( 'BooleanSelfWheight' ) == 0
  BooleanSelfWheight = 0 ; 
end


% ===  Conversion conec cell to matrix format... to improve in the future.... ===
if iscell( Conec )
  aux    = Conec ;
  nElems = size(  aux,    1 ) ; 
  Conec  = zeros( nElems, 9 ) ;
  for i=1:nElems
    aux2                   = aux{i,1} ;  
    auxnnodes              = length( aux2) - 5  ;
    Conec ( i,1:auxnnodes) = aux2( (5+1):(5+auxnnodes) ) ;
    Conec ( i,5:9        ) = aux2( (  1):(5          ) ) ;
  end
  clear aux aux2
end
% ====================================

[ Conec, nodalVariableLoads, nodalConstantLoads, nodalSprings ] = ...
  conversionLoadsSprings ( Nodes, Conec, ...
                          materialsParams, ...
                          elementsParams, ...
                          loadsParams, ...
                          crossSecsParams, ...
                          springsParams, ...
                          BooleanSelfWheight...
                        ) ;

nMats  = length( materialsParams     ) ;
nSecs  = size(   crossSecsParams, 1  ) ;
nNodes = size(   Nodes,           1  ) ;
nElems = size(   Conec,           1  ) ;

% -----------------------
% default values
 
if exist( 'prescribedDispsMat' ) == 0
  prescribedDispsMat = [] ; 
end

%% default variables
%~ if exist( 'Releases' ) == 0
  %~ Releases = [] ;
%~ end

if exist( 'plotsViewAxis' ) == 0
  plotsViewAxis = [] ;
end

if  exist( 'nodalDispDamping' ) == 0
  nodalDispDamping = 0 ;
end

if exist( 'loadFactorsFunc') == 0
  loadFactorsFunc = @(t) t ;
end

if exist( 'userLoadsFilename') == 0
  userLoadsFilename = '' ;
end
% -----------------------


% -----------------------
% analysis settings

if exist( 'stabilityAnalysisBoolean' ) == 0
  stabilityAnalysisBoolean = 0 ;
end

if ( exist( 'deformedScaleFactor' ) == 0 ) 
  deformedScaleFactor = 1 ;
end

if exist( 'analyticSolFlag' ) == 0
  analyticSolFlag = 0 ;
else
	if analyticSolFlag ~= 0
		if exist( 'analyticCheckTolerance' ) == 0
			analyticCheckTolerance = 1e-8 ;
		end
	end
	if analyticSolFlag == 1 || analyticSolFlag == 2
		if exist( 'analyticFunc' ) == 0
			error('analytic Function must be defined.')
		else
			analytSol = [] ;
		end
	end	
	if analyticSolFlag == 3 || analyticSolFlag == 4 || analyticSolFlag == 5 
		if exist( 'analytSol' ) == 0
			error('analytic Values analytSol must be defined.')
		else
			analyticFunc = [] ;
			loadFactors = 0 ;
			timesVec = 0 ;
		end
	end
end


if exist( 'plotParamsVector' ) == 0
  plotParamsVector = [ 1 ] ;
end

if exist( 'storeBoolean' ) == 0
  storeBoolean = 0 ;
end

if exist( 'printFlag' ) == 0
  printFlag = 0 ;
end

if exist( 'reportBoolean' ) == 0
  reportBoolean = 0 ;
end

if exist( 'consMatFlag' ) == 0
  global consMatFlag
  consMatFlag = 2 ;
end

if booleanScreenOutput
  fprintf(' done.\n');
end

% creates outputdir
outputDir = [ './output/' problemName '/' ] ;

if exist( './output/' ) ~= 7
  if booleanScreenOutput
    fprintf( '  - Creating directory ./output/ ...' );
  end

  mkdir('./', './output/' );

  if booleanScreenOutput
    fprintf( ' done. \n' );
  end
end

% -----------------
if exist( outputDir ) == 7 % problemName is a directory
  % the content is erased
  if booleanScreenOutput
    fprintf( ['|  - Cleaning output directory ...'] ) ;
  end
  if isThisOctave
    confirm_recursive_rmdir(0)
  end
  
  % delete
  [aux, msg] = rmdir( outputDir ,'s') ;
  
  % create empty
  mkdir( outputDir );

elseif exist( ['./' problemName '/' ] ) ~= 7 % problemName is not a directory
  % it is created
  if booleanScreenOutput
    fprintf( ['|  - Creating output directory ...'] ) ;
  end
  mkdir( outputDir );
end
if booleanScreenOutput
  fprintf( ' done. \n' );
end
% ------------------

if exist( 'nonHomogeneousInitialCondU0') ==0
  nonHomogeneousInitialCondU0 = [] ;
end

if exist( 'nonHomogeneousInitialCondUdot0') ==0 
  nonHomogeneousInitialCondUdot0 = [] ;
end


if length( controlDofs ) > 0
  controlDofsAndFactors = zeros( size( controlDofs,1 ) , 2 ) ;
  
  % control dof info
  for i=1:size(controlDofs,1)
    aux                = nodes2dofs( controlDofs(i,1), 6 ) ;
    controlDofsAndFactors(i,:) = [ aux( controlDofs(i, 2) ) controlDofs(i,3) ] ; 
  end
end


coordsElemsMat = zeros(nElems,4*6) ; % 6 dofs per node, maximum 4 nodes per element

for i = 1 : nElems
  % obtains nodes and dofs of element
  nodeselem = Conec(i, find(Conec(i,1:4)>0) )' ;
  dofselem  = nodes2dofs( nodeselem , 6 ) ;
  for j=1:length(nodeselem)
    coordsElemsMat( i, (j-1)*6+[1:2:5] ) = Nodes( nodeselem(j), : ) ;
  end
end

% ---------------- load vectors assembly -----------------------
variableFext = zeros( 6*nNodes , 1 );
constantFext = zeros( 6*nNodes , 1 );

if exist( 'nodalVariableLoads' ) ~= 0
  for i=1:size(nodalVariableLoads,1)
    aux = nodes2dofs ( nodalVariableLoads(i,1), 6 ) ;
    variableFext( aux ) = variableFext( aux ) + nodalVariableLoads(i,2:7)' ;
  end
end

if exist( 'nodalConstantLoads' ) ~= 0
  for i=1:size(nodalConstantLoads,1)
    aux = nodes2dofs ( nodalConstantLoads(i,1), 6 ) ;
    constantFext( aux ) = constantFext( aux ) + nodalConstantLoads(i,2:7)' ;
  end
end

tangentMatricesCell = cell(2,1) ;
contProgr           = 0 ; % counter for progress bar
itersPerTimeVec     = 0 ;


if numericalMethodParams(1) == 0
  finalTime = 1 ;
  nTimes    = 2 ;

elseif numericalMethodParams(1) <= 2
  finalTime = numericalMethodParams(5) ;
  nTimes = numericalMethodParams(6) +1 ;

elseif numericalMethodParams(1) > 2
  finalTime = numericalMethodParams(3);
  nTimes = round( numericalMethodParams(3) / numericalMethodParams(2) ) ;
end

if length( plotParamsVector ) > 1
  nplots = min( [ nTimes plotParamsVector(2) ] ) ;
else
  % default value: all
  nplots = nTimes ;
end

timesPlotsVec = round( linspace(1, nTimes, nplots ) ) ;


% -- conver to matrices to store M., E. & C. in the struct type model ---
materialsParamsMat = [] ;
for i = 1 : length( materialsParams )
  materialsParamsMat (i, 1:length( materialsParams{i} ) ) = materialsParams{i} ;
end
elementsParamsMat = [] ;
for i = 1 : length( elementsParams )
  elementsParamsMat (i, 1:length( elementsParams{i} ) ) = elementsParams{i} ;
end
crossSecsParamsMat = [] ;
for i = 1 : length( crossSecsParams )
  crossSecsParamsMat (i, 1:length( crossSecsParams{i} ) ) = crossSecsParams{i} ;
end

% ----------------------------------------------------------------------
