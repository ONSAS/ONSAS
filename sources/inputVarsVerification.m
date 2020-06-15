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

  fprintf([ '|=================================================|\n' ...
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
  fprintf('  - input variables verification ... ') ;
end

% --- verification of relevant variables ---
checkVarNamesList = { 'problemName', 'Nodes', 'Conec', 'dirOnsas', ...
                      'materialsParams', 'crossSecsParams', 'nodalSprings', ...
                      'numericalMethodParams' } ;
for j = 1:length(checkVarNamesList)
  varName = checkVarNamesList{j} ;
  if exist( varName, 'var' ) == 0,
    error([ varName ' variable was not defined.'] );
  end
end
% ------------------------------------------

nMats  = length( materialsParams    ) ;
nSecs  = size( crossSecsParams, 1  ) ;
nNodes = size( Nodes,           1  ) ;
nElems = size( Conec,           1  ) ;

% -----------------------
% default values

if exist( 'prescribedDisps' ) == 0
  prescribedDisps = [] ; 
elseif length(prescribedDisps)>0
  if size( prescribedDisps, 2) ~= 3
    error('The prescribedDisps matrix must have 3 columns') ; 
  end
end  

%% default variables
if exist( 'Releases' ) == 0
  Releases = [] ;
end

if exist( 'plotsViewAxis' ) == 0
  plotsViewAxis = [] ;
end

if  exist( 'nodalDamping' ) == 0
  nodalDamping = [] ;
end

if exist( 'sectPar' ) == 0
  sectPar = [0 0 ] ;
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

if exist( 'booleanConsistentMassMat' ) == 0
  booleanConsistentMassMat = 1 ;
end

if exist( 'analyticSolFlag' ) == 0
  analyticSolFlag = 0 ;
else
	if analyticSolFlag ~= 0
		if exist( 'analyticCheckTolerance' ) == 0
			error('analyticCheckTolerance must be defined.')
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

if exist( 'octaveBoolean' ) == 0
  octaveBoolean = 1 ;
end

if exist( 'printflag' ) == 0
  printflag = 0 ;
end

if exist( 'reportBoolean' ) == 0
  reportBoolean = 1 ;
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
    fprintf( ['  - Cleaning output directory ...'] ) ;
  end
  if octaveBoolean
    confirm_recursive_rmdir (0)
  end
  [aux,msg] = rmdir( problemName ,'s'); 

elseif exist( ['./' problemName '/' ] ) ~= 7 % problemName is not a directory
  % it is created
  if booleanScreenOutput
    fprintf( ['  - Creating output directory ...'] ) ;
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

tangentMatricesCell = cell(2,1) ;



coordsElemsMat = zeros(nElems,4*6) ; % 6 dofs per node, maximum 4 nodes per element

for i = 1 : nElems
  % obtains nodes and dofs of element
  nodeselem = Conec(i, find(Conec(i,1:4)>0) )' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  for j=1:length(nodeselem)
    coordsElemsMat( i, (j-1)*6+[1:2:5] ) = Nodes( nodeselem(j), : ) ;
  end
end

% ---------------- load vectors assembly -----------------------
variableFext = zeros( 6*nnodes , 1 );
constantFext = zeros( 6*nnodes , 1 );

if exist( 'nodalVariableLoads' ) ~= 0
  for i=1:size(nodalVariableLoads,1)
    aux = nodes2dofs ( nodalVariableLoads(i,1), ndofpnode ) ;
    variableFext( aux ) = variableFext( aux ) + nodalVariableLoads(i,2:7)' ;
  end
end

if exist( 'nodalConstantLoads' ) ~= 0
  for i=1:size(nodalConstantLoads,1)
    aux = nodes2dofs ( nodalConstantLoads(i,1), ndofpnode ) ;
    constantFext( aux ) = constantFext( aux ) + nodalConstantLoads(i,2:7)' ;
  end
end

% ------------------------------------------------------------
if exist( 'nodalConstantLoads' ) ~= 0 || exist( 'nodalVariableLoads' ) ~= 0
  [maxNorm2F, visualloadfactor] = visualLoadFac( strucSize( Nodes ) , variableFext, constantFext, nnodes) ;
else
  error( ' user loads not included yet in computation of visual load factor ') ;
end


%~ cellStress = [] ;
%~ matNts = [] ;
%~ matUts = [] ;

%~ contProgr = 0 ;

%~ if dynamicAnalysisBoolean == 0
  %~ deltaT    = numericalMethodParams(5)/nLoadSteps ;
  %~ finalTime = numericalMethodParams(5) ;
%~ else
  %~ deltaT = timeIncr;
  
%~ end
