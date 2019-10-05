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


%Script for verification of the input variables definition.

tic

if plotParamsVector(1)>0
  fprintf('  - input variables verification ... ');
end

nmats  = length( hyperElasParams ) ;
nsecs  = length( secGeomProps    ) ;
nnodes = size(Nodes,1) ;
nelems = size(Conec,1) ;


if exist( 'inputONSASversion' ) == 0
  error('The inputONSASversion string variable was not defined.') ;
elseif (~isequal( inputONSASversion, ONSASversion)),
  error('input for different ONSAS version.');
end

if exist( 'problemName' ) == 0
  error('problemName variable was not defined.');
end


% -----------------------
% structural properties

if exist('hyperElasParams') == 0
  error('hyperElasParams variable was not defined.');
end

if exist( 'secGeomProps' ) == 0
  error('secGeomProps variable was not defined.');
end

if exist( 'nodalSprings' ) == 0
  error('A nodalSprings matrix must be defined.');
end

if exist( 'prescribedDisps' ) == 0
  prescribedDisps = [] ; 
elseif length(prescribedDisps)>0
    if size( prescribedDisps ) (2) ~= 3
      error('The prescribedDisps matrix must have 3 columns') ; 
    end
end  

if exist( 'Releases' ) == 0
  Releases = [] ;
end

if exist( 'plotsViewAxis' ) == 0
  plotsViewAxis = [] ;
end

if exist( 'bendStiff' ) == 0
  bendStiff = [] ;
end

if exist( 'Nodes' ) == 0
  error('The Nodes matrix was not defined.') ;
end

if exist( 'Conec' ) == 0
  error('The Conec matrix was not defined.') ;
end


if exist( 'sectPar' ) == 0
  sectPar = [0 0 ] ;
end 

if exist( 'LinBuckModeImperf' ) == 0
  LinBuckModeImperf = 0 ;
end 

if exist( 'imperFactor' ) == 0
  imperfactor = 1e-4 ;
end

if exist( 'dim' ) == 0
	dim = 3 ;
end
% ---------------------------------


% ---------------------------------
% loading parameters

if exist( 'loadFactorsFunc') == 0
  loadFactorsFunc = @(t) t ;
end

if exist( 'selfWeightBoolean') == 0
  selfWeightBoolean = 0 ;
else
	if ~selfWeightBoolean == 0
		if exist( 'rho' ) == 0, error( 'Density was not defined.' ) ; end
	end
end

if exist( 'userLoadsFilename') == 0
  userLoadsFilename = '' ;
end

if exist( 'unifLoad' ) == 0
  unifLoad = [] ;
end



% -----------------------




% -----------------------
% analysis settings

if exist( 'nonLinearAnalysisBoolean' ) == 0
  nonLinearAnalysisBoolean  = 1 ; 
end

if exist( 'dynamicAnalysisBoolean' ) == 0
  dynamicAnalysisBoolean  = 0 ; 
end

if exist( 'LBAAnalyFlag' ) == 0
  LBAAnalyFlag = 0 ;
elseif LBAAnalyFlag == 0 && LinBuckModeImperf>0;
  error('LBAAnalyFlag must be set 1 to add an LBA mode imperfection.') ; 
end

if ( exist( 'numericalMethodParams' ) == 0 ) && ( nonLinearAnalysisBoolean ~= 0 || dynamicAnalysisBoolean ~= 0 )
  error( 'numericalMethodParams must be defined by the user if a nonlinear/dynamic analysis is performed.') ; 
end


if ( nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0 )
  if ( exist( 'linearDeformedScaleFactor' ) == 0 ) 
    linearDeformedScaleFactor = 1 ;
  end

else

  if exist( 'linearDeformedScaleFactor' ) ~= 0
    fprintf(' WARNING: linearDeformedScaleFactor set in input but not considered by ONSAS.\n');
  end
  linearDeformedScaleFactor = 1 ;

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


if exist( 'numericalMethodParams' ) == 0 && ( nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0 )
  numericalMethodParams = [] ;
end


if exist( 'plotParamsVector' ) == 0
  plotParamsVector = [1] ;
else
	if length(plotParamsVector) == 2
		plotParamsVector = [plotParamsVector 0] ;
	end
end

if exist( 'printflag' ) == 0
  printflag = 0 ;
end

if exist( 'reportBoolean' ) == 0
  reportBoolean = 1 ;
end

if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  timeIncr = 0 ;
end

if dynamicAnalysisBoolean == 1
  if exist('rho')== 0, error('Density was not defined.');end

  if exist( 'nodalDamping') == 0
    nodalDamping = 0 ;
  end
  
  if exist( 'deltamassMat' ) == 0
    if sum( Conec(:,7) == 2 ) > 0
      error(' a value for deltamassMat is required.');
    else
      deltamassMat = 0 ;
    end
  end

else
  if exist('rho') == 0
    rho = 0 ;
  end  
end


if plotParamsVector(1)>0
  fprintf(' done.\n');
end

tVarVer = toc ;



% creates outputdir
outputdir = [ './output/' problemName '/' ] ;

if exist( './output/' ) ~= 7
  fprintf( '  - Creating directory ./output/ ...' );
  mkdir('./', './output/' );
  fprintf( ' done. \n' );
end

cd output
if exist( ['./' problemName '/' ] ) ~= 7
  fprintf( ['  - Creating directory ./output/' problemName '/ ...'] );
elseif exist( ['./' problemName '/' ] ) == 7
  if plotParamsVector(1)>0
    fprintf( ['  - Cleaning directory ./output/' problemName '/ ...'] );
  end
  confirm_recursive_rmdir (0)
  [aux,msg] = rmdir( problemName ,'s'); 
end
mkdir('./', ['./' problemName '/' ] );
if plotParamsVector(1)>0
  fprintf( ' done. \n' );
end
cd ..
% -----------------------

