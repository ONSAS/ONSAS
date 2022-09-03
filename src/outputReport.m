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
 
% Function for creation of report file with analysis results.

function outputReport( outputDir, problemName )

% replace underscores in problemName
problemNameWithoutUnderscores                                      = problemName ;
problemNameWithoutUnderscores( find( problemNameWithoutUnderscores == '_') ) = [] ;


% -------------------------------------------------------------------------
% ------------------- generation of the report  main tex file -------------
fprintf(  '|  - writing report file ... ')
fileReport = fopen( [ outputDir  problemName '_Report.tex' ] ,'w') ;

fprintf(fileReport, [ '\\documentclass[a4paper,10pt]{article} \n'] ) ;
fprintf(fileReport, [ '\\usepackage[a4paper,margin=20mm]{geometry} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{longtable} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{float} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{adjustbox} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{graphicx} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{color} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{booktabs} \n'] ) ;
fprintf(fileReport, [ '\\aboverulesep=0ex \n'] ) ;
fprintf(fileReport, [ '\\belowrulesep=0ex \n'] ) ;
fprintf(fileReport, [ '\\usepackage{array} \n'] ) ;
fprintf(fileReport, [ '\\newcolumntype{?}{!{\\vrule width 1pt}}\n' ] ) ;
fprintf(fileReport, [ '\\usepackage{fancyhdr} \n'] ) ;
fprintf(fileReport, [ '\\pagestyle{fancy} \n' ] ) ;
fprintf(fileReport, [ '\\fancyhf{} \n' ] ) ;
%~ fprintf(fileReport, [ '\\lhead{Report ONSAS version %s} \n'], ONSASversion );
fprintf(fileReport, [ '\\rhead{Problem: %s } \n'], problemNameWithoutUnderscores );
fprintf(fileReport, [ '\\lfoot{Date: \\today} \n'] );
fprintf(fileReport, [ '\\rfoot{Page \\thepage} \n'] );
fprintf(fileReport, [ '\\renewcommand{\\footrulewidth}{1pt} \n'] ) ;
fprintf(fileReport, [ '\\renewcommand{\\headrulewidth}{1.5pt} \n'] ) ;
fprintf(fileReport, [ '\\setlength{\\parindent}{0pt} \n'] ) ;
fprintf(fileReport, [ '\\usepackage[T1]{fontenc} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{libertine} \n'] ) ;
fprintf(fileReport, [ '\\usepackage{arydshln} \n'] ) ;
fprintf(fileReport, [ '\\definecolor{miblue}{rgb}{0,0.1,0.38} \n' ] ) ;
fprintf(fileReport, [ '\\usepackage{titlesec} \n'] ) ;
fprintf(fileReport, [ '\\titleformat{\\section}{\\normalfont\\Large\\color{miblue}\\bfseries}{\\color{miblue}\\sectionmark\\thesection}{0.5em}{}[{\\color{miblue}\\titlerule[0.5pt]}] \n\n' ] ) ;


fprintf(fileReport, [ '\\begin{document} \n'] );


%~ fprintf(fileReport, [ '\\begin{center} \n'] );
%~ fprintf(fileReport, [ '\\textbf{ONSAS v.' ONSASversion  ' analysis report \\\\ Problem: ' problemNameWithoutUnderscores '} \n'] );
%~ fprintf(fileReport, [ '\\end{center} \n\n'] );

fprintf(fileReport, [ 'This is an ONSAS automatically-generated report with part of the results obtained after the analysis. The user can access other magnitudes and results through the GNU-Octave/MATLAB console. The code is provided AS IS \\textbf{WITHOUT WARRANTY of any kind}, express or implied.\n\n'] );

% ==============================================================================
% ----------------------------    Problem Data    ------------------------------

%~ fprintf(fileReport, [ '\\section{Problem data} \n\n'] ) ;
%~ fprintf(fileReport, [ 'This section contains the model data provided in the input. The information presented are the general mesh data and the constitutive and geometrical properties and the load cases. The nodes and conectivity matrix are presented in the \\textbf{Appendix} section. \n\n'] );


%~ fprintf(fileReport, [ '\\subsection{Mesh data} \n\n'] ) ;
%~ fprintf(fileReport, [ '\\subsubsection{Nodes data} \n\n'] ) ;

%~ % Nodes
%~ fprintf(fileReport, [ 'The number of nodes and fixed nodes are presented below. \n\n'] ) ;

%~ [enc,fin] = tablesFunc( 'Data & Number', 2, 'c|c', 'General data of nodes.' ) ;
%~ fprintf(fileReport, '%s', enc );
%~ fprintf(fileReport, [ 'Number of nodes & %i \\\\ \n' ], nNodes ) ;
%~ fprintf(fileReport, [ 'Number of nodes with disp. constraints & %i \\\\ \n' ], size(nodalSprings,1) ) ;
%~ fprintf(fileReport, [ 'Number of free nodes & %i \n' ], nNodes-size(nodalSprings,1) ) ;
%~ fprintf(fileReport, '%s', fin ) ;

%~ fprintf(fileReport, [ 'The fixed nodes are listed in the table below.' ] ) ;
%~ [enc,fin] = tablesFunc( 'Node & $u_x$ & $\theta_x$ & $u_y$ & $\theta_y$ & $u_z$ & $\theta_z$ ', 7, 'c|c|c|c|c|c|c', 'Fixed nodes data.' ) ;
%~ fprintf(fileReport, '%s', enc ) ;
%~ for i = 1:size(nodalSprings,1)
  %~ fprintf(fileReport, '%i', nodalSprings(i,1) ) ;
  %~ for j=2:7
    %~ fprintf(fileReport, ' & %i ', nodalSprings(i,j) ) ;
  %~ end
  %~ fprintf(fileReport, '\\\\ \n' ) ;
%~ end
%~ fprintf(fileReport, '%s', fin ) ;

%~ fprintf(fileReport, [ '\\subsubsection{Elements data} \n\n'] ) ;
%~ % Elements
%~ fprintf(fileReport, [ 'The number and type of elements are listed below. \n\n'] ) ;

%~ [enc,fin] = tablesFunc( 'Data & Number', 2, 'c|c', 'General data of elements.' ) ;

%~ ntruss = sum(Conec(:,7)==1) ;
%~ nbeam  = sum(Conec(:,7)==2) ;
%~ ntet   = sum(Conec(:,7)==3) ;
%~ nplate = sum(Conec(:,7)==4) ;

%~ fprintf(fileReport, '%s', enc );
%~ fprintf(fileReport, [ 'Number of elements & %i \\\\ \n' ], nElems ) ;
%~ fprintf(fileReport, [ 'Number of truss elements & %i \\\\ \n' ], ntruss ) ;
%~ fprintf(fileReport, [ 'Number of beam elements & %i \\\\ \n' ], nbeam ) ;
%~ fprintf(fileReport, [ 'Number of plate elements & %i \\\\ \n' ], nplate ) ;
%~ fprintf(fileReport, [ 'Number of tetrahedron elements & %i \n' ], ntet ) ;
%~ fprintf(fileReport, '%s', fin ) ;

%~ fprintf(fileReport, [ '\\clearpage\n\n']) ;

% Materials
%~ formatMaterial = '%12.3e' ;

%~ fprintf(fileReport, [ '\\subsection{Material properties} \n\n'] ) ;
%~ fprintf(fileReport, [ 'The mechanical properties and the constitutive model of the defined materials are listed in the corresponding tables below. The number of defined materials is: %i. \\\\ \n\n'], nMats ) ;

%~ [enc, fin] = tablesFunc( 'Parameter & Value', 2, 'c|c', 'Constitutive parameters of the material.') ;
%~ vecMat = { '$E_t$'; '$E_c$'; '$\nu$'; '$G$'; '$\rho$' ; '$Prestrain$' } ;
%~ vecVal = [] ;
%~ for i = 1:nMats
	%~ Et 	= materialsParams{i}(2) ;
  %~ if materialsParams{i}(1) == 1
    %~ nu = materialsParams{i}(3) ; Ec = Et ; G = Et/(2*(1+nu)) ; Prestrain = 0 ;
  %~ elseif materialsParams{i}(1) == 2
    %~ Ec = materialsParams{i}(3) ; nu = 0 ; G = 0 ; Prestrain = 0 ;
  %~ else
    %~ Ec = 0 ; nu = 0 ; G = 0 ; Prestrain = 0 ;
  %~ end

  %~ vecVal = [ Et, Ec, nu, G, rho, Prestrain ] ;
  %~ fprintf(fileReport, [ '\\textbf{Material %i}: material parameters are presented in the table.\n'], i ) ;
  %~ fprintf(fileReport, '%s', enc ) ;
  %~ for j = 1:size(vecMat,1)
    %~ fprintf(fileReport, ['%s & ' formatMaterial ], vecMat{j}, vecVal(j) );
    %~ fprintf(fileReport, '\\\\ \n' ) ;
  %~ end
  %~ fprintf(fileReport, '%s', fin ) ;
%~ end

%~ % Sections
%~ numberSections  = size( crossSecsParams, 1 ) ;
%~ formatSections = '%12.3e' ;

%~ fprintf(fileReport, [ '\\subsection{Sections properties}\n\n'] ) ;
%~ fprintf(fileReport, [ 'The geometrical properties of the defined sections are listed below. The number of defined sections is: %i. \\\\ \n\n'], numberSections ) ;

%~ [enc, fin] = tablesFunc( 'Geometrical property & Value ', 3, 'c|c', 'Geometrical properties of the section.') ;
%~ vecSec = { '$A$' ; '$I_y$' ; '$I_z$' ; '$J$' } ;
%~ for i = 1:numberSections
  %~ fprintf(fileReport, [ '\\textbf{Section %i}: geometry parameters are presented in the table.\n'], i ) ;
  %~ fprintf(fileReport, '%s', enc ) ;
  %~ for j = 1:size(vecSec,1)
    %~ fprintf(fileReport, ['%s & ' formatSections ], vecSec{j}, crossSecsParams(i,j) );
  %~ fprintf(fileReport, '\\\\ \n' ) ;
  %~ end
  %~ fprintf(fileReport, '%s', fin ) ;
%~ end

%~ % Load Cases
%~ fprintf(fileReport, [ '\\subsection{Load Cases}\n\n'] ) ;
%~ fprintf(fileReport, [ 'The load cases considered in the problem are listed in the table below. \n\n'] ) ;

%~ if norm(constantFext) > 0
	%~ constant = 'yes' ;
%~ else
	%~ constant = 'no' ;
%~ end

%~ selfWeightText = 'no' ;

%~ if norm(variableFext) > 0
	%~ variable = 'yes' ;
		%~ timesVec = 0:deltaT:finalTime ;
		%~ lambdamax = max(loadFactorsFunc(timesVec)) ; lambdamin = min(loadFactorsFunc(timesVec)) ; lambda0 = loadFactorsFunc(0) ;
%~ else
	%~ variable = 'no' ; lambdamax = 0 ; lambdamin = 0 ; lambda0 = 0 ;
%~ end

%~ if length(userLoadsFilename) > 0
	%~ userLoads = 'yes' ;
%~ else
	%~ userLoads = 'no' ;
%~ end

%~ fprintf(fileReport, [ '\\begin{table}[!htb]\n\\centering\n\\begin{tabular}{c|c}\n Load Case & Status \\\\ \\toprule \n' ] ) ;

%~ fprintf(fileReport, [ 'Self Weight & %s \\\\ \n' ], selfWeightText ) ;
%~ fprintf(fileReport, [ 'Constant Loads & %s \\\\ \n' ], constant ) ;
%~ fprintf(fileReport, [ 'Variable Loads & %s \\\\ \n' ], variable ) ;
%~ fprintf(fileReport, [ '$\\lambda_{0}$ & %12.3e \\\\ \n' ], lambda0 ) ;
%~ fprintf(fileReport, [ '$\\lambda_{min}$ & %12.3e \\\\ \n' ], lambdamin ) ;
%~ fprintf(fileReport, [ '$\\lambda_{max}$ & %12.3e \\\\ \n' ], lambdamax ) ;
%~ fprintf(fileReport, [ 'User Loads Func & %s\\\\ \n' ], userLoads ) ;
%~ fprintf(fileReport, [ '\\end{tabular}\n\\caption{Load cases definition.}\n\\end{table}\n' ] ) ;

%~ fprintf(fileReport, [ '\\newpage\n\n' ] ) ;
% ==============================================================================
% --------------------------    Analysis Output    -----------------------------

fprintf(fileReport, [ '\\section{Analysis results}\n\n'] ) ;

%~ fprintf(fileReport, [ '\\subsection{General parameters}\n\n'] ) ;

% Prints numerical methods and analysis parameters

%~ if length(numericalMethodParams) > 0
  %~ %
  %~ if numericalMethodParams(1) == 0
    %~ vecMethod = { 'Numerical method' ; '$\delta_t$' ; '$t_f$'} ;
    %~ vecParams = [ {'Linear variable force'} ; num2cell(numericalMethodParams(2:end-1)') ] ;
  %~ %
  %~ elseif numericalMethodParams(1) == 1
    %~ vecMethod = { 'Numerical method' ; 'Tol $\Delta_u$' ; 'Tol forces' ; 'Tol iterations' ; 'Target $\lambda(t)$' ; 'Load steps' } ;
    %~ vecParams = [ {'Newton Raphson'} ; num2cell(numericalMethodParams(2:end)') ] ;
  %~ %
  %~ elseif numericalMethodParams(1) == 2
    %~ vecMethod = { 'Numerical method' ; 'Tol $\Delta_u$' ; 'Tol forces' ; 'Tol iterations' ; 'Target $\lambda(t)$' ; 'Load steps' ; 'Arc Length increment'} ;
    %~ vecParams = [ {'Newton Raphson - Arc Length'} ; num2cell(numericalMethodParams(2:end)') ] ;
  %~ %
  %~ elseif numericalMethodParams(1) == 3
    %~ vecMethod = { 'Numerical method' ; '$\Delta_t$' ; '$t_f$' ; 'Tol $\Delta_u$' ; 'Tol forces' ; 'Tol iterations' ; '$\delta_{NW}$' ; '$\alpha_{NW}$' } ;
    %~ vecParams = [ {'Newmark'} ; num2cell(numericalMethodParams(2:end)') ] ;
  %~ %
  %~ elseif numericalMethodParams(1) == 4
    %~ vecMethod = { 'Numerical method' ; '$\Delta_t$' ; '$t_f$' ; 'Tol $\Delta_u$' ; 'Tol forces' ; 'Tol iterations' ; '$\alpha_{NW}$' } ;
    %~ vecParams = [ {'HHT'} ; num2cell(numericalMethodParams(2:end)') ] ;
  %~ end
%~ end


%~ if exist('nonHomogeneousInitialCondU0') ~= 0
  %~ initialCondU0 = 'yes' ;
%~ else
  %~ initialCondU0 = 'no' ;
%~ end

%~ if exist('nonHomogeneousInitialCondUdot0') ~= 0
  %~ initialCondUdot0 = 'yes' ;
%~ else
  %~ initialCondUdot0 = 'no' ;
%~ end

%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  %~ eigValuesM = eig(Kliblib) ;
  %~ valProp = 0 ;
  %~ for i = 1:size(Kliblib,1)
    %~ if Kliblib(i,i) <= 0
      %~ valProp = valProp + 1 ;
    %~ end
  %~ end
%~ end
%~ if nonLinearAnalysisBoolean == 1 && dynamicAnalysisBoolean == 0
  %~ fprintf(fileReport, [ 'Non Linear Analysis is performed.\n' ] ) ;
%~ elseif nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 1
  %~ fprintf(fileReport, [ 'Dynamic Analysis is performed.\n' ] ) ;
%~ elseif nonLinearAnalysisBoolean == 1 && dynamicAnalysisBoolean == 1
  %~ fprintf(fileReport, [ 'Non Linear and Dynamic Analysis is performed.\n' ] ) ;
%~ else
  %~ fprintf(fileReport, [ 'Linear Analysis is performed.\\\\\n' ] ) ;
  %~ if LBAAnalyFlag == 1
    %~ fprintf(fileReport, [ 'Linear buckling analysis is performed.\n' ]) ;
  %~ end
  %~ fprintf(fileReport, [ 'Number of negative eigenvalues: %3i'], valProp ) ;
%~ end

%~ [enc1,fin1] = tablesFunc( 'Optional parameters & Status', 2, 'c|c', 'Analysis optional parameters status.' ) ;
%~ [enc2,~] = tablesFunc( 'Node & Dof & Value', 3, 'c|c|c', 'Optional parameters values.' ) ;
%~ if nonLinearAnalysisBoolean == 1 || dynamicAnalysisBoolean == 1
  %~ fprintf(fileReport, '%s', enc1 ) ;
  %~ fprintf(fileReport, [ 'Non homogeneous initial cond. $u_0$ & %s  \\\\ \n' ], initialCondU0 ) ;
  %~ if dynamicAnalysisBoolean == 1
    %~ fprintf(fileReport, [ 'Non homogeneous initial cond. $\\dot{u_0}$ & %s \\\\ \n' ], initialCondUdot0 ) ;
		%~ fprintf(fileReport, [ 'Nodal damping $C$ & %8.2e' ], nodalDamping ) ;
  %~ end
  %~ fprintf(fileReport, '%s', fin1 ) ;
	%~ if strcmp(initialCondU0,'yes')
		%~ fprintf(fileReport, '%s', enc2 ) ;
		%~ for i = 1:size(nonHomogeneousInitialCondU0,1)
			%~ fprintf(fileReport, [ '%i & %i & %i' ], [nonHomogeneousInitialCondU0(i,1) nonHomogeneousInitialCondU0(i,2) nonHomogeneousInitialCondU0(i,3)] ) ;
			%~ fprintf(fileReport, [ '\\\\ \n' ] ) ;
		%~ end
		%~ fprintf(fileReport, '\n\\end{tabular}\n\\end{adjustbox}\n\\caption{Non homogeneous initial condition $u_0$ values.}\n\\end{table} \n\n' ) ;
  %~ end
  %~ if strcmp(initialCondUdot0,'yes')
		%~ fprintf(fileReport, '%s', enc2 ) ;
		%~ for i = 1:size(nonHomogeneousInitialCondUdot0,1)
			%~ fprintf(fileReport, [ '%i & %i & %i' ], [nonHomogeneousInitialCondUdot0(i,1) nonHomogeneousInitialCondUdot0(i,2) nonHomogeneousInitialCondUdot0(i,3)] ) ;
			%~ fprintf(fileReport, [ '\\\\ \n' ] ) ;
		%~ end
		%~ fprintf(fileReport, '\n\\end{tabular}\n\\end{adjustbox}\n\\caption{Non homogeneous initial condition $\\dot{u_0}$ values.}\n\\end{table} \n\n' ) ;
  %~ end
%~ end

%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  %~ if length(prescribedDisps)>0
    %~ fprintf(fileReport, '%s', enc1 ) ;
    %~ fprintf(fileReport, [ 'Prescribed disp. & %s \\\\ \n' ], prescDisps ) ;
    %~ fprintf(fileReport, '%s', fin1 ) ;
    %~ fprintf(fileReport, '%s', enc2 ) ;
    %~ for i = 1:size(prescribedDisps,1)
      %~ fprintf(fileReport, [ '%i & %i & %i' ], [prescribedDisps(i,1) prescribedDisps(i,2) prescribedDisps(i,3)] ) ;
      %~ fprintf(fileReport, [ '\\\\ \n' ] ) ;
    %~ end
    %~ fprintf(fileReport, '\n\\end{tabular}\n\\end{adjustbox}\n\\caption{Prescribed displacements values.}\n\\end{table} \n\n' ) ;
  %~ end
%~ end

%~ if length(numericalMethodParams) > 0
  %~ [enc, fin] = tablesFunc( 'Analysis parameter & Value', 2, 'l|c', 'Analysis parameters considered and numerical method.') ;
  %~ fprintf(fileReport, '%s', enc ) ;
  %~ fprintf(fileReport, ['%s & %s \\\\ \n' ], vecMethod{1}, vecParams{1} );
    %~ for j = 2:size(vecMethod,1)
      %~ fprintf(fileReport, ['%s & %12.3e' ], vecMethod{j}, vecParams{j} );
      %~ fprintf(fileReport, '\\\\ \n' ) ;
    %~ end
  %~ fprintf(fileReport, '%s', fin ) ;
%~ end

% Time report
%~ fprintf(fileReport, [ '\\subsection{Time performance}\n\n'] ) ;
%~ %
%~ fprintf(fileReport, [ '\\textbf{Reading and variables definition/verification}\n'] ) ;
%~ [enc, fin] = tablesFunc( 'Task & Time (s)', 2, 'c|c', 'Reading and variables definition/verification time performance.') ;
%~ fprintf(fileReport, '%s', enc )
%~ fprintf(fileReport, [ 'Reading input file: & %5.3e \\\\ \n'], tReadingInput) ;
%~ fprintf(fileReport, [ 'Variables verification: & %5.3e \\\\ \n'], tVarVer) ;
%~ fprintf(fileReport, [ 'Input auxiliar definitions: & %5.3e \\\\ \n'], tInputAuxDefs) ;
%~ fprintf(fileReport, [ '\\midrule\n'])
%~ fprintf(fileReport, [ 'Total elapsed time in reading and verification: & %5.3e \\\\ \n'], tReadingInput+tVarVer+tInputAuxDefs) ;
%~ fprintf(fileReport, '%s', fin )
%



%~ fprintf(fileReport, [ '\\textbf{Analysis}\n'] ) ;
%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
	%~ [enc, fin] = tablesFunc( 'Task & Time (s)', 2, 'c|c', 'Analysis time spent.') ;
	%~ fprintf(fileReport, '%s', enc )
	%~ fprintf(fileReport, [ 'Geometry computation: & %5.3e \\\\ \n'], tGeomReading) ;
	%~ fprintf(fileReport, [ 'Stiffness matrix assembly: & %5.3e \\\\ \n'], tStiffMatrix) ;
	%~ fprintf(fileReport, [ 'Loads assembly: & %5.3e \\\\ \n'], tLoadsAssembly) ;
	%~ fprintf(fileReport, [ 'System resolution: & %5.3e \\\\ \n'], tSystemResolution) ;
	%~ fprintf(fileReport, [ 'Elems. disps and solic.: & %5.3e \\\\ \n'], tSolicDisps) ;
	%~ fprintf(fileReport, [ '\\midrule\n'])
	%~ fprintf(fileReport, [ 'Total elapsed time in analysis: & %5.3e \\\\ \n'], tGeomReading+tStiffMatrix+tLoadsAssembly+tSystemResolution+tSolicDisps) ;
	%~ fprintf(fileReport, '%s', fin )
%~ else
	%~ fprintf(fileReport, [ '\\clearpage\n\n' ] ) ;
  %~ fprintf(fileReport, [ '\\begin{longtable}{cccc} \n'] );
  %~ fprintf(fileReport, [ '\\input{' problemName '_timePerformanceOutput.tex' '} \n'] ) ;
  %~ fprintf(fileReport, [ '\\caption{Incremental analysis time performance.}\n\\end{longtable}\n'] ) ;
%~ end





%
%~ if plotParamsVector(1)>0
	%~ fprintf(fileReport, [ '\\textbf{Plots}\n'] ) ;
	%~ [enc, fin] = tablesFunc( 'Task & Time (s)', 2, 'c|c', 'Plots time spent.') ;
	%~ fprintf(fileReport, '%s', enc )
	%~ if plotParamsVector(1) < 3
		%~ fprintf(fileReport, [ 'Deformed shape: & %5.3f \\\\ \n'], tDefShape) ;
		%~ fprintf(fileReport, [ 'Normal force: & %5.3f \\\\ \n'], tNormalForce) ;
		%~ if nonLinearAnalysisBoolean == 1 || dynamicAnalysisBoolean == 1
			%~ fprintf(fileReport, [ 'Load factor vs control disp:: & %5.3f \\\\ \n'], tLoadFac) ;
		%~ end
	%~ else
		%~ fprintf(fileReport, [ 'VTK ConecNodes function: & %5.3f \\\\ \n'], tVtkConecNodes) ;
		%~ fprintf(fileReport, [ 'VTK writer: & %5.3f \\\\ \n'], tVtkWriter) ;
	%~ end
	%~ if length(loadFactors)>1
		%~ fprintf(fileReport, [ 'Load vs Disps.: & %5.3f \\\\ \n'], tLoadDisps) ;
	%~ end
	%~ fprintf(fileReport, [ '\\midrule\n'])
	%~ if plotParamsVector(1) < 3
		%~ fprintf(fileReport, [ 'Total elapsed time in octave plots: & %5.3f \\\\ \n'], tDefShape+tLoadFac+tNormalForce+tLoadDisps) ;
	%~ else
		%~ fprintf(fileReport, [ 'Total elapsed time in vtk plots: & %5.3f \\\\ \n'], tVtkConecNodes+tVtkWriter) ;
	%~ end
  %~ fprintf(fileReport,'%s',fin);
%~ end

%


% Incremental analysis output
%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  %~ if norm(variableFext) > 0
    % Salida de datos para caso con load cases en linear analysis
  %~ end
%~ else
  %~ fprintf(fileReport, [ '\\clearpage\n\n' ] ) ;

  fprintf(fileReport, [ '\\begin{longtable}{cccccc} \n'] );
  fprintf(fileReport, [ '\\input{' problemName '_iterations.tex' '} \n'] ) ;
%  fprintf(fileReport, [ '\\caption{Output of incremental analysis.}\n\\end{longtable}\n\n'] ) ;
  fprintf(fileReport, [ '\\end{longtable}\n\n'] ) ;
  %~ fprintf(fileReport, [ '\\newpage \n\n' ] ) ;
%~ end



% ==============================================================================
% --------------------------    Plots and Tables    ----------------------------

%fprintf(fileReport, [ '\\section{Plots and Tables}\n\n'] ) ;

%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
	%~ auxv = max ( abs( matUs((3:6:end),:) ) ) ;  auxw = max ( abs( matUs((5:6:end), :) ) ) ;  auxf = max ( abs( matFs(:,:) ) ) ; auxstrain = max ( abs( trussStrain(:,:) ) ) ;
  %~ auxTruss = max ( abs( trussStrain(:,:,:) ) ) ;
%~ end

%-------------------------------- Deformed Shape -------------------------------

%~ if size(matUs,2) == 1 && size(matNts,2) == 1
	%~ matUs = [zeros(size(matUs,1),1) matUs] ;
	%~ matNts = [zeros(size(matNts,1),1) matNts] ;
%~ end

%~ nTimesTotal = size( matUs, 2 ) ;

%~ if (length(plotParamsVector)>1)
  %~ timesPlotsVec = round( linspace(1, nTimesTotal, plotParamsVector(2) ) ) ;
%~ else
  %~ timesPlotsVec = 1: size(matUs,2) ;
%~ end


%~ if plotParamsVector(1) < 3 && ( printFlag == 1 || printFlag == 2 )
  %~ fprintf(fileReport, [ '\\subsection{Plots}\n\n'] ) ;
  %%%fprintf(fileReport, [ 'The deformed shape of the structure is plotted in Octave as shown below. \n'] ) ;
  %~ fprintf(fileReport, [ '\\subsubsection{Deformed Shape}\n\n'] ) ;

	%~ if printFlag > 0
		%~ fprintf(fileReport, [ '\\begin{figure}[!htb] \n  \\centering \n ' ] ) ;
		%~ if printFlag == 1
			%~ fprintf(fileReport, [ '\\resizebox{.65\\textwidth}{!}{\\input{' problemName '_deform_' sprintf('%04i', length( timesPlotsVec )) '.tex}} \n'  ] ) ;
		%~ elseif printFlag == 2
			%~ fprintf(fileReport, [ '\\includegraphics[width =0.65\\textwidth]{' problemName '_deform_' sprintf('%04i', length( timesPlotsVec )) '}\n']) ;
		%~ end





		%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
			%~ if norm(variableFext) == 0
				%~ fprintf(fileReport, [ '\\caption{Deformed structure. Scale factor: ' sprintf('%8.2e', linearDeformedScaleFactor) '.}\n'] )
			%~ else
				%~ fprintf(fileReport, [ '\\caption{Deformed structure at time index %i, time %12.3e and load factor %12.3e. Scale factor: ' sprintf('%8.2e', linearDeformedScaleFactor) '.}\n'], timesPlotsVec(end), currTime, currLoadFactor )
			%~ end
		%~ else



			%~ fprintf(fileReport, [ '\\caption{Deformed structure at time index %i, time %12.3e and load factor %12.3e.}\n'], timesPlotsVec(end), currTime, currLoadFactor )
		%end
		%~ fprintf(fileReport,  '\\end{figure} \n\n' )
	%~ end

  %~ fprintf(fileReport, [ '\\subsubsection{Normal Force}\n\n'] ) ;

  %~ if printFlag > 0
		%~ fprintf(fileReport, [ '\\begin{figure}[!htb] \n  \\centering \n ' ] ) ;
		%~ if printFlag == 1
			%~ fprintf(fileReport, [ '\\resizebox{.65\\textwidth}{!}{\\input{' problemName '_normalForce_' sprintf('%04i', length( timesPlotsVec )) '.tex}} \n'  ] ) ;
		%~ elseif printFlag == 2
			%~ fprintf(fileReport, [ '\\includegraphics[width =0.65\\textwidth]{' problemName '_normalForce_' sprintf('%04i', length( timesPlotsVec )) '}\n']) ;
		%~ end


   % if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
			%if norm(variableFext) == 0
			%	fprintf(fileReport, [ '\\caption{Normal Force.}\n'] )
			%else
				%fprintf(fileReport, [ '\\caption{Normal Force at time index %i, time %12.3e and load factor %12.3e.}\n'], timesPlotsVec(end), currTime, currLoadFactor )
			%end
		%else


    	%~ fprintf(fileReport, [ '\\caption{Normal Force at time index %i, time %12.3e and load factor %12.3e.}\n'], timesPlotsVec(end), currTime, currLoadFactor )

    %end

    %~ fprintf(fileReport,  '\\end{figure} \n\n' )
	%~ end
%~ end
%~ fprintf(fileReport, [ '\\clearpage\n\n' ] ) ;
%--------------------------------- Control node --------------------------------




%~ if nonLinearAnalysisBoolean == 1 || dynamicAnalysisBoolean == 1
  %~ if length(loadFactors) > 1
    %~ fprintf(fileReport, [ '\\subsection{Control Dof}\n\n' ] ) ;
    %~ fprintf(fileReport, [ 'The control node information and displacement are shown below. \n' ] ) ;
    %~ [enc,fin] = tablesFunc( 'Node & Dof & Scale factor', 3, 'c|c|c', 'Control node information for incremental analysis.') ;
    %~ fprintf(fileReport, '%s', enc) ;
    %%%fprintf(fileReport, '%i & %i & %i \n', controlDofInfo(1), controlDofInfo(2), controlDofInfo(3)) ;
    %~ fprintf(fileReport, '%s', fin) ;
    %~ fprintf(fileReport, [ '\\begin{figure}[!htb] \n  \\centering \n ' ] ) ;
    %~ if printFlag == 1
      %~ fprintf(fileReport, [ '\\resizebox{.65\\textwidth}{!}{\\input{' problemName '_loadDisp.tex}} \n ' ] ) ;
    %~ elseif printFlag == 2
      %~ fprintf(fileReport, [ '\\includegraphics[width =0.65\\textwidth]{' problemName '_loadDisp}\n']) ;
    %~ end
    %~ fprintf(fileReport, [ '\\caption{Load vs Displacement of control node.}\n' ] );
    %~ fprintf(fileReport,  '\\end{figure} \n\n' );
  %~ end
%~ end


%-------------------------------- Solicitations --------------------------------


%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0 % PARA QUE NO ENTRE EL NONLINEAR ------------> VER COMO SACAR DATOS DE AHI

  %~ if nplate == 0 && ntet == 0
    %~ fprintf(fileReport, [ '\\subsection{Solicitations}\n\n'] ) ;
    %~ fprintf(fileReport, [ 'The solicitations of the elements are presented below. If exists the load case \\textbf{Variable Loads}, the results shown correspond to the final load factor.\n'] ) ;
  %~ end

  %~ if ntruss > 0
    %~ fprintf(fileReport, [ '\\subsubsection{Truss elements}' ] ) ;
    %~ fprintf(fileReport, [ 'The solicitations for truss elements are presented below. Compression force is considered negative.' ] ) ;
    %~ [enc,fin] = tablesFunc('Elem. & $A$ & $\varepsilon$ & $N$', 4, 'c?ccc', 'Table of normal forces for truss elements.' ) ;
    %~ fprintf(fileReport, '%s', enc ) ;
    %~ for m = 1:nElems
      %~ i = indexesElems(m) ;
      %~ auxStr = sum(trussStrain(i,:,:)) ;
      %~ if abs(auxStr) <= 1e-10*auxstrain
        %~ auxStr = 0 ;
        %~ formatStrain = '%2i' ;
      %~ else
        %~ formatStrain = '%10.2e' ;
      %~ end
      %~ auxFor = sum(normalForce(i,:,:)) ;
      %~ if abs(auxFor) <= 1e-10*auxf
        %~ auxFor = 0 ;
        %~ formatForce = '%2i' ;
      %~ else
        %~ formatForce = '%10.2e' ;
      %~ end
      %~ fprintf(fileReport, ['%3i & %12.3e & ' formatStrain ' & ' formatForce ' \\\\ \n ' ], [i crossSecsParams(Conec(m,6), 1) auxStr auxFor ] ) ;
    %~ end
    %~ fprintf(fileReport, '%s', fin) ;
  %~ elseif nbeam > 0
    %~ fprintf(fileReport, [ '\\subsubsection{Beam elements}' ] ) ;
    %~ fprintf(fileReport, [ 'The solicitations for beam elements are presented below. ' ] ) ;

    %~ [~,fin] = tablesFunc('','','', 'Table of solicitations for each element in local coordinates. Node 1 and Node 2 correspond to the first and second node in the connectivity of each element, respectively.') ;
    %~ fprintf(fileReport, [ '\\begin{table}[!htb]\n \\centering\n \\begin{adjustbox}{max width=\\textwidth} \\begin{tabular}{c?cccccc?cccccc} \n'] )
      %~ fprintf(fileReport, [ ' & \\multicolumn{6}{c?}{Node 1} & \\multicolumn{6}{c}{Node 2} \\\\ \\cmidrule{2-13} \n'] )
      %~ fprintf(fileReport, [ 'Element & ' ...
        %~ ' $F_x$ & $M_x$ & $F_y$ & $M_y$ & $F_z$ & $M_z$ & ' ...
        %~ ' $F_x$ & $M_x$ & $F_y$ & $M_y$ & $F_z$ & $M_z$ \\\\ \\toprule \n'] )
      %~ for i=1:nElems
        %~ if Conec(i,7) == 2
          %~ fprintf(fileReport, '%3i', i );
          %~ for j=1:12
            %~ aux = sum( ElemsSolic(i,j,:) ) ;
            %~ if abs(aux) <= 1e-10*auxf
              %~ aux = 0 ;
              %~ formatSolic='%2i' ;
            %~ else
              %~ formatSolic='%10.2e' ;
            %~ end
            %~ fprintf( fileReport, [ ' & ' formatSolic ] , aux ) ;
          %~ end
          %~ fprintf(fileReport, '\\\\ \n' ) ;
        %~ end
      %~ end
      %~ fprintf(fileReport, '%s', fin) ;
  %~ end

%~ %---------------------------------- Reactions ----------------------------------

  %~ fprintf(fileReport, [ '\\subsection{Reactions}\n\n'] ) ;
  %~ fprintf(fileReport, [ 'The reactions on the supports are presented below.\n'] ) ;


  %~ [enc,fin] = tablesFunc( 'Node & $F_x$ & $M_x$ & $F_y$ & $M_y$ & $F_z$ & $M_z$', 7, 'c?cccccc', 'Table of joint reactions in global coordinates.') ;
  %~ fprintf(fileReport, '%s ', enc ) ;
  %~ for i = 1:size(nodalSprings,1)
    %~ fprintf(fileReport, '%3i', i );
    %~ for j = 1:6
			%~ gdlnode = nodes2dofs(i,6) ;
			%~ aux = sum( Reactions(gdlnode(j),:) ) ;
      %~ if abs(aux)<= 1e-10*auxf
        %~ aux = 0 ;
        %~ formatReaction = '%2i' ;
      %~ else
        %~ formatReaction = '%10.2e' ;
      %~ end
      %~ fprintf( fileReport, [ ' & ' formatReaction ] , aux ) ;
    %~ end
    %~ fprintf(fileReport, '\\\\ \n' );
  %~ end
  %~ fprintf(fileReport, '%s ', fin ) ;

%~ %------------------------------ Nodal displacements ----------------------------

  %~ if nplate == 0 && ntet == 0
    %~ fprintf(fileReport, [ '\\subsection{Displacements}\n\n'] ) ;
    %~ fprintf(fileReport, [ 'The nodal displacements of the elements are presented below. If exists the load case \\textbf{Variable Loads}, the results shown correspond to the final load factor.\n'] ) ;
  %~ end
  %~ if nbeam > 0
    %~ fprintf(fileReport, [ '\\subsubsection{Beam elements} \n' ] ) ;

    %~ fprintf(fileReport, [ '\\begin{table}[!htb]\n \\centering\n \\begin{adjustbox}{max width=\\textwidth} \\begin{tabular}{c|cccccc|cccccc} \n'] )
    %~ fprintf(fileReport, [ ' & \\multicolumn{6}{c|}{Node 1} & \\multicolumn{6}{c}{Node 2} \\\\ \\cmidrule{2-13} \n'] )
    %~ fprintf(fileReport, [ 'Elem. & ' ...
        %~ ' $u_x$ & $\\theta_x$ & $u_y$ & $\\theta_y$ & $u_z$ & $\\theta_z$ & ' ...
        %~ ' $u_x$ & $\\theta_x$ & $u_y$ & $\\theta_y$ & $u_z$ & $\\theta_z$ \\\\ \\toprule \n'] )
    %~ for i=1:nElems
      %~ if Conec(i,7) == 2
        %~ fprintf(fileReport, '%3i', i );
        %~ for j=1:12
          %~ aux = sum( dispsElemsMat(i,j,:) ) ;
          %~ if abs(aux) <= 1e-10*(norm(matUs(:,end),1)/length(matUs(:,end)))
            %~ aux = 0 ;
            %~ formatDisps='%2i' ;
          %~ else
            %~ formatDisps='%8.2e' ;
          %~ end
          %~ fprintf( fileReport, [ ' & ' formatDisps ] , aux ) ;
        %~ end
        %~ fprintf(fileReport, '\\\\ \n' );
      %~ end
    %~ end
    %~ fprintf(fileReport, [ '\\end{tabular} \\end{adjustbox} \\caption{Table of displacements of both nodes of each element in global coordinates. Node 1 and Node 2 correspond to the first and second node in the connectivity of each element, respectively.} \\end{table} \n\n'] )

  %~ elseif ntruss > 0
    %~ fprintf(fileReport, [ '\\subsubsection{Truss elements} \n' ] ) ;
    %~ [enc,fin] = tablesFunc( 'Elem. & $u_{x,1}$ & $u_{x,2} $', 3, 'c?cc', 'Table of nodal displacements in local coordinates.') ;
    %~ fprintf(fileReport, '%s', enc ) ;
    %~ for i = 1:nElems
      %~ m = indexesElems(i) ;
      %~ auxDis1 = sum(trussDisps(m,1,:)) ;
      %~ if abs(auxDis1) <= 1e-10*auxTruss
        %~ auxDis1 = 0 ;
        %~ formatD1 = '%2i' ;
      %~ else
        %~ formatD1 = '%8.2e' ;
      %~ end
      %~ auxDis2  = sum(trussDisps(m,2,:)) ;
      %~ if abs(auxDis2) <= 1e-10*auxTruss
        %~ auxDis2 = 0 ;
        %~ formatD2 = '%2i' ;
      %~ else
        %~ formatD2 = '%8.2e' ;
      %~ end
      %~ fprintf(fileReport, ['%3i & ' formatD1 ' & ' formatD2 ' \\\\ \n' ] , [ i auxDis1 auxDis2] ) ;
    %~ end
    %~ fprintf(fileReport, '%s', fin ) ;
  %~ end


  %~ if nbeam > 0
    %~ fprintf(fileReport, [ '\\subsubsection{Interpolated beam elements displacements}' ] ) ;
    %~ fprintf(fileReport, [ '\\input{' problemName '_tables.tex' '} \n'] ) ;
  %~ end

%~ else
	%~ if nplate == 0 && ntet == 0
    %~ fprintf(fileReport, [ '\\subsection{Solicitations}\n\n'] ) ;
    %~ fprintf(fileReport, [ 'The solicitations of the elements are presented below. \n'] ) ;
  %~ end

	%~ fprintf(fileReport, [ '\\clearpage\n\n' ] ) ;
  %~ fprintf(fileReport, [ '\\begin{longtable}{ccccc} \n'] )
  %~ fprintf(fileReport, [ '\\input{' problemName '_incrementsNormalForceOutput.tex' '} \n'] ) ;
  %~ fprintf(fileReport, [ '\\caption{Output of incremental Normal Forces analysis.}\n\\end{longtable}\n'] ) ;
%~ end % endif

%~ fprintf(fileReport, [ '\\newpage\n\n' ] ) ;
%~ % ==============================================================================
%~ % -------------------------------    Appendix    -------------------------------

%~ fprintf(fileReport, [ '\\section{Appendix}\n\n'] ) ;

%~ % Nodes matrix
%~ fprintf(fileReport, [ '\\subsection{Nodes matrix}\n\n'] ) ;
%~ fprintf(fileReport, [ 'The nodes matrix is presented in the table below.\n'] ) ;
%~ fprintf(fileReport, [ '\\begin{longtable}{c|c|c|c|c|c|c}\n\\centering\n' ] ) ;
%~ fprintf(fileReport, [ 'Node & X & Y & Z \\\\ \\toprule \n' ] ) ;
%~ for i=1:nNodes
  %~ fprintf(fileReport, ['%2i & %4i & %4i & %4i'], [i Nodes(i,1) Nodes(i,2) Nodes(i,3)] ) ;
  %~ fprintf(fileReport, ' \\\\ \n' ) ;
%~ end
%~ fprintf(fileReport, [ '\\caption{Nodes coordinate matrix.}\n\\end{longtable}\n\n' ] ) ;

%~ % Conectivity matrix
%~ fprintf(fileReport, [ '\\subsection{Conectivity matrix}\n\n' ] ) ;
%~ fprintf(fileReport, [ 'The conectivity matrix is presented in the table below. \n' ] ) ;
%~ fprintf(fileReport, [ '\\begin{longtable}{c|c|c|c|c|c|c}\n\\centering\n' ] ) ;
%~ fprintf(fileReport, ['Node 1 & Node 2  & Node 3 & Node 4 & Mat. n$^o$ & Sec. n$^o$ & Elem. type \\\\ \\toprule \n'] ) ;
%~ for i = 1:nElems
  %~ for j = 1:7
    %~ if j == 7
      %~ fprintf(fileReport, '%i' , Conec(i,j)) ;
    %~ else
      %~ fprintf(fileReport, '%i & ' , Conec(i,j)) ;
    %~ end
  %~ end
  %~ fprintf(fileReport, '\\\\ \n' );
%~ end
%~ fprintf(fileReport, [ '\\caption{Conectivity matrix.}\n\\end{longtable}\n\n' ] ) ;



fprintf(fileReport, [ '\\end{document}'] ) ;
fclose(fileReport);


% ==============================================================================
% -------------------------------    Tables    ---------------------------------

% ==============================================================================


% -------------------------------------------------------------------------------
%~ % ------ generation of latex tables if linear analysis is performed -------------
%~ if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  %~ if nbeam > 0
  %~ fileTables = fopen( [ outputdir  problemName '_tables.tex' ] ,'w') ;

  %~ auxv = max ( abs( matUs((3:6:end),:) ) ) ;  auxw = max ( abs( matUs((5:6:end), :) ) ) ;  auxf = max ( abs( matFs(:,:) ) ) ;

  %~ npartelem = 9 ;
  %~ v = zeros(npartelem,1) ;
  %~ w = zeros(npartelem,1) ;
  %~ Rb = eye(4) ;
  %~ Rb(2,2) = -1 ;
  %~ Rb(4,4) = -1 ;

  %~ LocBendXYdofs = [ 3  6           3+ndofpnode 6+ndofpnode] ;
  %~ LocBendXZdofs = [ 5  4           5+ndofpnode 4+ndofpnode] ;

  %~ espacio = ' ' ;
  %~ sauxv = '' ;
  %~ sauxw = '' ;
  %~ xaux  = '' ;
  %~ stab = '' ;
  %~ stabaux = '' ;

  %~ for i = 1:npartelem
    %~ if i == npartelem
      %~ vi = sprintf('$v_{%i,L}$ \\\\ \\toprule \n', i) ;
      %~ wi = sprintf('$w_{%i,L}$ \\\\ \\toprule \n', i) ;
      %~ xi = sprintf('$x_{%i,L}$ \\\\ \\toprule \n', i) ;
    %~ else
      %~ vi = sprintf('$v_{%i,L}$ &', i) ;
      %~ wi = sprintf('$w_{%i,L}$ &', i) ;
      %~ xi = sprintf('$x_{%i,L}$ &', i) ;
    %~ end
    %~ sauxv = [ sauxv, vi, espacio ] ;
    %~ sauxw = [ sauxw, wi, espacio ] ;
    %~ xaux  = [ xaux , xi, espacio ] ;

    %~ tab   = sprintf('%s', 'c') ;
    %~ stab  = [ stab, tab ] ;
  %~ end
  %~ stab2 = [ '{c?',stab,'}' ] ;

  %~ % x
  %~ fprintf(fileTables,  '\\begin{table}[!htb] \\centering \\begin{adjustbox}{max width=\\textwidth} \\begin{tabular}%s \n', stab2) ;
  %~ fprintf(fileTables,['Elem. &  %s'], xaux) ;
  %~ for i = 1:nElems
		%~ if Conec(i,7) == 2
			%~ m = indexesElems(i) ;
			%~ fprintf(fileTables, '%3i', i) ;

			%~ x = linspace( 0, ElemLengths(m), npartelem ) ;

			%~ for j=1:npartelem
				%~ formatDisps='%8.1e' ;
				%~ fprintf(fileTables, [ ' & ' formatDisps ] , x(j) ) ;
			%~ end
			%~ fprintf(fileTables, ' \\\\ \n' );
    %~ end
  %~ end
  %~ fprintf(fileTables,  '\\end{tabular} \\end{adjustbox} \\caption{Local coordinates of data points.} \\end{table} \n\n' )

  %~ % -----------------------------------

  %~ % ---  v ----
  %~ % prints table environment beginning
  %~ fprintf(fileTables,  '\\begin{table}[!htb] \n \\centering \n \\begin{adjustbox}{max width=\\textwidth} \n \\begin{tabular}%s \n', stab2) ;
  %~ % prints headers
  %~ fprintf(fileTables,['Elem. &  %s'], sauxv) ;

  %~ for i = 1:nElems
		%~ if Conec(i,7) == 2
			%~ m = indexesElems(i) ;
			%~ fprintf(fileTables, '%3i', i) ;

			%~ R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;

			%~ localUelem = R' * dispsElemsMat(m,:,end)' ;

			%~ x = linspace( 0, ElemLengths(m), npartelem ) ;

			%~ for j=1:npartelem
				%~ N = bendingInterFuns( x(j), ElemLengths(m), 0) ;

				%~ v(j) = N * localUelem(LocBendXYdofs) ;

				%~ if abs(v(j)) <= 1e-10*auxv
					%~ v(j) = 0 ;
					%~ formatDisps='%2i' ;
				%~ else
					%~ formatDisps='%10.2e' ;
				%~ end
				%~ fprintf(fileTables, [ ' & ' formatDisps ] , v(j) ) ;
			%~ end
			%~ fprintf(fileTables, ' \\\\ \n' );
    %~ end
  %~ end
  %~ fprintf(fileTables,  '\\end{tabular} \\end{adjustbox} \\caption{Table of $v$ local displacements.} \\end{table} \n\n' )
  %~ % -----------------------------------

  %~ % w
  %~ fprintf(fileTables,  '\\begin{table}[!htb] \\centering \\begin{adjustbox}{max width=\\textwidth} \\begin{tabular}%s \n', stab2) ;
  %~ fprintf(fileTables,['Elem. &  %s'], sauxw) ;
  %~ for i = 1:nElems
		%~ if Conec(i,7) == 2
			%~ m = indexesElems(i) ;
			%~ fprintf(fileTables, '%3i', i) ;

			%~ R = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;

			%~ localUelem = R' * dispsElemsMat(m,:,end)' ;

			%~ x = linspace( 0, ElemLengths(m), npartelem ) ;

			%~ for j=1:npartelem
				%~ N = bendingInterFuns( x(j), ElemLengths(m), 0) ;

				%~ w(j) = N * Rb' * localUelem(LocBendXZdofs) ;

				%~ if abs(w(j)) <= 1e-10*auxw
					%~ w(j) = 0 ;
					%~ formatDisps='%2i' ;
				%~ else
					%~ formatDisps='%10.2e' ;
				%~ end
				%~ fprintf(fileTables, [ ' & ' formatDisps ] , w(j) ) ;
			%~ end
			%~ fprintf(fileTables, ' \\\\ \n' );
    %~ end
  %~ end
  %~ fprintf(fileTables,  '\\end{tabular} \\end{adjustbox} \\caption{Table of $w$ local displacements.} \\end{table} \n\n' )

  %~ fclose(fileTables) ;
  %~ % -----------------------------------
	%~ end % endif nbeam
%~ end

fprintf(  ' done.               |\n');
fprintf(    '|=================================================|\n' )
