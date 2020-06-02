% ==============================================================================
% --------     ONSAS: an Open Non-linear Structural Analysis System     --------
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

% ==============================================================================

% ==============================================================================
% -------------------------------    Main GUI    -------------------------------

function onsasGUI

  close all, clc
  
  addpath([pwd '/sources/'])
  addpath([pwd '/sourcesGUI/'])
  
  % Geometrical properties
  GUIgeometry
  
  ONSASversion = '0.1.8' ; 

  fig = figure('menubar', 'none', 'numbertitle', 'off','resize', 'off', 'position', vecFig, 'name', ['ONSAS: Open Non-linear Structural Analysis System. Version: ' sprintf('%s', ONSASversion)]) ;

  handlesFig = guidata(fig) ;
  % ------------------------------------------------------------------------------------------------------------------------
  
  % General handles
  handlesFig.nodes 	= [] 	;
  handlesFig.conec 	= [] 	;
  handlesFig.anaVec = [1] ;
  handlesFig.numVec = 0 	;
  % Structural properties handles
  handlesFig.matNumbers = '1' ;
  handlesFig.secNumbers = '1' ;
  handlesFig.matMatrix 	= zeros(str2double(handlesFig.matNumbers),5) ;
  handlesFig.secMatrix 	= zeros(str2double(handlesFig.secNumbers),4) ;
  % Boundary conditions handles
  handlesFig.loadRows 		= '1' ;
  handlesFig.suppsRows 		= '1' ;
  handlesFig.loadsMatrix 	= zeros(str2double(handlesFig.loadRows),7) 	;
  handlesFig.suppsMatrix 	= zeros(str2double(handlesFig.suppsRows),6) ;
  % Numerical method parameters handles
  handlesFig.prevValdeltaTLin 		= '0' ; 
  handlesFig.prevValfinalTimeLin 	= '0' ;
  handlesFig.prevValIndexTimeSol 	= '0' ;
  
  handlesFig.prevValTimeIncr 			= '0' ;
  handlesFig.prevValFinalTimeDyn 	= '0' ;
  handlesFig.prevValTolDeltaUDyn 	= '0' ;
  handlesFig.prevValTolForcesDyn 	= '0' ;
  handlesFig.prevValTolIterDyn 		= '0' ;
  handlesFig.prevValDeltaNW 			= '0' ;
  handlesFig.prevValAlphaNW 			= '0' ;
  
  handlesFig.prevValTolDeltaU = '0' ;
  handlesFig.prevValTolForces = '0' ;
  handlesFig.prevValTolIter 	= '0' ;
  handlesFig.prevValLoadFac 	= '0' ;
  handlesFig.prevValLoadSteps = '0' ;
  handlesFig.prevValIncremAL 	= '0' ;
  
  % Output handles
  handlesFig.reportBoolean 	= '0' ;
  handlesFig.plotOpt 				= '1' ;
  handlesFig.plotTimesNum 	= '1' ;
  handlesFig.printflagNum 	= '1' ;
  handlesFig.controldofVec 	= '1 1 -1' ;
  handlesFig.scaleNum 			= '1' ;
  handlesFig.plotViewVec 		= '1 1 1' ;

  % ------------------------------------------------------------------------------------------------------------------------ 
  
% => Input Panel
  handlesFig.inputPanel = uipanel ('parent', fig, "title", "Input options", 'titleposition', 'centertop', "position", vecInpPanel) ; 
  
% => Sub panel Analysis options
  analysisPanel
  
% => Action buttons Panel
  handlesFig.sysActions = uipanel 	(	'parent', fig, "position", vecSysActions) ;
  handlesFig.loadM 			= uicontrol (	"parent", handlesFig.sysActions, "string", "Load .m", ...
																			"position",[PPFromLeft PPFromBottomSys buttonSysWidth buttonSysHeight ], ...
																			'callback', {@loadFile, fig} ) ;
  handlesFig.saveM 			= uicontrol (	"parent", handlesFig.sysActions, "string", "Save .m", ...
																			"position",[PPFromLeft+incremSysHor PPFromBottomSys buttonSysWidth buttonSysHeight ], ...
																			'callback', {@saveFile, fig, ONSASversion} ) ;
  handlesFig.runONSAS 	= uicontrol (	"parent", handlesFig.sysActions, "string", "Run", ...
																			"position",[PPFromLeft+2*incremSysHor PPFromBottomSys buttonSysWidth buttonSysHeight], ...
																			'callback', {@runONSAS, fig} ) ; 
  handlesFig.GPLLicense = uicontrol (	"parent", handlesFig.sysActions, "string", "GPL License", ...
																			"position",[PPFromLeft+3*incremSysHor PPFromBottomSys buttonSysWidth buttonSysHeight], ...
																			'callback', {@license, fig} ) ; 
  handlesFig.quit 			= uicontrol (	"parent", handlesFig.sysActions, "string", "Quit", ...
																			"position",[PPFromLeft+4*incremSysHor PPFromBottomSys buttonSysWidth buttonSysHeight], ...
																			'callback', 'close (gcf)' ) ; 
  % ------------------------------------------------------------------------------------------------------------------------

% => Sub panel Geometry
  handlesFig.geomPanel = uipanel 	('parent', handlesFig.inputPanel, "position", vecGeomPanel) ;
  handlesFig.geomTitle = uicontrol('parent', handlesFig.geomPanel, 'style', 'text', 'string', 'Geometry options', 'HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft ref buttonWidth buttonHeight]) ; 
  % Geommetry buttons
  handlesFig.geomOpt = '' ;
  handlesFig.loadDxf = uicontrol('parent', handlesFig.geomPanel, 'string', 'Load .dxf', 'position', [PPFromLeft PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@loadFile, fig} ) ;
  
  handlesFig.loadMsh = uicontrol(	'parent', handlesFig.geomPanel, 'string', 'Load .msh', ...
																	'position', [PPFromLeft+incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ]	, 'callback', {@loadFile, fig} ) ;
  handlesFig.loadTxt = uicontrol(	'parent', handlesFig.geomPanel, 'string', 'Load .txt', ...
																	'position', [PPFromLeft+2*incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@loadFile, fig} ) ;
  handlesFig.view		 = uicontrol(	'parent', handlesFig.geomPanel, 'string', 'View', ...
																	'position', [PPFromLeft+3*incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@plot, fig} ) ;
  % ------------------------------------------------------------------------------------------------------------------------
  
% => Sub panel Structural Properties
  handlesFig.strucPanel = uipanel ('parent', handlesFig.inputPanel, "position", vecStrucPanel) ;
  handlesFig.strucTitle = uicontrol('parent', handlesFig.strucPanel, 'style', 'text', 'string', 'Structural Properties', 'HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft ref buttonWidth buttonHeight]) ; 
  % Structural buttons
  handlesFig.matsNumber 	= uicontrol('parent', handlesFig.strucPanel, 'style'	, 'edit', 'string'	, '# materials', 'tag', 'matNumbers', ...
																			'position', [PPFromLeft PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@editBoxes, fig} ) ;
  handlesFig.materialMat 	= uicontrol('parent', handlesFig.strucPanel, 'string'	, 'Material matrix'	, ...
																			'position', [PPFromLeft+incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@uitables, fig} ) ;
  handlesFig.secsNumber 	= uicontrol('parent', handlesFig.strucPanel, 'style'	, 'edit', 'string'	, '# sections', 'tag', 'secNumbers', ...
																			'position', [PPFromLeft+2*incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@editBoxes, fig} ) ;
  handlesFig.sectionMat 	= uicontrol('parent', handlesFig.strucPanel, 'string'	, 'Section matrix'	, ...
																			'position', [PPFromLeft+3*incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@uitables, fig} ) ;  
  % ------------------------------------------------------------------------------------------------------------------------

% => Sub panel Boundary conditions
  handlesFig.BCPanel = uipanel 	('parent', handlesFig.inputPanel, "position", vecBCPanel) ;
  handlesFig.BCTitle = uicontrol('parent', handlesFig.BCPanel, 'style', 'text', 'string', 'Boundary Conditions', 'HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft ref buttonWidth buttonHeight]) ; 
  % Geommetry buttons
  handlesFig.loadsMat 	= uicontrol('parent', handlesFig.BCPanel, 'string', 'Loads matrix', ...
																		'position', [PPFromLeft PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@uitables, fig} ) ;
  handlesFig.springsMat = uicontrol('parent', handlesFig.BCPanel, 'string', 'Springs matrix', ...
																		'position', [PPFromLeft+incremSysHor PPFromBottom buttonSysWidth buttonSysHeight ], 'callback', {@uitables, fig} ) ;
  % ------------------------------------------------------------------------------------------------------------------------
  
% => Sub panel output options
  handlesFig.outputPanel = uipanel 	('parent', handlesFig.inputPanel, "position", vecOutputPanel) ;
  handlesFig.outputTitle = uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Output Options', 'HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft 1.8*ref buttonWidth buttonHeight]) ; 
  % Outputbuttons
  % 1 row
  handlesFig.VTK 						= uicontrol('parent', handlesFig.outputPanel, 'style', 'radiobutton', 'string', 'VTK', ...
																				'position', [PPFromLeft PPFromBottom+buttonSysHeight buttonSysWidth buttonSysHeight ], ...
																				'callback', {@outputOpt, fig} ) ;
  handlesFig.Octave 				= uicontrol('parent', handlesFig.outputPanel, 'style', 'radiobutton', 'string', 'Octave', ...
																				'position', [PPFromLeft+editSep PPFromBottom+buttonSysHeight buttonSysWidth buttonSysHeight ], ...
																				'callback', {@outputOpt, fig} ) ;
  handlesFig.plotTimesText 	= uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Plot times:', ...
																				'position', [PPFromLeft+2*editSep PPFromBottom+buttonSysHeight buttonSysWidth buttonSysHeight ]) ;
  handlesFig.plotTimes 			= uicontrol('parent', handlesFig.outputPanel, 'style', 'edit', 'string', '1', 'tag', 'plotTimes', ...
																				'position', [PPFromLeft+3*editSep PPFromBottom+2+buttonSysHeight editW editH ], ...
																				'callback', {@editBoxes, fig}) ;
  handlesFig.report 				= uicontrol('parent', handlesFig.outputPanel, 'style', 'radiobutton', 'string', 'Report', 'tag', 'reportOpt', ...
																				'position', [PPFromLeft+4*editSep PPFromBottom+buttonSysHeight buttonSysWidth buttonSysHeight ], ...
																				'callback', {@outputOpt, fig} ) ;
  % 2 row
  handlesFig.printflagText 	= uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Print flag:', ...
																				'position', [PPFromLeft PPFromBottom buttonSysWidth buttonSysHeight ]) ;
  handlesFig.printflag 			= uicontrol('parent', handlesFig.outputPanel, 'style', 'edit', 'string', '1', 'tag', 'printflag', ...
																				'position', [PPFromLeft+editSep PPFromBottom+2 editW editH ], 'callback', {@editBoxes, fig}) ;
	% if nonLin
	handlesFig.controlDofText = uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Control dof:', 'position', [PPFromLeft+4*editSep PPFromBottom buttonSysWidth buttonSysHeight ], 'visible', 'off') ;
  handlesFig.controlDof 		= uicontrol('parent', handlesFig.outputPanel, 'style', 'edit', 'string', '1 1 -1', 'tag', 'controldof' ,'position', [PPFromLeft+5*editSep PPFromBottom+2 editW editH ], 'callback', {@editBoxes, fig}, 'visible', 'off') ;
	% if Lin
	handlesFig.scaleText 			= uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Scale:', ...
																				'position', [PPFromLeft+4*editSep PPFromBottom buttonSysWidth buttonSysHeight ], 'visible', 'off') ;
  handlesFig.scale 					= uicontrol('parent', handlesFig.outputPanel, 'style', 'edit', 'string', '1', 'tag', 'scale', ...
																				'position', [PPFromLeft+5*editSep PPFromBottom+2 editW editH ], 'callback', {@editBoxes, fig}, 'visible', 'off') ;
  
	handlesFig.plotViewText 	= uicontrol('parent', handlesFig.outputPanel, 'style', 'text', 'string', 'Plot view:', ...
																				'position', [PPFromLeft+2*editSep PPFromBottom buttonSysWidth buttonSysHeight ]) ;
  handlesFig.plotView 			= uicontrol('parent', handlesFig.outputPanel, 'style', 'edit', 'string', '1 1 1', 'tag', 'plotView', ...
																				'position', [PPFromLeft+3*editSep PPFromBottom+2 editW editH ], 'callback', {@editBoxes, fig}) ;
  % ------------------------------------------------------------------------------------------------------------------------
	
	% Close table button
	vecPos = get(handlesFig.anaPanel, 'position') ;
  set(handlesFig.geomPanel, 'units', 'pixels') 
  vecGeomPos = get(handlesFig.geomPanel, 'position') ;
  handlesFig.hideTable = uicontrol('parent', handlesFig.inputPanel, 'string', 'Close table', 'visible', 'off', 'position', [vecPos(3)+2*PPFromLeft+3*incremSysHor vecGeomPos(2)-5*vecGeomPos(4)+buttonSysHeight+PPFromBottom buttonSysWidth buttonSysHeight], 'tag', 'close table') ;

	% Tables handles
	
	% Column names
  colNamesMat 		= { 'Model', 'Es', 'Ec', 'nu', 'Prestrain' 						} ;
  colNamesSec 		= { 'A', 'Iy', 'Iz', 'It' 														} ;
  colNamesLoads 	= { 'Global flag', 'Fx', 'Mx', 'Fy', 'My', 'Fz', 'Mz' } ;
  colNamesSprings = { 'ux', 'thetax', 'uy', 'thetay', 'uz', 'thetaz' 		} ;
	
  handlesFig.matTable 	= uitable (	'parent', handlesFig.inputPanel, 'ColumnName', colNamesMat, 			'visible', 'off', ...
																		'position', [vecPos(3)+1.5 2 350 vecPos(4)-4.5*vecGeomPos(4)-1.], 'ColumnEditable', [true false false false false]) ;
	handlesFig.secTable 	= uitable (	'parent', handlesFig.inputPanel, 'ColumnName', colNamesSec, 			'visible', 'off', ...
																		'position', [vecPos(3)+1.5 2 350 vecPos(4)-4.5*vecGeomPos(4)-1], 'ColumnEditable', true ) ;
	handlesFig.loadsTable = uitable (	'parent', handlesFig.inputPanel, 'ColumnName', colNamesLoads, 		'visible', 'off', ...
																		'position', [vecPos(3)+1.5 2 350 vecPos(4)-4.5*vecGeomPos(4)-1], 'ColumnEditable', true ) ;
	handlesFig.suppsTable = uitable (	'parent', handlesFig.inputPanel, 'ColumnName', colNamesSprings, 	'visible', 'off', ...
																		'position', [vecPos(3)+1.5 2 350 vecPos(4)-4.5*vecGeomPos(4)-1], 'ColumnEditable', true ) ;
  
  
  % sets initial values
  set(handlesFig.matTable, 		'data', handlesFig.matMatrix) 	;
  set(handlesFig.secTable, 		'data', handlesFig.secMatrix) 	;
  set(handlesFig.loadsTable, 	'data', handlesFig.loadsMatrix) ;
  set(handlesFig.suppsTable, 	'data', handlesFig.suppsMatrix) ;
  

  fontsizeSet

  guidata(fig,handlesFig) ;
	% color gray 'backgroundcolor', [0.5 0.5 0.5] 
end


% ==============================================================================
% ----------------------------    Analysis Callbacks    ------------------------

function anaChoice(src, eventdata, fig, buttonHeight, buttonWidth)
  
  handles = guidata(fig) ;
  handles.buttonString = get(src, 'string') ;
  buttonString = handles.buttonString ;
  buttonVal = get(src, 'value') ;
  handles.anaVec = [] ;
  % Radiobtuttons of analysis options
  if strcmp(buttonString, 'Linear Analysis')
    set(handles.nonLinAna, 'value', 0), set(handles.dynNonLinAna, 'value', 0)
    handles.anaVec = [ 1 ] ;
    set(handles.scaleText			, 'visible', 'on')	, set(handles.scale			, 'visible', 'on'	)
    set(handles.controlDofText, 'visible', 'off')	, set(handles.controlDof, 'visible', 'off')
  elseif strcmp(buttonString, 'Dynamic NonLinear Analysis')
    set(handles.linAna, 'value', 0), set(handles.nonLinAna, 'value', 0)
    handles.anaVec = [ 2 ] ;
    set(handles.scaleText			, 'visible', 'off')	, set(handles.scale		 , 	'visible', 'off')
    set(handles.controlDofText, 'visible', 'on'	)	, set(handles.controlDof, 'visible', 'on'	)
  elseif strcmp(buttonString, 'NonLinear Analysis')  
    set(handles.linAna, 'value', 0) , set(handles.dynNonLinAna, 'value', 0) 
    handles.anaVec = [ 3 ] ;
    set(handles.scaleText			, 'visible', 'off')	, set(handles.scale			, 'visible', 'off')
    set(handles.controlDofText, 'visible', 'on')	, set(handles.controlDof, 'visible', 'on'	)
  end

  % Analysis numerical method choice
  if buttonVal == 1
    set(handles.numericalTitle, 'visible', 'on'), set(handles.numMethod, 'visible', 'on')
    if strcmp(buttonString, 'Linear Analysis')
      set(handles.numMethod, 'value', 1)
      set(handles.numMethod, 'string', {'Choose method', 'User Function'} )
      set(handles.tolDeltaUText			, 'visible', 'off'), set(handles.tolDeltaU			, 'visible', 'off') 
      set(handles.tolForcesText			, 'visible', 'off'), set(handles.tolForces			, 'visible', 'off') 
      set(handles.tolIterText				, 'visible', 'off'), set(handles.tolIter				, 'visible', 'off')
      set(handles.targetLoadFacText	, 'visible', 'off'), set(handles.targetLoadFac	, 'visible', 'off') 
      set(handles.loadStepsText			, 'visible', 'off'), set(handles.loadSteps			, 'visible', 'off')
      set(handles.incremALText			, 'visible', 'off'), set(handles.incremAL				, 'visible', 'off')
      
      set(handles.timeIncrText			, 'visible', 'off'), set(handles.timeIncr				, 'visible', 'off')
      set(handles.finalTimeText			, 'visible', 'off'), set(handles.finalTime			, 'visible', 'off')
      set(handles.tolDeltaUDynText	, 'visible', 'off'), set(handles.tolDeltaUDyn		, 'visible', 'off')
      set(handles.tolForcesDynText	, 'visible', 'off'), set(handles.tolForcesDyn		, 'visible', 'off')
      set(handles.tolIterDynText		, 'visible', 'off'), set(handles.tolIterDyn			, 'visible', 'off')
      set(handles.deltaNWText				, 'visible', 'off'), set(handles.deltaNW				, 'visible', 'off')
      set(handles.alphaNWText				, 'visible', 'off'), set(handles.alphaNW				, 'visible', 'off')
    elseif strcmp(buttonString, 'Dynamic NonLinear Analysis')
      set(handles.numMethod, 'value', 1)
      set(handles.numMethod, 'string', {'Choose method', 'Newmark'} )
      
      set(handles.deltaTText				, 'visible', 'off'), set(handles.deltaT					, 'visible', 'off')
      set(handles.finalTimeLinText	, 'visible', 'off'), set(handles.finalTimeLin		, 'visible', 'off')
      set(handles.indexTimeSolText	, 'visible', 'off'), set(handles.indexTimeSol		, 'visible', 'off')
      
      set(handles.tolDeltaUText			, 'visible', 'off'), set(handles.tolDeltaU			, 'visible', 'off') 
      set(handles.tolForcesText			, 'visible', 'off'), set(handles.tolForces			, 'visible', 'off') 
      set(handles.tolIterText				, 'visible', 'off'), set(handles.tolIter				, 'visible', 'off')
      set(handles.targetLoadFacText	, 'visible', 'off'), set(handles.targetLoadFac	, 'visible', 'off') 
      set(handles.loadStepsText			, 'visible', 'off'), set(handles.loadSteps			, 'visible', 'off')
      set(handles.incremALText			, 'visible', 'off'), set(handles.incremAL				, 'visible', 'off')
      
    elseif strcmp(buttonString, 'NonLinear Analysis')
      set(handles.numMethod, 'value', 1)
      set(handles.numMethod, 'string', {'Choose method', 'Newton Raphson', 'Newton Raphson-Arc Length'} )
    
      set(handles.deltaTText				, 'visible', 'off'), set(handles.deltaT					, 'visible', 'off')
      set(handles.finalTimeLinText	, 'visible', 'off'), set(handles.finalTimeLin		, 'visible', 'off')
      set(handles.indexTimeSolText	, 'visible', 'off'), set(handles.indexTimeSol		, 'visible', 'off')
      set(handles.timeIncrText			, 'visible', 'off'), set(handles.timeIncr				, 'visible', 'off')
      set(handles.finalTimeText			, 'visible', 'off'), set(handles.finalTime			, 'visible', 'off')
      set(handles.tolDeltaUDynText	, 'visible', 'off'), set(handles.tolDeltaUDyn		, 'visible', 'off')
      set(handles.tolForcesDynText	, 'visible', 'off'), set(handles.tolForcesDyn		, 'visible', 'off')
      set(handles.tolIterDynText		, 'visible', 'off'), set(handles.tolIterDyn			, 'visible', 'off')
      set(handles.deltaNWText				, 'visible', 'off'), set(handles.deltaNW				, 'visible', 'off')
      set(handles.alphaNWText				, 'visible', 'off'), set(handles.alphaNW				, 'visible', 'off')  
    end
  else
    set(handles.numMethod, 'value', 1)
    set(handles.numericalTitle		, 'visible', 'off'), set(handles.numMethod		, 'visible', 'off')
    
    set(handles.deltaTText				, 'visible', 'off'), set(handles.deltaT				, 'visible', 'off')
    set(handles.finalTimeLinText	, 'visible', 'off'), set(handles.finalTimeLin	, 'visible', 'off')
    set(handles.indexTimeSolText	, 'visible', 'off'), set(handles.indexTimeSol	, 'visible', 'off')
  
    set(handles.timeIncrText			, 'visible', 'off'), set(handles.timeIncr			, 'visible', 'off')
    set(handles.finalTimeText			, 'visible', 'off'), set(handles.finalTime		, 'visible', 'off')
    set(handles.tolDeltaUDynText	, 'visible', 'off'), set(handles.tolDeltaUDyn	, 'visible', 'off')
    set(handles.tolForcesDynText	, 'visible', 'off'), set(handles.tolForcesDyn	, 'visible', 'off')
    set(handles.tolIterDynText		, 'visible', 'off'), set(handles.tolIterDyn		, 'visible', 'off')
    set(handles.deltaNWText				, 'visible', 'off'), set(handles.deltaNW			, 'visible', 'off')
    set(handles.alphaNWText				, 'visible', 'off'), set(handles.alphaNW			, 'visible', 'off')
  
    set(handles.tolDeltaUText			, 'visible', 'off'), set(handles.tolDeltaU		, 'visible', 'off') 
    set(handles.tolForcesText			, 'visible', 'off'), set(handles.tolForces		, 'visible', 'off') 
    set(handles.tolIterText				, 'visible', 'off'), set(handles.tolIter			, 'visible', 'off')
    set(handles.targetLoadFacText	, 'visible', 'off'), set(handles.targetLoadFac, 'visible', 'off') 
    set(handles.loadStepsText			, 'visible', 'off'), set(handles.loadSteps		, 'visible', 'off')
    set(handles.incremALText			, 'visible', 'off'), set(handles.incremAL			, 'visible', 'off')
    
    set(handles.scaleText					, 'visible', 'off'), set(handles.scale				, 'visible', 'off')
    set(handles.controlDofText		, 'visible', 'off'), set(handles.controlDof		, 'visible', 'off')
  end
  guidata(fig, handles) ;
end  

% Callback editbox analysis numerical method params
function editBoxes(src, eventdata, fig)
  handles = guidata(fig);
  a = get(src, 'tag') ;
  inp = get(src, 'string') ;
  if strcmp(a, 'delta T')
    handles.prevValdeltaTLin 			= num2str(inp) ;
  elseif strcmp(a, 'finalTime')  
    handles.preValfinalTimeLin 		= num2str(inp) ;
  elseif strcmp(a, 'indexTimeSol')
    handles.prevValIndexTimeSol 	= num2str(inp) ;
  elseif strcmp(a, 'timeIncr')
    handles.prevValTimeIncr 			= num2str(inp) ;
  elseif strcmp(a, 'finalTimeDyn')
    handles.prevValFinalTimeDyn 	= num2str(inp) ;
  elseif strcmp(a, 'tolDeltaUDyn')
    handles.prevValTolDeltaUDyn 	= num2str(inp) ;
  elseif strcmp(a, 'tolForcesDyn')
    handles.prevValTolForcesDyn 	= num2str(inp) ;
  elseif strcmp(a, 'tolIterDyn')
    handles.prevValTolIterDyn 		= num2str(inp) ;
  elseif strcmp(a, 'deltaNW')
    handles.prevValDeltaNW 				= num2str(inp) ;
  elseif strcmp(a, 'alphaNW')
    handles.prevValAlphaNW 				= num2str(inp) ;
  elseif strcmp(a,'tolDeltaU')
    handles.prevValTolDeltaU 			= num2str(inp) ;
  elseif strcmp(a,'tolForces')
    handles.prevValTolForces 			= num2str(inp) ;
  elseif strcmp(a,'tolIter')
    handles.prevValTolIter 				= num2str(inp) ;
  elseif strcmp(a,'loadFac')
    handles.prevValLoadFac 				= num2str(inp) ;  
  elseif strcmp(a,'loadSteps')
    handles.prevValLoadSteps 			= num2str(inp) ;
  elseif strcmp(a, 'incremAL')  
    handles.prevValIncremAL 			= num2str(inp) ;
  elseif strcmp(a, 'matNumbers')
    handles.matNumbers 						= num2str(inp) ;
  elseif strcmp(a, 'secNumbers')
		handles.secNumbers 						= num2str(inp) ;
  elseif strcmp(a, 'plotTimes')
    handles.plotTimesNum 					= num2str(inp) ;
  elseif strcmp(a, 'printflag')
		handles.printflagNum 					= num2str(inp) ;
	elseif strcmp(a, 'controldof')
		handles.controldofVec 				= num2str(inp) ;
	elseif strcmp(a, 'plotView')
		handles.plotViewVec 					= num2str(inp) ;
	elseif strcmp(a, 'scale')
		handles.scaleNum							= num2str(inp) ;
  end  

  set(src,'string', inp) ;  
  guidata(fig, handles) ;
end  

function numMethodChoice(src, eventdata, fig)
  
  handles = guidata(fig) ;
  name = get(src, 'string') ;
  val = get(src, 'value') ;

  if strcmp(handles.buttonString, 'Linear Analysis')
    if val == 2
      set(handles.deltaTText			, 'visible', 'on'), set(handles.deltaT				, 'visible', 'on')
      set(handles.finalTimeLinText, 'visible', 'on'), set(handles.finalTimeLin	, 'visible', 'on')
      set(handles.indexTimeSolText, 'visible', 'on'), set(handles.indexTimeSol	, 'visible', 'on')
      handles.numVec = 1 ;
    else
      set(handles.deltaTText			, 'visible', 'off'), set(handles.deltaT				, 'visible', 'off')
      set(handles.finalTimeLinText, 'visible', 'off'), set(handles.finalTimeLin	, 'visible', 'off')
      set(handles.indexTimeSolText, 'visible', 'off'), set(handles.indexTimeSol	, 'visible', 'off')
      handles.numVec = 0 ;
    end
  elseif strcmp(handles.buttonString, 'Dynamic NonLinear Analysis')
    if val == 2
      set(handles.timeIncrText			, 'visible', 'on'), set(handles.timeIncr		, 'visible', 'on')
      set(handles.finalTimeText			, 'visible', 'on'), set(handles.finalTime		, 'visible', 'on')
      set(handles.tolDeltaUDynText	, 'visible', 'on'), set(handles.tolDeltaUDyn, 'visible', 'on')
      set(handles.tolForcesDynText	, 'visible', 'on'), set(handles.tolForcesDyn, 'visible', 'on')
      set(handles.tolIterDynText		, 'visible', 'on'), set(handles.tolIterDyn	, 'visible', 'on')
      set(handles.deltaNWText				, 'visible', 'on'), set(handles.deltaNW			, 'visible', 'on')
      set(handles.alphaNWText				, 'visible', 'on'), set(handles.alphaNW			, 'visible', 'on')
      handles.numVec = 1 ;
    else
      set(handles.timeIncrText			, 'visible', 'off'), set(handles.timeIncr			, 'visible', 'off')
      set(handles.finalTimeText			, 'visible', 'off'), set(handles.finalTime		, 'visible', 'off')
      set(handles.tolDeltaUDynText	, 'visible', 'off'), set(handles.tolDeltaUDyn	, 'visible', 'off')
      set(handles.tolForcesDynText	, 'visible', 'off'), set(handles.tolForcesDyn	, 'visible', 'off')
      set(handles.tolIterDynText		, 'visible', 'off'), set(handles.tolIterDyn		, 'visible', 'off')
      set(handles.deltaNWText				, 'visible', 'off'), set(handles.deltaNW			, 'visible', 'off')
      set(handles.alphaNWText				, 'visible', 'off'), set(handles.alphaNW			, 'visible', 'off')
      handles.numVec = 0 ;
    end
  elseif strcmp(handles.buttonString, 'NonLinear Analysis')
    if val == 2
      set(handles.tolDeltaUText			, 'visible', 'on')	, set(handles.tolDeltaU			, 'visible', 'on') 
      set(handles.tolForcesText			, 'visible', 'on')	, set(handles.tolForces			, 'visible', 'on') 
      set(handles.tolIterText				, 'visible', 'on')	, set(handles.tolIter				, 'visible', 'on')
      set(handles.targetLoadFacText	, 'visible', 'on')	, set(handles.targetLoadFac	, 'visible', 'on') 
      set(handles.loadStepsText			, 'visible', 'on')	, set(handles.loadSteps			, 'visible', 'on')
      set(handles.incremALText			, 'visible', 'off')	, set(handles.incremAL			, 'visible', 'off') 
      handles.numVec = 1 ;
    elseif val == 3
			set(handles.tolDeltaUText			, 'visible', 'on')	, set(handles.tolDeltaU			, 'visible', 'on') 
      set(handles.tolForcesText			, 'visible', 'on')	, set(handles.tolForces			, 'visible', 'on') 
      set(handles.tolIterText				, 'visible', 'on')	, set(handles.tolIter				, 'visible', 'on')
      set(handles.targetLoadFacText	, 'visible', 'on')	, set(handles.targetLoadFac	, 'visible', 'on') 
      set(handles.loadStepsText			, 'visible', 'on')	, set(handles.loadSteps			, 'visible', 'on')
      set(handles.incremALText			, 'visible', 'on')	, set(handles.incremAL			, 'visible', 'on')
      handles.numVec = 2 ;
    else
      set(handles.tolDeltaUText			, 'visible', 'off'), set(handles.tolDeltaU			, 'visible', 'off') 
      set(handles.tolForcesText			, 'visible', 'off'), set(handles.tolForces			, 'visible', 'off') 
      set(handles.tolIterText				, 'visible', 'off'), set(handles.tolIter				, 'visible', 'off')
      set(handles.targetLoadFacText	, 'visible', 'off'), set(handles.targetLoadFac	, 'visible', 'off') 
      set(handles.loadStepsText			, 'visible', 'off'), set(handles.loadSteps			, 'visible', 'off')
      set(handles.incremALText			, 'visible', 'off'), set(handles.incremAL				, 'visible', 'off')
      handles.numVec = 0 ;
    end
  end  
  guidata(fig, handles);
end  

% ==============================================================================
% ----------------------------   Tables callback   -----------------------------

function uitables(src, eventdata, fig)
  
  
  handles = guidata(fig) ; 
  
  % Tag
  a = get(src, 'string') ;
  set(handles.hideTable, 'callback', {@closeTable, fig, a})
	% Hide table button
  set(handles.hideTable, 'visible', 'on')
  % Matrix data
  if strcmp(a, 'Material matrix')
		% Material mat
		handles.matMatrixAux = zeros(str2double(handles.matNumbers),5) ;
		[m,n] = size(handles.matMatrix);
		handles.matMatrixAux(1:m,:) = handles.matMatrix ;
		set(handles.matTable, 'data', handles.matMatrixAux)
		set(handles.matTable, 'visible', 'on')
		
		set(handles.matTable, 'celleditcallback', {@uiTableEdit, fig, a})
		set(handles.secTable	, 'visible', 'off') 
		set(handles.loadsTable, 'visible', 'off')
		set(handles.suppsTable, 'visible', 'off')  
  elseif strcmp(a, 'Section matrix')
		% Section mat
		handles.secMatrixAux = zeros(str2double(handles.secNumbers),4) ;
		[m,n] = size(handles.secMatrix);
		handles.secMatrixAux(1:m,:) = handles.secMatrix ;
		set(handles.secTable, 'data', handles.secMatrixAux)
		set(handles.secTable, 'visible', 'on') 
		set(handles.secTable, 'celleditcallback', {@uiTableEdit, fig, a})
		set(handles.matTable	, 'visible', 'off') 
		set(handles.loadsTable, 'visible', 'off')
		set(handles.suppsTable, 'visible', 'off')
  elseif strcmp(a, 'Loads matrix')
		% Loads mat
		handles.loadsMatrixAux = zeros(handles.loadRows,7) ;
		[m,n] = size(handles.loadsMatrix);
		handles.loadsMatrixAux(1:m,:) = handles.loadsMatrix ;
		set(handles.loadsTable, 'data', handles.loadsMatrixAux)
		set(handles.loadsTable, 'visible', 'on')
		set(handles.loadsTable, 'celleditcallback', {@uiTableEdit, fig, a})
		set(handles.matTable	, 'visible', 'off') 
		set(handles.secTable	, 'visible', 'off')
		set(handles.suppsTable, 'visible', 'off')
  elseif strcmp(a, 'Springs matrix')
		% Springs mat
		handles.suppsMatrixAux = zeros(handles.suppsRows,6) ;
		[m,n] = size(handles.suppsMatrix);
		handles.suppsMatrixAux(1:m,:) = handles.suppsMatrix ;
		set(handles.suppsTable, 'data', handles.suppsMatrixAux)
		set(handles.suppsTable, 'visible', 'on')
		set(handles.suppsTable, 'celleditcallback', {@uiTableEdit, fig, a})
		set(handles.matTable	, 'visible', 'off') 
		set(handles.secTable	, 'visible', 'off')
		set(handles.loadsTable, 'visible', 'off')
  end
  % ----------
  guidata(fig, handles) ;
  
end  

function uiTableEdit(src, eventdata, fig, a)
  
  handles = guidata(fig) ;
  
  if strcmp(a, 'Material matrix')
		handles.matMatrix = get(handles.matTable, 'data') ;
		set(src, 'Data', handles.matMatrix)
		rows = size(handles.matMatrix,1) ;
		for i = 1:rows
			if handles.matMatrix(i,1) == 1
				set(handles.matTable, 'ColumnEditable', [true true false true false])
			elseif
				handles.matMatrix(i,1) == 2
				set(handles.matTable, 'ColumnEditable', [true true true false false])
			elseif
				handles.matMatrix(i,1) == 3
				set(handles.matTable, 'ColumnEditable', [true true true false true ])
			end 
		end 
  elseif strcmp(a, 'Section matrix')
		handles.secMatrix = get(handles.secTable, 'data') ;
		set(src, 'Data', handles.secMatrix) 
  elseif strcmp(a, 'Loads matrix')
		handles.loadsMatrix = get(handles.loadsTable, 'data') ;
		set(src, 'Data', handles.loadsMatrix)
  elseif strcmp(a, 'Springs matrix')
		handles.suppsMatrix = get(handles.suppsTable, 'data') ;
		set(src, 'Data', handles.suppsMatrix)
  end
  
  guidata(fig, handles) ;
  
end  

function closeTable(src, eventdata, fig, a)

	handles = guidata(fig) ; 
	
	set(handles.hideTable	, 'visible', 'off')
	set(handles.matTable	, 'visible', 'off')
	set(handles.secTable	, 'visible', 'off') 
	set(handles.loadsTable, 'visible', 'off')
	set(handles.suppsTable, 'visible', 'off')  
	
	guidata(fig, handles) ;
	
end

% ==============================================================================
% ----------------------------   Output Callbacks    ---------------------------

function outputOpt(src, eventdata, fig)

  handles = guidata(fig) ;
  
  buttonString = get(src, 'string') ;
  buttonVal = get(src, 'value') ;
  % Radiobuttons of output options
  if strcmp(buttonString, 'VTK')
    handles.plotOpt = '3' ;
    set(handles.Octave, 'value', 0) ;
  elseif strcmp(buttonString, 'Octave')
    handles.plotOpt = '2' ;
    set(handles.VTK, 'value', 0) ;
  else
    if get(handles.report, 'value')
      handles.reportBoolean = '1' ;
    else  
      handles.reportBoolean = '0' ;
    end
  end
  guidata(fig, handles) ;
end

% ==============================================================================
% ----------------------------    Actions Buttons    ---------------------------
 
%~ Load button 
function loadFile(src, eventdata, fig)

  handles = guidata(fig) ;
  buttonString = get(src, 'string') ;
  handles.geomOpt = buttonString ;
  
  if strcmp(buttonString, 'Load .dxf')
    [fname,fpath] = uigetfile('./input/*.dxf') ;
    if fname == 0
			return
    else
			[nodesMat, conecMat] = dxf2ONSAS(fname) ;
		end
  elseif strcmp(buttonString, 'Load .msh')
    [fname,fpath] = uigetfile('./input/*.msh') ;
    if fname == 0
			return
		else	
			[nodesMat, conecMat] = mshReader(fname) ;
		end	
  elseif strcmp(buttonString, 'Load .m')
    [fname,fpath] = uigetfile('./input/*.m') ;
    if fname == 0
			return
		else		
			acdir = pwd ; cd(fpath) ;
			% Variable to run ONSAS
			filesList = readdir('./') ;
			previouslyDefinedSelectedFileVar = 0 ;
			for i = 1:length(filesList)
				if strcmp(filesList(i,1),fname)  
					previouslyDefinedSelectedFileVar = i ;
				end  
			end
			handles.selectedFile = previouslyDefinedSelectedFileVar ;
			handles.fileName = fname ;
			run(fname);
			cd(acdir) ;
			nodesMat = [] ;
			conecMat = [] ;
			
      

			%~ loadMFileSets
			
		end	
    % ----------
  elseif strcmp(buttonString, 'Load .txt' )
    [fname,fpath] = uigetfile('./*.txt') ;
    if fname == 0
			return
    else
			[nodesMat, conecMat] = txt2ONSAS(fname) ; 
    end
  end
  handles.geomFileName = fname ;
  handles.geomFilePath = fpath ;
  handles.nodesMat = nodesMat ;
  handles.conecMat = conecMat ;
  
  suppTag = [] ;
	loadTag = [] ;
	
	for i = 1:size(nodesMat,1)
		if nodesMat(i,5)  > 0 
			if ~ismember(nodesMat(i,5), suppTag)
				suppTag = [suppTag ; nodesMat(i,5)] ;
			end
		end
		if nodesMat(i,4) > 0
			if ~ismember(nodesMat(i,4), loadTag)
				loadTag = [loadTag ; nodesMat(i,4)] ;
			end
		end
	end

	for i = 1:size(conecMat,1)
		if conecMat(i,9)  > 0 
			if ~ismember(conecMat(i,9), suppTag)
				suppTag = [suppTag ; conecMat(i,9)] ;
			end
		end
		if conecMat(i,8) > 0
			if ~ismember(conecMat(i,8), loadTag)
				loadTag = [loadTag ; conecMat(i,8)] ;
			end
		end
	end

	handles.suppsRows = length(unique(suppTag)) ;
	handles.loadRows = length(unique(loadTag)) ;

  guidata(fig, handles) ;

end

%~ Save button
function saveFile(src, eventdata, fig, version)

  handles = guidata(fig) ;
  cd './input' ;
  
  fileName1 = inputdlg('File Name') ;
  fileName = [fileName1{:} '.m'] ;
  handles.fileName = fileName ;
  
  fileM = fopen(fileName, 'w+') ;
  % Version
  fprintf(fileM, 'inputONSASversion = ''%s'' ; \n', version) ;
  fprintf(fileM, 'problemName = ''%s'' ; \n\n', fileName1{:}) ;
  % Material
  fprintf(fileM, '%% Constitutive properties\n' );
  fprintf(fileM, ['hyperElasParams = cell(%i,1) ;\n'], str2num(handles.matNumbers) );
  matsData = get(handles.matTable, 'data') ;
  
  for i = 1:str2num(handles.matNumbers)
		if matsData(i,1) == 1
			fprintf(fileM, ['hyperElasParams{%i} = [ 1 %12.3e %12.3e ] ;\n'], i, matsData(i,2), matsData(i,4) );
    elseif matsData(i,1) == 2
			fprintf(fileM, ['hyperElasParams{%i} = [ 2 %12.3e %12.3e ] ;\n'], i, matsData(i,2), matsData(i,3) );
    elseif matsData(i,1) == 3
			fprintf(fileM, ['hyperElasParams{%i} = [ 3 %12.3e %12.3e %12.3e ] ;\n'], i, matsData(i,2), matsData(i,3), matsData(i,5) );
    end
  end
  
  % Section
  fprintf(fileM, '%% Geometrical properties\n' );
  fprintf(fileM, ['secGeomProps = [ '], handles.secNumbers );
  secsData = get(handles.secTable, 'data') ;
  
  for i = 1:str2num(handles.secNumbers)
    fprintf(fileM, ['%12.3e %12.3e %12.3e %12.3e'], secsData(i,1), secsData(i,2), secsData(i,3), secsData(i,4) );
		if i == str2num(handles.secNumbers)
			fprintf(fileM, ' ] ;') ;
		else
			fprintf(fileM, ' ; ... \n') ;
		end
		fprintf(fileM, '\n') ;	
  end
  
  % Matrices
  auxSup = get(handles.suppsTable, 'data') ;
  auxLoad = get(handles.loadsTable, 'data') ;

	[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, ~, ~, nodalSprings] = inputFormatConversion ( handles.nodesMat, handles.conecMat, auxLoad, auxSup ) ;

  % Nodes
  fprintf(fileM, '%% Nodes\n' );
  nnodes = size(Nodes,1) ;
  fprintf(fileM, 'Nodes = [ ') ;
  for i = 1:nnodes
    fprintf(fileM, ['%4i %4i %4i'], [Nodes(i,1) Nodes(i,2) Nodes(i,3)] ) ;
    if i == nnodes
      fprintf(fileM, ' ] ;') ;
    else  
      fprintf(fileM, ' ; ...') ; 
    end
    fprintf(fileM, '\n' ) ;
  end

  % Conec
  fprintf(fileM, '%% Conec\n' );
  nelems = size(Conec,1) ;
  fprintf(fileM, 'Conec = [ ') ;
  for i = 1:nelems
    for j = 1:7
      fprintf(fileM, ['%4i %4i %4i %4i %4i %4i %4i %4i %4i'], [Conec(i,j)] ) ;
    end
    if i == nelems
      fprintf(fileM, ' ] ;') ;
    else  
      fprintf(fileM, ' ; ...')  
    end
    fprintf(fileM, '\n' ) ;
  end
  % Springs matrix
	fprintf(fileM, 'nodalSprings = [ ') ;
	for i = 1:size(nodalSprings,1)
		fprintf(fileM, '%i', nodalSprings(i,1) ) ;
		for j = 2:7
			fprintf(fileM, ' %3f ', nodalSprings(i,j) ) ;
		end
		if j == 7
			if i == size(nodalSprings,1)
				fprintf(fileM, ' ] ;') ;
			else 
				fprintf(fileM, ' ; ...') ;
			end	
		end
		fprintf(fileM, '\n') ;		
	end
	% Load matrix
	if size(nodalConstantLoads,1) == 0
		nodalConstBool = 0 ;
	else
		nodalConstBool = 1 ;
	end
	if size(nodalVariableLoads,1) == 0
		nodalVarBool = 0 ;
	else
		nodalVarBool = 1 ;
	end
	if nodalConstBool	
		fprintf(fileM, 'nodalConstantLoads = [ ');
		for i = 1:size(nodalConstantLoads,1)
			fprintf(fileM, '%i', nodalConstantLoads(i,1)) ;
			for j = 2:7
				fprintf(fileM, [' %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e '], nodalConstantLoads(i,j)) ;
			end
			if j == 7
				if i == size(nodalConstantLoads,1) 
					fprintf(fileM, ' ] ;') ;
				else
					fprintf(fileM, ' ; ...') ;
				end
			end
			fprintf(fileM, '\n') ;		
		end
	elseif nodalVarBool
		fprintf(fileM, 'nodalVariableLoads = [ ') ;
		for i = 1:size(nodalVariableLoads,1)
			fprintf(fileM, '%i', nodalVariableLoads(i,1)) ;
			for j = 2:7
				fprintf(fileM, [' %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e '], nodalVariableLoads(i,j)) ;
			end
			if j == 7
				if i == size(nodalVariableLoads,1) 
					fprintf(fileM, ' ] ;') ;
				else
					fprintf(fileM, ' ; ...') ;
				end
			end
			fprintf(fileM, '\n') ;		
		end
	end
	
  % Analysis options
  fprintf(fileM, '%% Analysis options\n' );
  if handles.anaVec == 1
    fprintf(fileM, ['nonLinearAnalysisBoolean = 0 ;\n'] );
    fprintf(fileM, ['dynamicAnalysisBoolean = 0 ;\n'] );
    if handles.numVec == 1
      fprintf(fileM, ['numericalMethodParams = [ %i %4i %4i %4i] ;\n'], [0 str2num(handles.prevValdeltaTLin) str2num(handles.prevValfinalTimeLin) str2num(handles.prevValIndexTimeSol)] ); 
    end 
  elseif handles.anaVec == 2
    fprintf(fileM, ['nonLinearAnalysisBoolean = 1 ;\n'] );
    fprintf(fileM, ['dynamicAnalysisBoolean = 1 ;\n'] );
    if handles.numVec == 1
      fprintf(fileM, ['numericalMethodParams = [ %i %12.3e %12.3e %12.3e %12.3e %12.3 %12.3e %12.3e ] ;\n'], [3  str2num(handles.prevValTimeIncr) str2num(handles.prevValFinalTimeDyn) str2num(handles.prevValTolDeltaUDyn) str2num(handles.prevValTolForcesDyn) str2num(handles.tolIterDyn) str2num(handles.prevValDeltaNW) str2num(handles.prevValAlphaNW) ] ); 
    end
  elseif handles.anaVec == 3
    fprintf(fileM, ['nonLinearAnalysisBoolean = 1 ;\n'] );
    fprintf(fileM, ['dynamicAnalysisBoolean = 0 ;\n'] );
    if handles.numVec == 1
			fprintf(fileM, [ 'stopTolDeltau = %12.3e ;\n' ], str2num(handles.prevValTolDeltaU)) ;
			fprintf(fileM, [ 'stopTolForces = %12.3e ;\n' ], str2num(handles.prevValTolForces)) ;
			fprintf(fileM, [ 'stopTolIts = %12.3e ;\n' ], str2num(handles.prevValTolIter)) ;
			fprintf(fileM, [ 'targetLoadFactr = %12.3e ;\n' ], str2num(handles.prevValLoadFac)) ;
			fprintf(fileM, [ 'nLoadSteps = %12.3e ;\n' ], str2num(handles.prevValLoadSteps)) ;
      fprintf(fileM, [ 'numericalMethodParams = [ %i %12.3e %12.3e %12.3e %12.3e %12.3e ] ;\n'], [1 str2num(handles.prevValTolDeltaU) str2num(handles.prevValTolForces) str2num(handles.prevValTolIter) str2num(handles.prevValLoadFac) str2num(handles.prevValLoadSteps)] ); 
    elseif handles.numVec == 2
      fprintf(fileM, [ 'numericalMethodParams = [ %i %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e ] ;\n'], [2 str2num(handles.prevValTolDeltaU) str2num(handles.prevValTolForces) str2num(handles.prevValTolIter) str2num(handles.prevValLoadFac) str2num(handles.prevValLoadSteps) str2num(handles.prevValIncremAL)] );
    end
  end
  % Output options
  fprintf(fileM, '%% Output options\n' );
  fprintf(fileM, ['plotParamsVector = [%i %3i] ;\n'], [ str2num(handles.plotOpt) str2num(handles.plotTimesNum) ] );
  fprintf(fileM, ['reportBoolean = %i ;\n'], str2num(handles.reportBoolean) );
  fprintf(fileM, ['printflag = %i ;\n'], str2num(handles.printflagNum)) ;
  if handles.anaVec == 1
		fprintf(fileM, ['linearDeformedScaleFactor = %12.3e ;\n'], str2num(handles.scaleNum) ); 
  elseif handles.anaVec == 2 || handles.anaVec == 3
		fprintf(fileM, ['controlDofInfo = [ %s ] ;\n'], handles.controldofVec ); 
  end
  fprintf(fileM, ['plotsViewAxis = [ %s ] ;\n'], handles.plotViewVec ); 
  fclose(fileM) ;
  
  cd ..
  guidata(fig, handles) ;

end

%~ Plot button
function plot(src, eventdata, fig) 

  handles = guidata(fig) ;
  figPlot = figure ;
  
  %~ Nodes and Conec
  Nodes = handles.nodesMat(:,1:3) ;
  Conec = [ handles.conecMat(:, 2:5) handles.conecMat(:,7) handles.conecMat(:,6) handles.conecMat(:,1) ] ;
  %~ ---------------------------------------------------------
  plateElems = find(Conec(:,7)==3) ;
  if length(plateElems) > 0
		errordlg('View of solid elements is not available yet', 'View error') ;
		return
  end
  % coordsElemsMat
  nelems = size(Conec,1) ;
  nnodes = size(Nodes,1) ;
  ndofpnode =  6 ;
  coordsElemsMat = zeros(nelems,2*ndofpnode) ;
  for i=1:nelems
		if Conec(i,7) == 1 || 2
			% obtains nodes and dofs of element
			nodeselem = Conec(i,1:2)' ;
			dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
			coordsElemsMat( i, (1:2:11) ) = [ Nodes( nodeselem(1),:)  Nodes( nodeselem(2),:) ] ;
		end	
  end
  %~ ---------------------------------------------------------
  
  %~ Plot options
  lw  = 2   ; ms  = 5.5 ;
  lw2 = 3.2 ; ms2 = 23 ;
  plotfontsize = 22 ;
  hold on, grid on
  for i=1:nelems
    plot3( coordsElemsMat(i,[1 7]), coordsElemsMat(i,[3 9]), coordsElemsMat(i,[5 11]),'b-','linewidth',lw*0.8,'markersize',ms*0.8) ;
  end
  for j = 1:nnodes
    % Nodes circle 
    plot3( Nodes(j,1), Nodes(j,2), Nodes(j,3), 'color', 'r', '.', 'markersize', ms2*0.75 ) ;
    % Nodes label
    text( Nodes(j,1), Nodes(j,2), Nodes(j,3), sprintf( '%i', j ), 'fontsize', plotfontsize*0.75, 'position', [Nodes(j,1)+0.05, Nodes(j,2)+0.05, Nodes(j,3)+0.05]  ) ;
  end
  %~ ---------------------------------------------------------
  
  guidata(fig, handles) ;
end

%~ Run button
function runONSAS(src, eventdata, fig)
  
  handles = guidata(fig) ;
  %~ previouslyDefinedSelectedFileVar = handles.selectedFile ;
  %~ handles.fileName
  
  run(['./input/' handles.fileName]);
  environInputVars = 1 ;
  ONSAS

end

%~ License button
function license(src, eventdata, fig)
  h = msgbox('Read COPYING.txt file for more details.', 'License directory.') ;
end 
