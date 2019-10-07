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

% ----------------------------------------------------------------------------
% GUI geometry
% ----------------------------------------------------------------------------

% User screensize
anchoPantalla = get(0,'screensize')(3) ;
altoPantalla  = get(0,'screensize')(4) ;

% Main window size
anchoMenu = .8 * anchoPantalla ;
altoMenu  = .8 * altoPantalla ;
posXMenu = .1 * anchoPantalla ;
posYMenu = .1 * altoPantalla ;
vecFig = [posXMenu posYMenu anchoMenu altoMenu] ; % pixels

% Uipanel position variables
sysActionsFromBottom = 0 ;
sysActionsFromLeft = 0 ;
sysActionsWidth = 1 ;
sysActionsHeight = 0.1 ;
vecSysActions = [sysActionsFromLeft, sysActionsFromBottom, sysActionsWidth, sysActionsHeight] ; % normalized

% Input panel
inpPanelFromBottom = sysActionsHeight ;
inpPanelFromLeft = 0 ;
inpPanelWidth = 1 ;
inpPanelHeight = 1 - inpPanelFromBottom ;
vecInpPanel = [inpPanelFromLeft, inpPanelFromBottom, inpPanelWidth, inpPanelHeight] ; % normalized

% Analysis panel
anaPanelFromBottom = 0 ;
anaPanelFromLeft = 0 ;
anaPanelWidth = 0.3 ;
anaPanelHeight = 1 ;
vecAnaPanel = [anaPanelFromLeft, anaPanelFromBottom, anaPanelWidth, anaPanelHeight] ; % normalized

% Geommetry panel
geomPanelFromBottom = 0.865 ;
geomPanelFromLeft = anaPanelWidth ;
geomPanelWidth = 1-anaPanelWidth ;
geomPanelHeight = 1-geomPanelFromBottom ;
vecGeomPanel = [geomPanelFromLeft, geomPanelFromBottom, geomPanelWidth, geomPanelHeight] ; % normalized

% Structural properties panel
strucPanelFromBottom = geomPanelFromBottom - geomPanelHeight ;
strucPanelFromLeft = geomPanelFromLeft ; 
strucPanelWidth = geomPanelWidth ;
strucPanelHeight = geomPanelHeight ; 
vecStrucPanel = [strucPanelFromLeft, strucPanelFromBottom, strucPanelWidth, strucPanelHeight] ; % normalized


% Boundary Conditions panel
BCPanelFromBottom = geomPanelFromBottom - 2 * geomPanelHeight ;
BCPanelFromLeft = geomPanelFromLeft ; 
BCPanelWidth = geomPanelWidth ;
BCPanelHeight = geomPanelHeight ; 
vecBCPanel = [BCPanelFromLeft, BCPanelFromBottom, BCPanelWidth, BCPanelHeight] ;

% Output panel
outputPanelFromBottom = geomPanelFromBottom - 3.5 * geomPanelHeight ;
outputPanelFromLeft = geomPanelFromLeft ; 
outputPanelWidth = geomPanelWidth ;
outputPanelHeight = 1.5*geomPanelHeight ; 
vecOutputPanel = [outputPanelFromLeft, outputPanelFromBottom, outputPanelWidth, outputPanelHeight] ;

% Pixels alignment
PPFromLeft = 0.015 * anchoPantalla ;
PPFromBottom = 0.01 * altoPantalla ;

% Uicontrol position variables

ratioSizeButtons = 0.25 ;

% Input Panel buttons
buttonWidth = .08 * anchoPantalla ;
buttonHeight = ratioSizeButtons * buttonWidth ;
incremHeight = 25 ;
%~ incremHor = 115 ;
incremHor = 0.085 * anchoPantalla ;

% Action panel button
buttonSysHeight = buttonHeight ;
buttonSysWidth  = buttonWidth ;
incremSysHeight = 25 ;
%~ incremSysHor = 115 ; 

% Editboxes size
editW = buttonWidth ;
editH = buttonHeight ;
textW = buttonWidth*1.5 ;
textH = buttonHeight ;
%~ editSep = 0.12 * altoPantalla ;
editSep = 0.2 * altoPantalla ;

% Sys Action
incremSysHor = 0.085 * anchoPantalla ;
PPFromBottomSys = 0.0225 * altoPantalla ;

% Titles ref from bottom pannel 
ref = 0.0257 * altoPantalla + buttonHeight;

fontSize = 10 * (anchoPantalla<1500) + 10.5 * (anchoPantalla>=1500)*(anchoPantalla<2500) + 11 * (anchoPantalla>=2500) ;
fontSizeTitle = 11 * (anchoPantalla<1500) + 11.5 * (anchoPantalla>=1500)*(anchoPantalla<2500) + 12.5 * (anchoPantalla>=2500) ;
fontSizeInpTitle = 16 * (anchoPantalla<1500) + 17 * (anchoPantalla>=1500)*(anchoPantalla<2500) + 18 * (anchoPantalla>=2500) ;

% ------------------------------------------------------------------------------------------------------------------------
