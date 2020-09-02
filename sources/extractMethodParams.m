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


% ======================================================================
function [ solutionMethod, stopTolDeltau,   stopTolForces, ...
           stopTolIts,     targetLoadFactr, nLoadSteps,    ...
           incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
           = extractMethodParams( numericalMethodParams ) 

  solutionMethod   = numericalMethodParams(1) ;
  
  if solutionMethod == 0

    % ----- resolution method params -----
    stopTolDeltau    = 0     ;
    stopTolForces    = 1e-10 ;
    stopTolIts       = 2     ;
    targetLoadFactr  = 1     ;
    nLoadSteps       = 1     ; 

    deltaT = targetLoadFactr/nLoadSteps ; finalTime = targetLoadFactr ;

    incremArcLen = [] ; deltaNW = []; AlphaNW = [] ; alphaHHT = [] ;

  elseif ( solutionMethod == 1) || ( solutionMethod == 2)

    % ----- resolution method params -----
    stopTolDeltau    = numericalMethodParams(2) ;
    stopTolForces    = numericalMethodParams(3) ;
    stopTolIts       = numericalMethodParams(4) ;
    targetLoadFactr  = numericalMethodParams(5) ;
    nLoadSteps       = numericalMethodParams(6) ;
  

    if solutionMethod ==2
      incremArcLen     = numericalMethodParams(7) ;
    else
      incremArcLen = [] ;
    end
    
    deltaT = targetLoadFactr/nLoadSteps ; finalTime = targetLoadFactr ;
    
    deltaNW = []; AlphaNW = [] ; alphaHHT = [] ;
  
  elseif solutionMethod == 3
    deltaT         = numericalMethodParams(2)        ;
    finalTime      = numericalMethodParams(3)        ;
    stopTolDeltau  = numericalMethodParams(4)        ;
    stopTolForces  = numericalMethodParams(5)        ;
    stopTolIts     = numericalMethodParams(6)        ;
    deltaNW        = numericalMethodParams(7)        ;
    AlphaNW        = numericalMethodParams(8)        ;

    alphaHHT = [] ;
    targetLoadFactr = [] ; nLoadSteps = []; incremArcLen = [] ;    

  elseif solutionMethod == 4
  
    deltaT         = numericalMethodParams(2)        ;
    finalTime      = numericalMethodParams(3)        ;
    stopTolDeltau  = numericalMethodParams(4)        ;
    stopTolForces  = numericalMethodParams(5)        ;
    stopTolIts     = numericalMethodParams(6)        ;
    alphaHHT       = numericalMethodParams(7)        ;
    
    deltaNW = []; AlphaNW = [] ;
    targetLoadFactr = [] ; nLoadSteps = []; incremArcLen = [] ;    
    
  end
  
% ======================================================================
