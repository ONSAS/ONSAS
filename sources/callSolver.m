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


%This functions calls the corresponding solver according to the anlysis settings and the numerical method provided by the user.

function  [ modelCurrState, BCsNextState, auxIO ]  = callSolver( modelCurrState, BCsNextState, auxIO ) ;

modelExtract

if dynamicAnalysisBoolean == 0

  currDeltau = zeros(length(neumdofs),1) ;

  [ nextLoadFactor, itersPerTime, stopCritPar, factorCrit, nKeigpos, nKeigneg, Ut, FintGt, Stresst, Strainst , dsigdepst ] ...
  = analysisNRAndNRAL ( ...
    Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat,  ...
    numericalMethodParams, constantFext, variableFext, KS, userLoadsFilename, ...
    Ut, Stresst, Strainst, dsigdepst, FintGt, currLoadFactor, nextLoadFactor, ...
    convDeltau ) ;

  modelCompress

else
  
  deltaT         = numericalMethodParams(2)        ;
  finalTime      = numericalMethodParams(3)        ;
  stopTolDeltau  = numericalMethodParams(4)        ;
  stopTolForces  = numericalMethodParams(5)        ;
  stopTolIts     = numericalMethodParams(6)        ;
  deltaNW        = numericalMethodParams(7)        ; % add this to documentation
  AlphaNW        = numericalMethodParams(8)        ; % add this to documentation
  nextLoadFactor = loadFactorsFunc(currTime+deltaT);
  

  if currTime == 0 ;
  
    % assemble M and C
    massMatAssembly                                ;
    dampingMat = eye(size(massMat)) * nodalDamping   ;
    % podriai ir en el initial y agregarse al model compress
    
    %inital and constant parameters  
    a0NW = 1/(AlphaNW*(deltaT)^2)        ;
    a1NW = deltaNW/(AlphaNW*deltaT)      ;
    a2NW = 1/(AlphaNW*deltaT)            ;
    a3NW = 1/(2*AlphaNW)-1               ;
    a4NW = deltaNW/AlphaNW -1            ;
    a5NW = (deltaT/2)*(deltaNW/AlphaNW-2);
    a6NW = deltaT*(1-deltaNW)            ;
    a7NW = deltaNW*deltaT                ;
        
    
    if strcmp( userLoadsFilename , '')
      FextUser = zeros(size(constantFext));
    else
      FextUser = feval( userLoadsFilename, currTime)  ;
    end

    %Fext mas, damping and stifness matrix reduced
    massMatred        = massMat(neumdofs,neumdofs)                    ;
    dampingMatred     = dampingMat(neumdofs,neumdofs)                 ; 
    currLoadFactor    = loadFactorsFunc(currTime)                     ;
    FextGt            = variableFext * currLoadFactor + constantFext +  FextUser  ;
    FextGtred         = FextGt(neumdofs)                              ;
  
    FintGt            = assemblyFintVecTangMat( Conec, secGeomProps, ...
                                      coordsElemsMat, hyperElasParamsMat, KS, Ut,1 ); %el ultimo 1 saca la fuerza
                                             
    FintGtred         = FintGt(neumdofs)                              ;
  

    Udotdott(neumdofs)= massMatred \ (FextGtred-dampingMatred*Udott(neumdofs)-FintGtred) ;

  end

  % assemble M and C
  massMatAssembly                                ;
  dampingMat = eye(size(massMat))* nodalDamping   ;


  %inital and constant parameters  
  a0NW = 1/(AlphaNW*(deltaT)^2)        ;
  a1NW = deltaNW/(AlphaNW*deltaT)      ;
  a2NW = 1/(AlphaNW*deltaT)            ;
  a3NW = 1/(2*AlphaNW)-1               ;
  a4NW = deltaNW/AlphaNW -1            ;
  a5NW = (deltaT/2)*(deltaNW/AlphaNW-2);
  a6NW = deltaT*(1-deltaNW)            ;
  a7NW = deltaNW*deltaT                ;
  
  % la matriz y lo ai de newmark podrian ir en initial y definirse en model compress
  [ Utp1, Udottp1, Udotdottp1, FintGtp1, dispIter, Strainst, Stresst ] ...
    = analysisNM ( ...
    % --- inputs ---
    % constant data
    Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat, ...
    constantFext, variableFext, KS, ...
    %
    massMat, dampingMat, a0NW, ...
    a1NW, a2NW, a3NW, a4NW,...
    a5NW, a6NW,   a7NW, ...
    % model variable data
    dispsElemsMat, Ut, Udott, Udotdott, nextLoadFactor, stopTolDeltau, stopTolForces, stopTolIts, userLoadsFilename, currTime + deltaT ) ;

  % Releases displacements velocity and aceleration
  Ut       = Utp1              ;
  Udott    = Udottp1           ;
  Udotdott = Udotdottp1        ;
  currTime = currTime + deltaT ;    
  
  modelCompress
  

end
