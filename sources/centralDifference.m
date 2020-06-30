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

function ...
%  outputs ---
[ Utp1, Udott, Udotdott ] ...
  = centralDifference ( ...
% inputs ---
  % constant data
  Conec, secGeomProps, coordsElemsMat, neumdofs, hyperElasParams, ...
  constantFext, variableFext, KS, ...
  %
  massMat, dampingMat, Utm1, ...
  a0CD, a1CD, a2CD, a3CD, ...
  % model variable data
  dispsElemsMat, Ut, currLoadFactor ) ;
% ----------------


  ndofpnode = 6;
  nelems = size(Conec,1) ;

  FintGt = internalForce( Conec, secGeomProps, coordsElemsMat, dispsElemsMat, hyperElasParams, KS, Ut );
  
  FextGt  = variableFext * currLoadFactor + constantFext ;

  Fhat  = FextGt - FintGt + a2CD * massMat * Ut - ( a0CD * massMat - a1CD * dampingMat )* Utm1 ;
  Mhat  = a0CD * massMat + a1CD * dampingMat ;
  
  Mhatred = Mhat(neumdofs,neumdofs) ;
  Fhatred = Fhat(neumdofs) ;

  Utp1red = Mhatred \ Fhatred ;

  Utp1 = zeros(size(Ut));
  Utp1(neumdofs) = Utp1red ;
  
  Udotdott = a0CD * ( Utp1 -2*Ut + Utm1 ) ;
  Udott    = a1CD * ( Utp1       - Utm1 ) ;
  
