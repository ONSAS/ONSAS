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

% ======================================================================
% KT and Finte computation
% ======================================================================

global ne

KTe = zeros(4,4) ;
finte = zeros(4,1) ;

% Elem Gauss points
[xge, we] = gaussPointsAndWeights (ne) ;

a 	= 0 ;
b 	= l ;
p1 	= (b-a)/2 ;
p2 	= (b+a)/2 ;

pgeVec = (p1  * xge' + p2 ) ;	

for j = 1:length(we)
	secFint 	= 0 ;
	secKTe 		= 0 ;
	pge 			= pgeVec(j) ; 
	
	% Bending intern functions second derivative
	B = bendingInterFuns(pge, l, 2) ;
	
	if intBool == 1
		% Tangent stiffness matrix
		secKTe = quadv('secKT', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], ...
														 elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), ...
														 Ut(LocBendXZdofs), hyperElasParams, hyperElasModel) ;
		KTe = p1*( B'*secKTe*B*we(j) ) + KTe ;	
	end
	
	secFint = quadv('secFint', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], ...
															elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), ...
															Ut(LocBendXZdofs), hyperElasParams, hyperElasModel) ;
	finte = p1*secFint*we(j) + finte ;	
													
end % endfor we

%~ if intBool == 1
	%~ KTe = quadv('KTint', 0, l, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), R(LocBendXZdofs,LocBendXZdofs), ...
											 %~ Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, l ) ;
%~ end

% quadv int in the elment and section 
%~ finte = quadv('finteInt', 0, l, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), ...
													%~ R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, ...
													%~ hyperElasModel, l) ;
						
KbendXZ = KTe ;
