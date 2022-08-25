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

function [ b ] = secWidth(p, type, coef, h)
	
  % sec types
  % ------------------------------------------------------
  % Type 1
  % Section conformed by discrete changes in the width (i.e section composed by arrangement of rectangles)
  % coef struc:
  % coef = [ bn hn ; bn-1 hn-1 bn-2 hn-2 ... b2 h2 b1 h1 ] ;
  
  
  % Type 2
	% Section width is described by a polinomial function 
	% coef of length n struc:
	% coef = [ an*y^n-1 + an-2^yn-2 ... a2*y^2 + a1*y^1 + a0 ] ;
	
  
     
  if type == 1
    coef = reshape(coef,3,length(coef)/3)' ;
    secElements = size(coef,1) ;
    stopParam = 0 ;
		secElem = 1 ;
    while stopParam == 0 && secElem <= secElements

			h1 = coef(secElem,2) ;
			h2 = coef(secElem,3) ;
			hinf = min(h1,h2) ;
			hsup = max(h1,h2) ;
			if (p<=hsup) && (p>=hinf)
				b = coef(secElem,1) ;
				stopParam = 1 ;
			end
			secElem = secElem+1 ;
    end  
    
  elseif type == 2 
		%~ coef = coef(2:end) ;
    bool = length(coef) ;
  
    switch bool
    
      case 1
        b = coef * ones(length(p),1) ;
      case 2
        b = coef(1)*p + coef(2) ;
    end	
	end 

end
