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
 
function [Dh]=dinvTs(t,v);

nt = norm(t) ;

if nt==0
  Dh=-1/2*skew(v);
else
  a=nt/2;
  eta=(sin(a)-a*cos(a))/(nt^2*sin(a));
  miu=(nt*(nt+sin(nt))-8*sin(a)^2)/(4*nt^4*sin(a)^2);
  I3=eye(3);
  M=skew(t);
  M1=skew(v);
  M2=t*v'-2*v*t'+(t'*v)*I3;
  M3=M*M*v*t';
  Dh=eta*M2+miu*M3-1/2*M1;
end
