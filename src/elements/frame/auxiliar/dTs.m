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
 
function [Dk]=dTs(t,v);

nt=norm(t);

if nt==0
  Dk=1/2*skew(v);
else
  e=t/nt;
  ev=cross(e,v);

  a1=(cos(nt)-sin(nt)/nt)/nt;
  M1=v*e'-(e'*v)*e*e';

  a2=(1-sin(nt)/nt)/nt;
  M2=e*v'-2*(e'*v)*e*e'+(e'*v)*eye(3);

  a3=sin(nt)/nt-(2*sin(nt/2)/nt)^2;
  M3=ev*e';

  a4=2*(sin(nt/2)/nt)^2;
  M4=skew(v);

  Dk=a1*M1+a2*M2-a3*M3+a4*M4;
end
