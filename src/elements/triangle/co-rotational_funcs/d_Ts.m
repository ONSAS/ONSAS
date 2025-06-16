% Copyright 2024, ONSAS Authors (see documentation)
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
%
function [d_Ts] = d_Ts(t, v)



  % nt = norm(t);

  % if nt == 0
  %   d_Ts = 1 / 2 * skew(v);
  % else
  %   e = t / nt;
  %   ev = cross(e, v);

  %   a1 = (cos(nt) - sin(nt) / nt) / nt;
  %   M1 = v * e' - (e' * v) * e * e';

  %   a2 = (1 - sin(nt) / nt) / nt;
  %   M2 = e * v' - 2 * (e' * v) * e * e' + (e' * v) * eye(3);

  %   a3 = sin(nt) / nt - (2 * sin(nt / 2) / nt)^2;
  %   M3 = ev * e';

  %   a4 = 2 * (sin(nt / 2) / nt)^2;
  %   M4 = skew(v);

  %   d_Ts = a1 * M1 + a2 * M2 - a3 * M3 + a4 * M4;
  % end

% fprintf('=== \n')
  % d_Ts
  % fprintf('--------------------------------- \n')
  psi = norm(t);
  I = eye(3, 3);

  if psi == 0
    d_Ts = 1 / 2 * skew(v);
  else
    u = t / psi ;
    a = sin(psi) / psi ;
    b = ( sin(psi/2) / (psi/2) ) ;
    sum1 = -(a - b^2) * cross(u,v) * u' ;
    sum2 = 1/2*b^2 * skew(v) ;
    sum3 = ( cos(psi) -a ) * 1/psi * (v*u'- (u'*v)*u*u' ) ;
    sum4 = (1-a)* 1/psi * ( u*v' -2*(u'*v)*u*u' + (u'*v)*I ) ; 
    d_Ts = sum1 + sum2 + sum3 + sum4 ;
  end
  % d_Ts
  % fprintf('=== \n')
end
