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

function [F1, F2] = matrixF(n, flag_second_mod) % ok
  % Eq. (30) of 10.1016/j.cma.2006.10.006
  n1 = n(1:3) ;
  n2 = n(4:6) ;
  %
  n3 = n(7:9) ;
  n4 = n(10:12) ;
  %
  n5 = n(13:15) ;
  n6 = n(16:18) ; 
  if flag_second_mod == 1
    F1 = [ skew(n1)(1:2,:)' zeros(3,3) skew(n3)(1:2,:)' zeros(3,3) skew(n5)(1:2,:)' zeros(3,3) ]' ;
  else
    F1 = [ skew(n1)' zeros(3,3) skew(n3)' zeros(3,3) skew(n5)' zeros(3,3) ]' ;
    % F1 = [ skew(n1) ; zeros(3,3) ; skew(n3) ; zeros(3,3) ; skew(n5) ; zeros(3,3) ] ;
  end
  F2 = [ skew(n1)' skew(n2)' skew(n3)' skew(n4)' skew(n5)' skew(n6)' ]' ;
    % F2 = [ skew(n1) ; skew(n2) ; skew(n3) ; skew(n4) ; skew(n5) ; skew(n6) ] ;
end