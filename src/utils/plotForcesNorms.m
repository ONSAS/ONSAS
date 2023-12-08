% Copyright 2023, ONSAS Authors (see documentation)
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

% This functions reads a text file with the norms of the forces at
%  each iteation of the code execution and generates a series of plots

function plotForcesNorms(outputDir, problemName )

filename = [ outputDir  problemName '_forcesNorms.txt' ] ;

matriz = load(filename);

times = matriz(:,1);
[~,inds_last_iters] = unique( times,'last');

iters = matriz(:,2);

times = matriz(:,end);


indexes = inds_last_iters ;

figure

subplot(5,1,1)
plot(times(indexes), matriz(indexes,4),'b-x')
grid on
title('external forces')

subplot(5,1,2)
plot(times(indexes),matriz(indexes,5),'b-x')
grid on
title('internal forces')

subplot(5,1,3)
plot( times(indexes), matriz(indexes,7) ,'b-x')
grid on
title('inertial forces')

subplot(5,1,4)
plot( times(indexes), matriz(indexes,3) ,'b-x')
grid on
title('residual forces')

subplot(5,1,5)
plot( times(indexes), iters(indexes) ,'b-x')
grid on
title('iterations')

print([ outputDir  problemName '_forcesNorms_plot_converged.png'], '-dpng');


figure


finalIters = min([ 300, length(times)-1]) ;

subplot(4,1,1)
plot(matriz((end-finalIters):end,4),'b-x')
grid on
title('external forces')

subplot(4,1,2)
plot(matriz((end-finalIters):end,5),'b-x')
grid on
title('internal forces')

subplot(4,1,3)
plot(matriz((end-finalIters):end,7),'b-x')
grid on
title('inertial forces')

subplot(4,1,4)
plot(matriz((end-finalIters):end,3),'b-x')
grid on
title('residual forces')

print([ outputDir  problemName '_forcesNorms_plot_all.png'], '-dpng');

