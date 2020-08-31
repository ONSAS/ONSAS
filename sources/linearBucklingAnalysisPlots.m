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

% --------------------------------------------------------------------------------------------------

% ----------- plots ----------------------------

lw = 1.5;
ms = 3;

linfig = figure ;
title('Linear analysis results')

subplot(2,2,1)
hold on, grid on

displacmat = reshape(Uklin,2,nnodes)' / norm( Uklin, 'inf' ) * (0.1 * strucsize) ;
Nodesdefaux = Nodes + displacmat  ;

for i=1:nelems
  plot(xelems(i,:),yelems(i,:),'b-o','linewidth',lw,'markersize',ms)
end

for i=1:nelems
  xelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 1);
  yelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 2);
end

for i=1:nelems
  plot(xelemsaux(i,:),yelemsaux(i,:),'r-s','linewidth',lw,'markersize',ms)
end

axis equal
title('Deformed structure (scaled)')

subplot(2,2,2); hold on, grid on

tolstress = 1e-7*max(abs(stresslink));

for i=1:nelems
  if     stresslink(i) > tolstress, colornormalforce='b';
  elseif stresslink(i) <-tolstress, colornormalforce='r';
  else                             colornormalforce='w'; end

  plot(xelems(i,:),yelems(i,:),[colornormalforce '-o'],'linewidth',lw,'markersize',ms)
end

axis equal
title('Normal forces: >0 blue, <0 red, ==0 white.')

% ----------- plots ----------------------------
subplot(2,2,3), hold on, grid on

modeaux = 1;

displacmat = reshape(UsModesLinBuck(:,modeaux),2,nnodes)' * (0.1*strucsize) ;
Nodesdefaux = Nodes + displacmat  ;


for i=1:nelems
  plot(xelems(i,:),yelems(i,:),'b-o','linewidth',lw,'markersize',ms)
end

for i=1:nelems
  xelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 1);
  yelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 2);
end

for i=1:nelems
  plot(xelemsaux(i,:),yelemsaux(i,:),'r-s','linewidth',lw,'markersize',ms)
end

axis equal, title('LBA mode 1')

% ----------- plots 2 ----------------
if length(lambdas)>1

  subplot(2,2,4), hold on, grid on

  modeaux = 2;
  displacmat = reshape(UsModesLinBuck(:,modeaux),2,nnodes)' * ( 0.1*strucsize);
  Nodesdefaux = Nodes + displacmat  ;
  
  for i=1:nelems
    plot(xelems(i,:),yelems(i,:),'b-o','linewidth',lw,'markersize',ms)
  end
  
  for i=1:nelems
    xelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 1);
    yelemsaux(i,:) = Nodesdefaux( Conec(i,1:2) , 2);
  end
  
  for i=1:nelems
    plot(xelemsaux(i,:),yelemsaux(i,:),'r-s','linewidth',lw,'markersize',ms)
  end
  axis equal, title('LBA mode 2')
end
