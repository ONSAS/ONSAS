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


%script for ploting load factor versus control displacement 

lw = 1.5; ms = 3;

flagMarkersLoadDisp = 1;

loadFactors 

timeVals 
stop
figure , hold on, grid on
plot( timeVals, loadFactors        , 'b-x', 'linewidth', lw,'markersize',ms)
labx = xlabel('step/time'); laby = ylabel('Load factors');
set(gca, 'linewidth', 1, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize);

cd( outputdir )
if printflag == 1
  print( [ problemName '_timeLoad'] ,'-depslatex') ;
elseif printflag == 2
  print( [ problemName '_timeLoad'] ,'-dpng') ;
end
cd(currdir)

% ===================================================

figure , hold on, grid on
plot( timeVals , controlDisps      , 'b-x', 'linewidth', lw,'markersize',ms)
labx = xlabel('step/time'); laby = ylabel('Displacement');
set(gca, 'linewidth', 1, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize);

cd( outputdir )
if printflag == 1
  print( [ problemName '_timeDisp'] ,'-depslatex') ;
elseif printflag == 2
  print( [ problemName '_timeDisp'] ,'-dpng') ;
end
cd(currdir)


% ===================================================
figure, hold on, grid on
plot( controlDisps, loadFactors , 'b-x', 'linewidth', lw,'markersize',ms)

if flagMarkersLoadDisp == 1
  for indplot = 1 : length( timesPlotsVec ) ;
    text( controlDisps( timesPlotsVec( indplot)), loadFactors( timesPlotsVec( indplot) ), sprintf('%6.2f', timesPlotsVec(indplot)/nTimesTotal) ) 
    plot( controlDisps( timesPlotsVec( indplot)), loadFactors( timesPlotsVec( indplot) ), 'ro','markersize',ms*3 ) 
  end
end

labx=xlabel('Displacements'); laby=ylabel('Load factors');
set(gca, 'linewidth', 1, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize);

cd( outputdir )
if printflag == 1
  print( [ problemName '_loadDisp'] ,'-depslatex') ;
elseif printflag == 2
  print( [ problemName '_loadDisp'] ,'-dpng') ;
end
cd(currdir)
% ===================================================
