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
