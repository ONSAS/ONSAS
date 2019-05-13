%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


%script for generation of plots of deformed structure.

currdir = pwd ;
lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

xs = Nodes(:,1) ;
ys = Nodes(:,2) ;
zs = Nodes(:,3) ;

marginAxis = 0.1*strucsize ;

minxdef = min( xs + min( linearDeformedScaleFactor*matUts(1:6:end,:)' )' ) - marginAxis ;
minydef = min( ys + min( linearDeformedScaleFactor*matUts(3:6:end,:)' )' ) - marginAxis ;
minzdef = min( zs + min( linearDeformedScaleFactor*matUts(5:6:end,:)' )' ) - marginAxis ;

maxxdef = max( xs + max( linearDeformedScaleFactor*matUts(1:6:end,:)' )' ) + marginAxis ;
maxydef = max( ys + max( linearDeformedScaleFactor*matUts(3:6:end,:)' )' ) + marginAxis ;
maxzdef = max( zs + max( linearDeformedScaleFactor*matUts(5:6:end,:)' )' ) + marginAxis ;

nelems = size(Conec,1) ;
nnodes = size(Nodes,1) ;
ndofpnode = 6 ;

for indplot = 1 : length( timesPlotsVec ) ;

  figdef = figure ;
  if size(matUts,2) > 1
		if nonLinearAnalysisBoolean && dynamicAnalysisBoolean ~= 0
			subplot(3,2,1:4)
		end	
  end
  
  hold on, grid on

  for i=1:nelems
    plot3( coordsElemsMat(i,[1 7]), coordsElemsMat(i,[3 9]), coordsElemsMat(i,[5 11]),'b--o','linewidth',lw*0.8,'markersize',ms*0.8);
  end


  Utplot = matUts ( :, timesPlotsVec( indplot) ) *linearDeformedScaleFactor ;

  dispsElemsMat = zeros(nelems,2*ndofpnode) ;
  for i=1:nelems
    % obtains nodes and dofs of element
    nodeselem = Conec(i,1:2)' ;
    dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
    dispsElemsMat( i, : ) = Utplot(dofselem)' ;
  end
	
  %~ NodesDef = Nodes + reshape()
  
	aux = zeros( nelems,2*ndofpnode ) ;
	for i = 1:nelems
		for j = 1:12
			aux(i,j) = sum ( dispsElemsMat(i,j,:) )'; 
		end
	end
	
	
  for i=1:nelems
		if Conec(i,end) == 1 || Conec(i,end) == 2 

  		[~, locglomat] = beamParameters(Nodes(Conec(i,1:2),:)) ;
	
  		if (nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0)
				[ xselemdef, yselemdef, zselemdef ] = outputFrameElementPlotLin ( coordsElemsMat(i,:)', aux(i,:)', Conec(i,end), locglomat ) ;

			else
				[ xselemdef, yselemdef, zselemdef ] = outputFrameElementPlot ( coordsElemsMat(i,:)', dispsElemsMat(i,:)', Conec(i,end) ) ;  

			end
			plot3( xselemdef, yselemdef, zselemdef, 'b-', 'linewidth', lw, 'markersize', ms );

		else
			% Deformed shape of plate
    end
  end

  FG = variableFext ;
  UG = Utplot ;

  quiver3( Nodes(:,1)+UG(1:6:end), Nodes(:,2)+UG(3:6:end), Nodes(:,3)+UG(5:6:end) , ...
      FG(1:6:end)*visualloadfactor , ...
      FG(3:6:end)*visualloadfactor , ...
      FG(5:6:end)*visualloadfactor , ...
      0,'c',"filled",'linewidth',lw2)

  FG = constantFext ;
  UG = Utplot ;

  quiver3( Nodes(:,1)+UG(1:6:end), Nodes(:,2)+UG(3:6:end), Nodes(:,3)+UG(5:6:end) , ...
      FG(1:6:end)*visualloadfactor , ...
      FG(3:6:end)*visualloadfactor , ...
      FG(5:6:end)*visualloadfactor , ...
      0,'m',"filled",'linewidth',lw2)
	
	if size(matUts,2) == 1
		tit = title(['Deformed shape' ] );
	else
		tit = title(['Deformed increment: ' sprintf('%04i', timesPlotsVec( indplot)) '/' sprintf('%04i', nTimesTotal) ] );
	end
  labx=xlabel('x'); laby=ylabel('y'); labz=zlabel('z') ;
  set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize*0.5 ) ;
  set(tit, "FontSize", plotfontsize) ;
  set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize) ; set(labz, "FontSize", plotfontsize) ;

  axis equal

  % sets global or coordinate-wise axis limits

  axis( [ minxdef maxxdef   minydef maxydef   minzdef maxzdef ] );

  if length( plotsViewAxis ) == 0
    view(3);
  elseif sum( plotsViewAxis ~= [0 0 1] ) == 0
    view(2);
  else
    view(plotsViewAxis);
  end

  if size(matUts,2) > 1
		if nonLinearAnalysisBoolean && dynamicAnalysisBoolean ~= 0
			subplot(3,2,5)
			plot( controlDisps(1:timesPlotsVec(indplot)) , loadFactors(1:timesPlotsVec(indplot))  , 'b-x', 'linewidth', lw*0.4,'markersize',ms*0.4)
			axis( [ min( controlDisps) max(controlDisps) min( loadFactors) max( loadFactors) ] ); grid on
			labx = xlabel('control displacement'); laby = ylabel('load factor');
			
			subplot(3,2,6)
			plot( loadFactors(1:timesPlotsVec(indplot))        , 'b-x', 'linewidth', lw*0.4,'markersize',ms*0.4)
			axis( [ 0 length(loadFactors) min(loadFactors) max(loadFactors) ] ); grid on
			labx=xlabel('step'); laby=ylabel('load factor');
		end
  end
  % ---------------------------------------------------------------------

  cd(outputdir )
  if printflag == 1
		if size(matUts,2) == 1
			print( [ problemName '_deform' ] ,'-depslatex') ;
		else
			print( [ problemName '_deform_' sprintf('%04i', indplot)  ] ,'-depslatex') ;
		end	
  elseif printflag == 2
		if size(matUts,2) == 1
			print( [ problemName '_deform' ] ,'-dpng') ;
		else
			print( [ problemName '_deform_' sprintf('%04i', indplot) ] ,'-dpng') ;
		end
  end
  cd(currdir)

  if printflag > 0  
    close(figdef);
  end
  % ---------------------------------------------------------------------


  % ----------------------------------------
  % ---------- Axial force plots  ----------
  if length(plotParamsVector) > 0 && plotParamsVector(1)>0
    
    figAxial = figure ;
    hold on, grid on

    cmap = flipud( colormap('jet') ) ; % other good options: hot
    colormap(cmap);


    axis equal
    axis( [ minxdef maxxdef   minydef maxydef   minzdef maxzdef ] );

    if length( plotsViewAxis ) == 0
      view(3);
    elseif sum( plotsViewAxis ~= [0 0 1] ) == 0
      view(2);
    else
      view(plotsViewAxis);
    end

    normalForce    = matNts(:,timesPlotsVec( indplot)) ;
    minNormalForce = min( normalForce);
    maxNormalForce = max( normalForce);
    
    for i = 1:nelems
 
      if Conec(i,end) == 1 || Conec(i,end) == 2
        nodeselem = Conec(i,1:2) ;
        [lengths, ~] = beamParameters(Nodes(nodeselem,:)) ;
        offsetText = min(lengths) / 15 ;
        
        posText = ( Nodes(nodeselem(2),:)+Nodes(nodeselem(1),:) ) / 2 ;
        
        if abs(maxNormalForce - minNormalForce) < 1e-10
          cmapi = cmap( 1 ,: );
        else
          cmapi = cmap( max( [ 1 ceil( (normalForce(i)-minNormalForce) / abs( maxNormalForce-minNormalForce) * length(cmap) ) ] ) ,: );
        end

        % --- plot of each element
        if Conec(i,end) == 1 || Conec(i,end) == 2 
    
          [~, locglomat] = beamParameters(Nodes(Conec(i,1:2),:)) ;
      
          if (nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0)
            [ xselemdef, yselemdef, zselemdef ] = outputFrameElementPlotLin ( coordsElemsMat(i,:)', aux(i,:)', Conec(i,end), locglomat ) ;
    
          else
            [ xselemdef, yselemdef, zselemdef ] = outputFrameElementPlot ( coordsElemsMat(i,:)', dispsElemsMat(i,:)', Conec(i,end) ) ;  
    
          end
          plot3( xselemdef, yselemdef, zselemdef, 'color',cmapi,'linewidth',lw*0.7 );
    
        else
          % Deformed shape of plate
        end
        % ----------

        text( posText(1)+offsetText, posText(2)+offsetText, posText(3)+offsetText, sprintf( '%5.1e', normalForce(i)), 'color', cmapi, 'fontsize', 9 )
        
      else
        fprintf('Missing: normal forces plot in octave for beam elements.\n')
      
      end
    
    end

    colorbar('title','Normal Force')
    if minNormalForce ~= maxNormalForce
      caxis([minNormalForce maxNormalForce])
    end
    
    if size(matUts,2) == 1
      tit = title(['Deformed shape' ] );
    else
      tit = title(['Step/increment: ' sprintf('%04i', timesPlotsVec( indplot)) '/' sprintf('%04i', nTimesTotal) ] );
    end
    labx=xlabel('x'); laby=ylabel('y'); labz=zlabel('z') ;
    set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize*0.5 ) ;
    set(tit, "FontSize", plotfontsize) ;
    set(labx, "FontSize", plotfontsize); set(laby, "FontSize", plotfontsize) ; set(labz, "FontSize", plotfontsize) ;
    
    
  end 
end %endfor indplot
