
function [ nodesMat, conecMat ] = dxf2ONSAS(nomArch) 

  [c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = f_LectDxf(nomArch);


  ndofpnode = 6 ;
  nelems = size(c_Line,1) ;
  
  Nodes = [] ;
  conecMat = [] ;

  pos = 1 ;
  for i = 1:nelems
    aux1 = [] ;
    aux2 = [] ;  
    nod1 = c_Line{i}(1:3) ;
    nod2 = c_Line{i}(4:6) ;
    nnodes = size(Nodes,1) ;
    
    elemEntity = str2num(c_Line{i,2}(1:2)) ;
    secEntity = str2num(c_Line{i,2}(4:5)) ;
    matEntity = str2num(c_Line{i,2}(7:8)) ;
    loadEntity = str2num(c_Line{i,2}(10:11)) ;
    suppEntity = str2num(c_Line{i,2}(13:14)) ;
    entityVec = [ secEntity matEntity loadEntity suppEntity ] ;
    
    if i == 1
      Nodes = [ Nodes ; nod1 ; nod2 ] ;
      conecMat = [ conecMat ; elemEntity pos pos+1 0 0 entityVec] ;
      pos = pos+2 ;
    else
      for j = 1:nnodes
        if norm(Nodes(j,:)-nod1) < 1e-10 
          aux1 = [ j ] ;
        elseif norm(Nodes(j,:)-nod2) < 1e-10
          aux2 = [ j ] ;
        end
      end
      if length(aux1) == 0 && length(aux2) == 0 
          Nodes = [ Nodes ; nod1 ; nod2 ] ;
          conecMat = [ conecMat ; elemEntity pos pos+1 0 0 entityVec] ;
          pos = pos+2 ;
      elseif length(aux1) == 0 && length(aux2) == 1
          Nodes = [ Nodes ; nod1 ] ;
          conecMat = [ conecMat ; elemEntity pos aux2 0 0 entityVec] ;
          pos = pos+1 ;
      elseif length(aux1) == 1 &&  length(aux2) == 0
          Nodes = [ Nodes ; nod2 ] ;
          conecMat = [ conecMat ; elemEntity aux1 pos 0 0 entityVec] ;
          pos = pos+1 ;
      elseif length(aux1) == 1 && length(aux2) == 1

          conecMat = [ conecMat ; elemEntity aux1 aux2 0 0 entityVec ] ;

      end
         
    end
    
        
  end   

  nnodes = size(Nodes,1) ;
  nnodesEntity = size(c_Poi,1) ;
  auxLength = length(c_Poi{1}) ;
  nodesMat = [ Nodes zeros(nnodes,2) ] ;
  if auxLength > 0
    for i = 1:nnodesEntity
      nod = c_Poi{i}(1:3) ;
      aux = [] ;
      loadEntity = str2num(c_Poi{i,2}(10:11)) ;
      suppEntity = str2num(c_Poi{i,2}(13:14)) ;
      for j = 1:nnodes
        if norm(Nodes(j,:)-nod) < 1e-10 
          nodesMat(j,4) = loadEntity ;
          nodesMat(j,5) = suppEntity ;
        end  
      end
    end
  end

end

function [c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = f_LectDxf(nomArch)
%% Read entities information of dxf file
    %author = SebaTM 
    %Jul 27 2009
    %Based in dxf2coord 1.1 matrix of lukas wischounig, but is not dependent of the Code 
    %Group relative position. That is better way to read dxf format. Don't fail if the 
    %polyline has arcs (budges), but yet don't read them. Don't read arcs as circles. Read 
    %properties (see case 'LINE' by examples of modifications). Group Codes and Associated 
    %Values is read in accurately way (accept associated values with space).
    %
    %Use cell2mat(cell(:,1)) to acquire geometry data in matrix, 
    %by example cell2mat(c_Cir(:,1))
 
%% Read file
		fileName = auxReadDXF(nomArch);
		cd sources
    fId = fopen(fileName);    
    c_ValAsoc = textscan(fId,'%d%s','Delimiter', '\n');
    fclose(fId);
    %~ delete(fileName);
    % Code Group Matrix
    m_GrCode = c_ValAsoc{1};
    % Associated value String Cell
    c_ValAsoc = c_ValAsoc{2};
    %[m_GrCode,c_ValAsoc] = c_ValAsoc{:};
    
%% Entities
    m_PosCero = find(m_GrCode==0);
    %Is searched by (0,SECTION),(2,ENTITIES)
    indInSecEnt = strmatch('ENTITIES',c_ValAsoc(m_PosCero(1:end-1)+1),'exact');
    %(0,ENDSEC)
    m_indFinSecEnt = strmatch('ENDSEC',c_ValAsoc(m_PosCero(indInSecEnt:end)),'exact');
    % Entities Position
    m_PosCero = m_PosCero(indInSecEnt:indInSecEnt-1+m_indFinSecEnt(1));
    % Variable initiation
    %accelerate?
%     c_Line = cell(sum(strcmp('LINE',c_ValAsoc(m_PosCero))),2);
%     c_Poly = cell(sum(strcmp('LWPOLYLINE',c_ValAsoc(m_PosCero))),2);
%     c_Cir = cell(sum(strcmp('CIRCLE',c_ValAsoc(m_PosCero))),2);
%     c_Arc = cell(sum(strcmp('ARC',c_ValAsoc(m_PosCero))),2);
%     c_Poi = cell(sum(strcmp('POINT',c_ValAsoc(m_PosCero))),2);
    c_Line = cell(1,2);
    c_Poly = cell(1,2);
    c_Cir = cell(1,2);
    c_Arc = cell(1,2);
    c_Poi = cell(1,2);
    % 
    iLine = 1;
    iPoly = 1;
    iCir = 1;  
    iArc = 1;
    iPoi = 1;
    % Loop on the Entities
    for iEnt = 1:length(m_PosCero)-2
        m_GrCodeEnt = m_GrCode(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
        c_ValAsocEnt = c_ValAsoc(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
        nomEnt = c_ValAsocEnt{1};  %c_ValAsocEnt{m_PosCero(iEnt+1)}
        %In the entitie's name is assumed uppercase
        switch nomEnt            
            case 'LINE'
                % (Xi,Yi,Zi,Xj,Yj,Zj) start and end points
                c_Line{iLine,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(30,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(11,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(21,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(31,m_GrCodeEnt,c_ValAsocEnt))];
                % Layer
                c_Line(iLine,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
                % Color
                %if no exist is ByLayer (256)
                %c_Line{iLine,3} = str2double(f_ValGrCode(62,m_GrCodeEnt,c_ValAsocEnt));
                % XData
                %c_Line(iLine,4) = f_XData(GroupCode,'XDataName',m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iLine = iLine+1;  
            case 'LWPOLYLINE'
                % (X,Y) vertexs
                %Is not take into account the budge (group code 42, arc in the polyline).
                m_Coord = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                        str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt))];
                if strcmp(f_ValGrCode(70,m_GrCodeEnt,c_ValAsocEnt),'1')&&...
                        any(m_Coord(1,:)~=m_Coord(end,:))
                    %Close polyline
                    c_Poly{iPoly,1} = [m_Coord;m_Coord(1,:)];
                else
                    c_Poly{iPoly,1} = m_Coord;
                end
                % Layer
                c_Poly(iPoly,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iPoly = iPoly+1;   
            case 'CIRCLE'
                % (X Center,Y Center,Radius)
                c_Cir{iCir,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(40,m_GrCodeEnt,c_ValAsocEnt))];
                % Layer
                c_Cir(iCir,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iCir = iCir+1;
            case 'ARC'
                % (X Center,Y Center,Radius,Start angle,End angle)
                c_Arc{iArc,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(40,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(50,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(51,m_GrCodeEnt,c_ValAsocEnt))];
                % Layer
                c_Arc(iArc,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iArc = iArc+1;
            case 'POINT'
                % (X,Y,Z) Position
                c_Poi{iPoi,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(30,m_GrCodeEnt,c_ValAsocEnt))];
                % Layer
                c_Poi(iPoi,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iPoi = iPoi+1;
            %case Add Entities
        end        
    end 
		cd ..
%%   
end
%%
function c_Val = f_ValGrCode(grCode,m_GrCode,c_ValAsoc)
    c_Val = c_ValAsoc(m_GrCode==grCode);
end
%%
function c_XData = f_XData(grCode,XDatNom,m_GrCode,c_ValAsoc)
    m_PosXData = find(m_GrCode==1001);
    if ~isempty(m_PosXData)
        indInXData = m_PosXData(strmatch(upper(XDatNom),c_ValAsoc(m_PosXData),'exact'));
        m_indFinXData = find(m_GrCode(indInXData+2:end)==1002)+indInXData+1;
        m_IndXData = indInXData+2:m_indFinXData(1)-1;
        c_XData = f_ValGrCode(grCode,m_GrCode(m_IndXData),c_ValAsoc(m_IndXData));
    else
        c_XData = {[]};
    end
end
