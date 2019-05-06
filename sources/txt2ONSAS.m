function [nodesMat conecMat] = txt2ONSAS(fname)
  fid = fopen(fname, 'r') ;

  % Nodes title
  title = fgets(fid) ;
  nodesVec = fscanf(fid, '%g'); 
  % Conec title
  title = fgets(fid) ;
  conecVec = fscanf(fid, '%g');
  fclose(fid) ;
  
  aux = length(nodesVec);
  nnodes = aux/5 ;
  nodesMat = [] ;
  for i = 1:nnodes
    nodesMat = [nodesMat ; nodesVec((i-1)*5+1:i*5)'] ;
  end
  
  aux = length(conecVec) ;
  nelems = aux/7 ;
  conecMat=[];
  for i = 1:nelems
    conecMat = [conecMat ; conecVec((i-1)*9+1:i*9)'] ;
  end

end
