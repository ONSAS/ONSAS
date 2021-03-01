% auxiliar function to convert a cell with vectors to a zeros-filled-matrix

function Mat = myCell2Mat( Cell )

if iscell( Cell )
  nCellRows = size(  Cell,      1 ) ; 
  Mat       = zeros( nCellRows, 1 ) ;
  for i = 1:nCellRows
    aux                   = Cell{i,1} ;
    Mat ( i,1:length(aux)) = aux ;
  end
else
  error('argumento not cell')
end
