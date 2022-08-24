%script for extracting variables from struct data sets.

arrfie    = cell(6,1);
arrfie{1} = modelCurrSol;
arrfie{2} = BCsData ;
arrfie{3} = modelProperties  ;
arrfie{4} = 'modelCurrSol';
arrfie{5} = 'BCsData' ;
arrfie{6} = 'modelProperties'  ;

for k=1:3
  fieldNames = fieldnames( arrfie{k} ) ;
  for i=1:length(fieldNames)
    eval( [ fieldNames{i} '= getfield(' arrfie{k+3} ',''' fieldNames{i} ''');' ] )
  end
end
