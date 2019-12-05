% run all examples from examples folder.

fileslist = readdir('../examples');

for ind = 1:length(fileslist)

  if length( strfind( fileslist{ind},'TEST') ) > 0  
    run( [ '../examples/' fileslist{ind} ] ) ;
    pause(1)
  end
  
end
