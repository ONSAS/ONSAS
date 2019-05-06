
for previouslyDefinedSelectedFileVar = 1:length(readdir('./input'))
  ind = previouslyDefinedSelectedFileVar ;
  
  fileslist = readdir('./input');

  if length( strfind( fileslist{ind},'TEST') ) > 0  
    previouslyDefinedSelectedFileVar = ind ;
    ONSAS
  
    if noErrorsOccurred == 1
      fprintf('No errors occurred. Press any key to continue...');
      pause,
      fprintf('\n\n'); clear all, close all
    else
      error('Errors occurred');
    end
  end
end
