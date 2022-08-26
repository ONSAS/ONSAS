% Recursive function for updating headers of source files in src and downwards ...

function updateLicenseHeaders( folder )

folder
files = dir( folder ) ;

for i = 1:length( files )
  i
  if files(i).isdir
    if ~strcmp(files(i).name(1),'.')
      files(i).name
      updateLicenseHeaders( [ folder '/' files(i).name ] )
    end 
  else
    completeFilename = [ folder '/' files(i).name ] 
    showHeaderAndReplace( completeFilename )
  end
end


function showHeaderAndReplace( filename )

lengthOfCurrentHeader = 17 ;

system( [ 'head -n ' num2str(lengthOfCurrentHeader) ' ' filename ] );

reply = input('    replace yes or no? (y/n):','s' ) ;

if strcmp( reply, 'y')
  system( [ 'more currentFileHeader.txt > aux.txt' ] ) ;
  system( [ 'tail  -n +' num2str(lengthOfCurrentHeader+1) ' ' filename ' >> aux.txt'] ) ;
  system( [ 'mv aux.txt ' filename ] ) ;
  disp('file updated.')
end
pause
