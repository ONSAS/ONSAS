
%
% export_ONSAS_docs_example('staticVonMisesTruss/onsasExample_staticVonMisesTruss.m' ,'~/.julia/dev/ONSAS_docs/docs/src/tutorials/StaticVonMisesTruss/staticVonMisesTruss.md' )

function export_ONSAS_docs_example( fileIn, fileOut )

fidIn  = fopen( fileIn, 'r' );
fidOut = fopen( fileOut,'w' );
isInCodeBlock = false ;

while ~feof(fidIn)
  currentLine = fgetl( fidIn ) ;

  if length(currentLine)>=2 && strcmp( currentLine(1:2), '%#' )

    if isInCodeBlock % closes code block before writing comment
      fprintf( fidOut,'```\n' );
      isInCodeBlock = false ;
    end
    fprintf( fidOut,'%s\n', currentLine(3:end) );

  else
    if ~isInCodeBlock % open code block
      fprintf( fidOut,'```\n' );
      isInCodeBlock = true ;
    end
    fprintf( fidOut,'%s\n', currentLine );
  end  

end

fclose(fidIn);
fclose(fidOut);



