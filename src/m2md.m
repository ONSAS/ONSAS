function m2md( fileIn, fileOut, includeCodeBoolean )

fidIn  = fopen( fileIn, 'r' );
fidOut = fopen( fileOut,'w' );

currentLine = fgetl( fidIn ) ;
isInCodeBlock = ~( length( currentLine )>=2 && strcmp( currentLine(1:2), '%#' ) ) ;
lineCount = 0 ;

while ~feof( fidIn )
  
  lineCount = lineCount + 1;
  if lineCount != 1
    currentLine = fgetl( fidIn ) ;
  end
  
  if length( currentLine )>=2 && strcmp( currentLine(1:2), '%#' ) % not code

    if isInCodeBlock % closes code block before writing comment
      if includeCodeBoolean
        fprintf( fidOut,'```\n' );
      end
      isInCodeBlock = false ;
    end
    fprintf( fidOut,'%s\n', currentLine(3:end) );

  else
    if ~isInCodeBlock % open code block
      if includeCodeBoolean
        fprintf( fidOut,'```\n' );
      end
      isInCodeBlock = true ;
    end
    if includeCodeBoolean
      fprintf( fidOut,'%s\n', currentLine );
    end
  end  
end

fclose(fidIn);
fclose(fidOut);
