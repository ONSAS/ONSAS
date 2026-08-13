function m2md( fileIn, fileOut, includeCodeBoolean, iniLine )

fidIn  = fopen( fileIn, 'r' );
fidOut = fopen( fileOut,'w' );

for i=1:iniLine
  currentLine = fgetl( fidIn ) ;
end

isInCodeBlock = ~( length( currentLine )>=3 && strcmp( currentLine(1:3), '%md' ) ) ;
lineCount = 0 ;

while ~feof( fidIn )

  lineCount = lineCount + 1;
  if lineCount ~= 1
    currentLine = fgetl( fidIn ) ;
  end

  if length( currentLine )>=3 && strcmp( currentLine(1:3), '%md' ) % not code

    if isInCodeBlock % closes code block before writing comment
      if includeCodeBoolean
        fprintf( fidOut,'```\n' );
      end
      isInCodeBlock = false ;
    end
    fprintf( fidOut,'%s\n', currentLine(4:end) );

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
