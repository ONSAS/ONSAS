function sep = dirSep()
if isunix,
  sep = '/';
else 
  sep = '\';
end