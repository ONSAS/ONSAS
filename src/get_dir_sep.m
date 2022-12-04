function dirSep = get_dir_sep()
if isunix
  dirSep = '/';
else
  dirSep = '\';
end
