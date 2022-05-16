function addcopyright(target,update,year_or_string,organization,prompt)
%  ADD COPYRIGHT info to M-file.
%  Syntax: ADDCOPYRIGHT(TARGET,UPDATE,YEAR_OR_STRING,ORG,INFODISP)  
%   ADDCOPYRIGHT adds the copyright info to an M-file or all M-files in a 
%   designated folder as well as subfolders.  
%
%   The copyright line that is built up looks like the following:
%         %   Copyright <YEAR> <ORGANIZATION>
%   As for example:
%         %   Copyright 2003 The MathWorks, Inc.
%
%   The copyright info is place directly above the Revision line (containing
%   "$Revision") or after the first blank line that is not commented after
%   the first commented section. 
%
%   TARGET is either an M-file name or a folder name. 
%       If empty (default) a folder dialog box is opened.
%       If 0 or ' ' then the current directory is used.
%       If no extension is provided, the program looks for TARGET.M before
%           checking if TARGET is a folder.
%   UPDATE is either Boolean (TRUE,FALSE), a case insensitive string that
%       uniquely denotes one of {true, false, replace}, or a number {1, 0, -1}.
%       If copyright info does not exist, it is added irregarless of UPDATE.
%       If UPDATE is FALSE or 0, copyright info is not modified. (Default)
%       If UPDATE is TRUE or 1 then the copyright info is modified.
%           If the inputed year, YEAR2, is older than the previous year (YEAR1),
%           or older than the 1st year in a range of years, then the copyright
%           year has the format: YEAR1-YEAR2.
%       If UPDATE is REPLACE or -1 then the copyright year or range of
%       years is replaced by the inputed year or range of years.
%   YEAR_OR_STRING is either the year as a number, year or range of years as
%       a string (e.g. '2000-2003') or a string replacement, COPYRIGHT_STRING,
%       for the entire copyright line. See UPDATE.
%       If empty, it defaults to the current year.
%       If blank & UPDATE is TRUE, then the copyright line is removed.
%   ORG is a string with the organization's name. (Optional) 
%       Edit DEF_ORGANIZATION to put in your own organization as default.
%   INFODISP is a string argument for bypassing screen display and the user confirmation prompt.
%       If '-noconfirm' then the confirmation prompt is not displayed.
%       If '-suppress' then all screen display and confirmation is suppressed.    
%   The order of ORG and INFODISP may be interchanged.
%
%  EXAMPLES:
%   Note that the command line syntax is supported!
%   ADDCOPYRIGHT(FOLDER) 
%   Adds the copyright info to all M-files in the folder FOLDER & its subfolders.
%   If copyright info exists it is not modified. 
%   If FOLDER is an M-file then adds copyright info only to the file.
%   The copyright year is the current year.
%   The organization is set to its default value.
%
%   ADDCOPYRIGHT(FOLDER,[],COPYRIGHT_YEAR) 
%   The copyright year is set to COPYRIGHT_YEAR if no copyright exists.
%   COPYRIGHT_YEAR may either be a character string (e.g. '2003', '2001-2003',
%   or '2001:2003') or a numeric (e.g. 2006).
%
%   ADDCOPYRIGHT '' T 
%   Calls up a dialog box to select the directory and
%   updates the copyright info in the selected folder and all subfolders.
%
%   ADDCOPYRIGHT ' ' T 
%   Updates the copyright in the current folder & subfolders.
%
%   ADDCOPYRIGHT(FOLDER,UPDATE,COPYRIGHT_STRING)
%   Uses COPYRIGHT_STRING instead of '   Copyright <YEAR> <ORGANIZATION>'.
%   This can be used to replace all instances of the copyright info.
%
%   ADDCOPYRIGHT ' ' T ' ' (note the blank space in the second quotes)
%   Deletes the copyright line in the current folder.
%
%   ADDCOPYRIGHT ' ' R 2002-2006
%   Using the current folder, replaces the year in the copyright info with 
%   the range, 2002-2006, ignoring all prior year information.
%
%   ADDCOPYRIGHT FILENAME T 2005 -suppress
%   Update the year in the copyright info in FILENAME. 
%   If the prior year is 2001 then the new year will be 2001-2005.
%   If the prior year is 2002-2003 then the new year will be 2002-2005.
%   All screen output and confirmation prompt is suppressed.
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
% Copyright 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. 

% This file is part of ONSAS 

% ONSAS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 

% ONSAS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
%   Copyright 2006-2008 Mirtech, Inc.
%   Copyright 2003 The MathWorks, Inc.
%   $Revision: 1.2 $ $Date: 2003/09/11 21:25:11 $
%   Author: Raymond Norris (rayn@mathworks.com)
%   Revised 10/05/2006  Mirko Hrovat (mhrovat@email.com) 
%   Revised 08/18/2009  Mirko Hrovat
%       update eol character to be based upon os instead of computer type.
%       changed '' to be equivalent to empty, and ' ' to be the current directory
%       allowed screen display to be suppressed
%       added routine getfilenames, so that a list of all m-files is created before modification
%
%   NOTES:
%   The following summarizes some of the differences between the
%       original version by Raymond Norris and this version.
%   1. Argument list order is changed, and number of arguments is reduced.
%   2. Command line syntax is supported.
%   3. If folder argument is empty, a dialog box selects the folder.
%   4. Copyright info may be added to a single file.
%   5. Copyright info is added either before $Revision keyword or after a
%       blank line.
%   6. Number of spaces between % and Copyright is immaterial.
%   7. The organization is checked and if different a new copyright line is
%       added.
%   8. A Def_organization constant is added.
%   9. Will not attempt to add copyright info to contents.m or files that
%       are not M-files.
%   10. If copyright_string is blank and update is true then copyright line
%       is deleted.
%   11. The correct eol character for the OS is used to rewrite the files.
%   12. More file checks are performed to prevent leaving temporary files.
%   13. Summary output has been changed to be easier to read.
%   14. Added features to work with range of years as specified by year1-year2.
% ------- CONSTANTS -------
% Modify def_organization to put in your orgnaization.
%def_organization   = 'The MathWorks, Inc.';
def_organization    = "2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, \n \
% Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini. \
\n\n \
% \This file is part of ONSAS \
\n\n \
% ONSAS is free software: you can redistribute it and/or modify \n \
% it under the terms of the GNU General Public License as published by \n \
% the Free Software Foundation, either version 3 of the License, or \n \
% (at your option) any later version. \
\n\n \
% ONSAS is distributed in the hope that it will be useful, \n \
% but WITHOUT ANY WARRANTY; without even the implied warranty of\n \
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n \
% GNU General Public License for more details.\
\n\n \
% You should have received a copy of the GNU General Public License\n \
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.\n " ;
confirm         = true;         % set to false to bypass confirmation prompt.
scrndisp        = true;         % set to false to stop all screen display
% Parse input arguments
switch nargin,
    case 5
    case 4
        prompt=[];
    case 3
        prompt=[];  organization=[];
    case 2
        prompt=[];  organization=[];    year_or_string=[];
    case 1
        prompt=[];  organization=[];    year_or_string=[];  update=[];
    case 0
        prompt=[];  organization=[];    year_or_string=[];  update=[];  target=[];
    otherwise
        error ('  Too many arguments!');
end
switch true,
    case strcmpi(prompt,'-noconfirm')
        confirm = false;
    case strcmpi(prompt,'-suppress')
        scrndisp = false;
        confirm = false;
    case strcmpi(organization,'-noconfirm')
        confirm = false;
        organization = prompt;
    case strcmpi(organization,'-suppress')
        scrndisp = false;
        confirm = false;
        organization = prompt;
    otherwise
        % do nothing
end
if isempty(organization),   organization=def_organization;      end
if isempty(update),         update=0;                           end
if ischar(update),    
    if any(strcmpi(update,{'t','tr','tru','true'})),
        update = 1;
    elseif any(strcmpi(update,{'r','re','rep','relp','repla','replac','replace'})),
        update = -1;
    else
        update = 0;
    end
elseif islogical(update),
    update=double(update);
elseif all(update ~= [-1,0,1]),
    update = 0;
end
% check if year_or_string is an expression or a year
year_flg = 0;
if isempty(year_or_string),
    year_or_string=datestr(date,10);
elseif isnumeric(year_or_string),
    year_or_string=num2str(year_or_string); % it is a number, convert to string
elseif isempty(str2num(year_or_string)),    %#ok
% in this case year_or_string remains unchanged, just make sure it is not numeric
    year_flg = 3;   % not numeric, use for UPDATE to indicate replacement of entire string
end
if isempty(target),     % have user select the folder
    target = uigetdir('','Select directory for adding Copyright stamp.');
    if isempty(target)||isequal(target,0),
        error('  No directory folder was selected!')
    end
end
curdir = pwd;
if (~isempty(target)&&all(target==' ')) || all(target==0),
    target = curdir;
end
switch true
    case exist(target,'file')==2,
        % target or target.m is a file, just do the file
        status = ac2f(target,update+year_flg,year_or_string,organization);
        if scrndisp,
            fprintf('  ... %s\n',status)
        end
    case exist(target,'dir')==7,
        % target is a folder do all files in the folder and subfolders
        % Use "what" to get full path name to folder, since "dir" does not use
        % matlab search paths. Do use "dir" to see everything in the folder.
        whatlst = what(target);
        folder = whatlst(1).path;
        if scrndisp,
            fprintf(' Starting Path: %s \n',folder)
        end
        if confirm,
            ok = input('  Please confirm, is this folder ok? (y/n)[y] :', 's');
            if isempty(ok),     ok = 'y';       end
        else
            ok = 'y';
        end
        if strcmpi(ok,'y'),
            mfiles = getfilenames(folder);  % get all mfiles and paths
            oldfolder = [];
            for n = 1:length(mfiles)
                curfolder = mfiles(n).path;
                if ~isempty(curfolder) && ~strcmp(curfolder,oldfolder),
                    oldfolder = curfolder;
                    if scrndisp,
                        fprintf('\n    Directory: %s\n',curfolder)
                    end
                end
                fname = fullfile(folder,curfolder,mfiles(n).name);
                status = ac2f(fname,update+year_flg,year_or_string,organization);
                if scrndisp,
                    fprintf('  %-20s... %s\n',mfiles(n).name,status)
                end
            end
        end
    otherwise
        fprintf('\n *** File or folder %s does not exist.\n',target)
end
end  % ----------------- AddCopyRight ---------------------
% ----------------- getfilenames -------------------
function filelist = getfilenames(dirname)
% gets the names of all of the m-files in the directory and subdirectories
% filelist is a struct array with the fields:
%   .name - name of the file
%   .path - relative path for the file, starting after "dirname"
dlist = dir(dirname);
% Order the list with files first then directories.
idx = [dlist.isdir];
newdlist = [dlist(~idx); dlist(idx)];
filelist = struct('name',{},'path',{});
for n = 1:length(newdlist),
    dname = newdlist(n);
    switch true
        case strcmp(dname.name,'.') || strcmp(dname.name,'..') || strcmp(dname.name,'CVS'),
            % do nothing
        % found directory, call getfilenames again
        case dname.isdir==true
            nextfolder = fullfile(dirname,dname.name);
            subdirlist = getfilenames(nextfolder);
            % correct patch by adding subdirectory
            for k = 1:length(subdirlist),
                subdirlist(k).path = fullfile(dname.name,subdirlist(k).path);
            end
            if ~isempty(subdirlist),
                filelist = [filelist,subdirlist];       %#ok
            end
        case length(dname.name)>2 && strcmp(dname.name(end-1:end),'.m')...
                && ~strcmpi(dname.name,'contents.m') 
            filelist(end+1).name = dname.name;          %#ok
            filelist(end).path = [];                  
        otherwise
            % do nothing
    end
end
end% ----------------- getfilenames -------------------
% ----------------- ac2f -------------------
function statstr = ac2f(fname,update,year_or_string,organization)
% Add Copyright to(2) File
% This function searches line by line for the copyright info and writes it
%   to the file pointed to by fname.
% This function was modified to insure that the organization matches.
% It also supports multiyear formats.
% If the organization does not match then a new Copyright line is
%   prepended.
% Update has 6 possibilites, first 3 are with a valid year date
%       second 3 are with a replacement copyright string
%   -1: Replace year        0: Add only     1: Update year
%    2: Replace string      3: Add only     4: Replace string
% A string is returned that describes the action that was taken on fname.
% define platform dependent eol
cr=13;
lf=10;
switch true
    case ispc
        eol = char([cr,lf]);
    case isunix
        eol = char(lf);
    otherwise       % if not a pc or unix os, assume it is an old mac os
        eol = char(cr);
end
[filepath,file,ext] = fileparts(fname);
% need to check for blank extensions since "exist" in the calling function
% will find M-files without the extension explicitly given.
if strcmp(ext,''),
    ext='.m';
end
sfname1 = [file,ext];
% skip over any file that is not an M-file or is "contents.m".
if strcmpi(file,'contents') || ~strcmpi(ext,'.m'),
    return          
end
filename1 = fullfile(filepath,sfname1);
fid_1 = fopen(filename1,'r+');
if fid_1<3
    statstr = sprintf('*** Failed to open %s for copyrighting.',sfname1);
    return
end
sfname2 = ['new_' sfname1];
filename2 = fullfile(filepath,sfname2);
fid_2 = fopen(filename2,'w');
if fid_2<3
    statstr = sprintf('*** Failed to open %s for writing.',sfname2);
    fclose(fid_1)
    return
end
blank_copyright_string  = strcmp(deblank(year_or_string),'');
found_token             = false;
token_already_embedded  = false;
first_blank_line        = false;
another_org             = false;
comment_section         = false;
next_line = fgetl(fid_1);
while ischar(next_line)
    if ~found_token
        % modified next lines to make case insensitive,
        % ignore spaces between % and Copyright, look for inclusion of date,
        % & check the organization
        if regexpi(strrep(next_line,' ',''),'%Copyright\d\d\d\d','ONCE')==1,
            % found Copyright line
            found_token = true;
            if isempty(strfind(next_line,organization)),
                % organization does not match, prepend copyright info
                WriteCopyrightLine(year_or_string)
                another_org = true;
            else
                % organization matches
                token_already_embedded = true;
                switch update,
                    case {0,3}
                        % no update, do nothing.
                    case 1
                        % update but check year range
                        year1=str2double(regexp(next_line,'\d\d\d\d','match','once'));
                        temp=regexp(year_or_string,'\d\d\d\d','match');
                        year2=str2double(temp{end});
                        if year2>year1
                            WriteCopyrightLine([num2str(year1),'-',num2str(year2)])
                        else
                            WriteCopyrightLine(year_or_string)
                        end
                        next_line = fgetl(fid_1);
                        continue
                    case {-1,2,4}
                        % for 2 & 4 update but replace entire string
                        % for -1 simply replace year or year range
                        WriteCopyrightLine(year_or_string)
                        next_line = fgetl(fid_1);
                        continue
                end  % switch
            end  % if isempty(strfind(next_line,organization)),
        elseif ~found_token && (strncmp(next_line,'%   $Revision',13)||first_blank_line),
                WriteCopyrightLine(year_or_string)
                found_token = true;
        elseif ~first_blank_line && strncmp(next_line,'%',1),
            comment_section = true;
        elseif strcmp(deblank(next_line),'') && comment_section,
            first_blank_line = true;
            comment_section = false;
        end
    end  % if ~found_token
    fprintf(fid_2,'%s',[next_line,eol]);    % copy fid_1 into fid_2
    next_line = fgetl(fid_1);
end
fclose(fid_1);
fclose(fid_2);
% set output string and cleanup
switch true
    case ~found_token
        statstr = 'Could not find place to add Copyright info.';
        delete(filename2)
    case token_already_embedded==true && any(update==[0,3])
        statstr = 'Copyright info already found. File not updated.';
        delete(filename2)
    case blank_copyright_string
        if token_already_embedded==true && any(update==[-1,1,2,4])
            moved=movefile(filename2,filename1,'f');
            if ~moved, 
                delete(filename2)
                statstr = 'Could not move file!!! Copyright info is unchanged.';
            else
                statstr = 'Copyright info removed.';
            end 
        else
            statstr = 'File not updated. Blank Copyright line.';
            delete(filename2)
        end
    case another_org   
        moved = movefile(filename2,filename1,'f');
        if ~moved, 
            delete(filename2)
            statstr = 'Could not move file!!! Copyright info is unchanged.';
        else
            statstr = 'Another organization found! Copyright info added.';
        end
    otherwise
        moved = movefile(filename2,filename1,'f');
        if ~moved, 
            delete(filename2)
            statstr = 'Could not move file!!! Copyright info is unchanged.';
        else
            statstr = 'Copyright info added/modified.';
        end
end
    % ----------------- WriteCopyrightLine -------------------
    function WriteCopyrightLine(year_or_string)
    % This Nested! Function writes the copyright line to the file
        if ~blank_copyright_string
            if update>1;
                copyright_string = ['%   ' year_or_string];
            else
                copyright_string=['%   Copyright ' year_or_string ' ' organization];
            end
            fprintf(fid_2,'%s',[copyright_string,eol]);
        end
    end  % ----------------- WriteCopyrightLine -------------------
end  % ----------------- ac2f -------------------
