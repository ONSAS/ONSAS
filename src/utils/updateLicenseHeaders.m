% Recursive function for updating headers of source files in src and examples

function updateLicenseHeaders(folder, root)

lengthOfCurrentHeader = 17; % Define expected header length

% Only generate newFileHeader.txt once when root is true
if root == true
    if exist('newFileHeader.txt', 'file') ~= 2
        system(['head -n ' num2str(lengthOfCurrentHeader) ' updateLicenseHeaders.m > newFileHeader.txt']);
    end
    % Process both 'src' and 'examples' from root
    updateLicenseHeaders('../src', false);

    % Remove the temporary file at the end
    if exist('newFileHeader.txt', 'file') == 2
        system('rm newFileHeader.txt');
    end

    return; % Exit after processing both main folders
end

% Process the given folder recursively
disp(['Processing folder: ', folder]);
files = dir(folder);

for i = 1:length(files)
    if files(i).isdir
        % Skip hidden folders like .git
        if ~strcmp(files(i).name(1), '.')
            subfolder = fullfile(folder, files(i).name);
            disp(['Entering directory: ', subfolder]);
            updateLicenseHeaders(subfolder, false);
        end
    else
        completeFilename = fullfile(folder, files(i).name);
        if endsWith(completeFilename, '.m') % Only process .m files
            showHeaderAndReplace(completeFilename, lengthOfCurrentHeader);
        end
    end
end
end

function showHeaderAndReplace(filename, lengthOfCurrentHeader)

% Display current file header
system(['head -n ' num2str(lengthOfCurrentHeader) ' ' filename]);

% Ask user whether to replace the header
reply = input('Replace header? (y/n): ', 's');

if strcmp(reply, 'y')
    % Create a new file with the updated header
    system('cp newFileHeader.txt aux.txt');
    system(['tail -n +' num2str(lengthOfCurrentHeader + 1) ' ' filename ' >> aux.txt']);
    system(['mv aux.txt ' filename]);
    disp(['Updated: ', filename]);
else
    disp(['Skipped: ', filename]);
end

pause; % Wait before moving to the next file
end
