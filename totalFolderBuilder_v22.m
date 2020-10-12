%% totalFolderBuilder_v22.m
% This function builds the folder structure for storing the converted nc files
% for total data and the full filename of the converted nc total file,
% according to the structure YYYY/YYYY_MM/YYYY_MM_DD/.

% INPUT:
%         mainPath: root path for total nc folder.
%         ts: timestamp of the total file

% OUTPUT:
%         tFB_err: error flag (0 = correct, 1 = error)
%         fullPath: full path of the folder where to store the converted nc
%                   total file.


% Author: Lorenzo Corgnati
% Date: September 28, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [tFB_err,fullPath] = totalFolderBuilder_v22(mainPath,ts)

disp(['[' datestr(now) '] - - ' 'totalFolderBuilder_v22.m started.']);

tFB_err = 0;

fullPath = '';

warning('off', 'all');

% Retrieve the time-related elements of the path
try
    [yMDF_err,yearFolder,monthFolder,dayFolder] = yearMonthDayFolder(ts);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tFB_err = 1;
    return
end

% Create the folder structure
% Version v2.2
try
    if (exist([mainPath filesep 'v2.2' filesep yearFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.2' filesep yearFolder]);
    end
    
    if (exist([mainPath filesep 'v2.2' filesep monthFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.2' filesep monthFolder]);
    end
    
    if (exist([mainPath filesep 'v2.2' filesep dayFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.2' filesep dayFolder]);
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tFB_err = 1;
    return
end

% Build the full filename for the nc total file
fullPath = [mainPath filesep 'v2.2' filesep dayFolder];

if(tFB_err==0)
    disp(['[' datestr(now) '] - - ' 'totalFolderBuilder_v22.m successfully executed.']);
end

return

