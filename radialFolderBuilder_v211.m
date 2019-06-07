%% radialFolderBuilder_v211.m
% This function builds the folder structure for storing the converted nc files
% for radial data and the full filename of the converted nc radial file,
% according to the structure YYYY/YYYY_MM/YYYY_MM_DD/.

% INPUT:
%         mainPath: root path for total nc folder.
%         siteCode: code of the radial site.
%         ts: timestamp of the total file.

% OUTPUT:
%         rFB_err: error flag (0 = correct, 1 = error)
%         fullPath: full path of the folder where to store the converted nc
%                   radial file.


% Author: Lorenzo Corgnati
% Date: May 30, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [rFB_err,fullPath] = radialFolderBuilder_v211(mainPath,siteCode,ts)

disp(['[' datestr(now) '] - - ' 'radialFolderBuilder_v211.m started.']);

rFB_err = 0;

fullPath = '';

warning('off', 'all');

% Retrieve the time-related elements of the path
try
    [yMDF_err,yearFolder,monthFolder,dayFolder] = yearMonthDayFolder(ts);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rFB_err = 1;
    return
end

% Create the folder structure
% Version 2.1.1
try
    if (exist([mainPath filesep 'v2.1.1' filesep siteCode filesep yearFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.1.1' filesep siteCode filesep yearFolder]);
    end
    
    if (exist([mainPath filesep 'v2.1.1' filesep siteCode filesep monthFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.1.1' filesep siteCode filesep monthFolder]);
    end
    
    if (exist([mainPath filesep 'v2.1.1' filesep siteCode filesep dayFolder], 'dir') ~= 7)
        mkdir([mainPath filesep 'v2.1.1' filesep siteCode filesep dayFolder]);
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rFB_err = 1;
    return
end

% Build the full filename for the nc total file
fullPath = [mainPath filesep 'v2.1.1' filesep siteCode filesep dayFolder];

if(rFB_err==0)
    disp(['[' datestr(now) '] - - ' 'radialFolderBuilder_v211.m successfully executed.']);
end

return