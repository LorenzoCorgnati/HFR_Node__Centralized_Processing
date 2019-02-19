%% yearMonthDayFolder.m
% This function builds the strings for year folder, month folder and day folder
% according to the folder structure YYYY/YYYY_MM/ YYYY_MM _DD.


% INPUT:
%         ts: timestamp

% OUTPUT:
%         yMDF_err: error flag (0 = correct, 1 = error)
%         yearFolder: string for year folder
%         monthFolder: string for month folder
%         dayFolder: string for day folder


% Author: Lorenzo Corgnati
% Date: July 3, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [yMDF_err,yearFolder,monthFolder,dayFolder] = yearMonthDayFolder(ts)

disp(['[' datestr(now) '] - - ' 'yearMonthDayFolder.m started.']);

yMDF_err = 0;

warning('off', 'all');

% Retrieve the time-related elements of the path
try
    tsCell = strsplit(ts);
    tsVec = [];
    for tsCell_idx=1:length(tsCell)
        tsVec = [tsVec str2double(tsCell{tsCell_idx})];
    end
    tsNum = datenum(tsVec);
    [year,month,day,hr,minutes,sec] = datevec(tsNum);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    yMDF_err = 1;
    return
end

try
    % Build the folder strings
    yearFolder = num2str(year);
    monthFolder = [yearFolder filesep num2str(year) '_' num2str(month,'%02d')];
    dayFolder = [monthFolder filesep num2str(year) '_' num2str(month,'%02d') '_' num2str(day,'%02d')];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    yMDF_err = 1;
end

if(yMDF_err==0)
    disp(['[' datestr(now) '] - - ' 'yearMonthDayFolder.m successfully executed.']);
end

return

