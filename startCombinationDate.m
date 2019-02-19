%% startCombinationDate.m
% This function builds the datetime of the starting day of the processing
% period in the MySQL format yyyy-mm-dd from the present timestamp.

% INPUT:
%         tsNow: present timestamp in numeric format

% OUTPUT:
%         sCD_err: error flag (0 = correct, 1 = error)
%         startDatetime: string for datetime in the MySQL format yyyy-mm-dd


% Author: Lorenzo Corgnati
% Date: July 18, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function startDatetime = startCombinationDate(tsNow)

disp(['[' datestr(now) '] - - ' 'startCombinationDate.m started.']);

sCD_err = 0;

warning('off', 'all');

% Evaluate the starting date as 8 days ago
try
    tsStart = tsNow - 8;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sCD_err = 1;
    return
end

try
    % Build the string for datetime in the MySQL format YYYY-MM-DD
    formatOut='yyyy-mm-dd';
    startDatetime = datestr(tsStart,formatOut);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sCD_err = 1;
end

if(sCD_err==0)
    disp(['[' datestr(now) '] - - ' 'startCombinationDate.m succesfully executed.']);
end

return

