%% UScompareTimeStamps.m
% This function compares the timestamps of the total files already
% converted into the standard format for HFR-US networks with the
% timestamps read from HFR-US TDS catalogs, selects the ones yet to 
% be converted and returns their indices related to the netCDF data
% variable.

% INPUT:
%         DBdatetime: timestamps of the total files already converted into 
%                     the standard format for HFR-US network
%         OpenDAPdatetime: timestamps read from HFR-US TDS catalog

% OUTPUT:
%         UcT_err: error flag (0 = correct, 1 = error)
%         toBeConvertedTS: indices of data from HFR-US network to be converted


% Author: Lorenzo Corgnati
% Date: April 10, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [UcT_err, toBeConvertedTS] = UScompareTimeStamps(DBdatetime,OpenDAPdatetime)

disp(['[' datestr(now) '] - - ' 'UScompareTimeStamps.m started.']);

UcT_err = 0;

warning('off', 'all');

%%

try
    % Find the indices of the timestamps already converted
    [~,idx] = ismember(OpenDAPdatetime,datenum(DBdatetime));
    % Retrieve the indices of the timestamps to be converted
    toBeConvertedTS = find(idx==0);    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    UcT_err = 1;
end

%%

if(UcT_err==0)
    disp(['[' datestr(now) '] - - ' 'UScompareTimeStamps.m successfully executed.']);
end

return

