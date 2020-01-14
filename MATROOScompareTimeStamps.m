%% MATROOScompareTimeStamps.m
% This function compares the timestamps of the total files already
% converted into the standard format for HFR-MATROOS network with the
% timestamps read from HFR-MATROOS TDS catalog, selects the ones yet to 
% be converted and returns their indices related to the netCDF data
% variable.

% INPUT:
%         DBdatetime: timestamps of the total files already converted into 
%                     the standard format for HFR-MATROOS network
%         OpenDAPdatetime: timestamps read from HFR-MATROOS TDS catalog

% OUTPUT:
%         McT_err: error flag (0 = correct, 1 = error)
%         toBeConvertedTS: indices of data from HFR-MATROOS network to be converted


% Author: Lorenzo Corgnati
% Date: January 11, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [McT_err, toBeConvertedTS] = MATROOScompareTimeStamps(DBdatetime,OpenDAPdatetime)

disp(['[' datestr(now) '] - - ' 'MATROOScompareTimeStamps.m started.']);

McT_err = 0;

warning('off', 'all');

%%

try
    % Find the indices of the timestamps already converted
    [~,idx] = ismember(OpenDAPdatetime,datenum(DBdatetime));
    % Retrieve the indices of the timestamps to be converted
    toBeConvertedTS = find(idx==0);    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    McT_err = 1;
end

%%

if(McT_err==0)
    disp(['[' datestr(now) '] - - ' 'MATROOScompareTimeStamps.m successfully executed.']);
end

return

