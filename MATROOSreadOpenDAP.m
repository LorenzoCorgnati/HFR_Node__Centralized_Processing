%% MATROOSreadOpenDAP.m
% This function reads data from the HFR-MATROOS network via OpenDAP and
% selects the data subset to be converted according to the processing time
% interval.

% INPUT:
%         procStart: processing start date
%         networkData: cell array containing information about the network
%                      (metadata)
%         networkFields: field names of the cell array containing
%                       information about the network.

% OUTPUT:
%         MrO_err: error flag (0 = correct, 1 = error)
%         nc: data from HFR-MATROOS network to be converted


% Author: Lorenzo Corgnati
% Date: January 11, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [MrO_err, nc] = MATROOSreadOpenDAP(procStart,networkData,networkFields)

disp(['[' datestr(now) '] - - ' 'MATROOSreadOpenDAP.m started.']);

MrO_err = 0;

warning('off', 'all');

%% Find the processing time interval

try
    % Find the TDS_root_url field from network data
    TDS_root_urlIndex = find(not(cellfun('isempty', strfind(networkFields, 'TDS_root_url'))));
    TDS_root_url = networkData{TDS_root_urlIndex};
    
    % Read time and convert it to Matlab time
    nc.time = ncread_cf_time(TDS_root_url,'time');
    
    % Select the data range according to the processing time span
    iTime = find(nc.time>=datenum(procStart));
    nc.time = nc.time(iTime);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    MrO_err = 1;
end

%%

%% Read variables

try
    % Coordinate variables
    nc.latitude = ncread(TDS_root_url,'latitude');
    nc.longitude = ncread(TDS_root_url,'longitude');
    nc.depth = 0;
    
    % Data variables
    nc.ewct = ncread(TDS_root_url,'ewct',[1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(iTime)]);
    nc.nsct = ncread(TDS_root_url,'nsct',[1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(iTime)]);
    nc.uacc = ncread(TDS_root_url,'ewct_error',[1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(iTime)]);
    nc.vacc = ncread(TDS_root_url,'nsct_error',[1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(iTime)]);
    nc.gdopX = ncread(TDS_root_url,'gdopx');
    nc.gdopY = ncread(TDS_root_url,'gdopy');
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    MrO_err = 1;
end

%%

if(MrO_err==0)
    disp(['[' datestr(now) '] - - ' 'MATROOSreadOpenDAP.m successfully executed.']);
end

return

