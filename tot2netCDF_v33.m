%% tot2netCDF_v33.m
% This function converts hourly total files combined starting from radial
% files in netCDF-4 format according to the European common data and metadata
% model integrating CMEMS-INSTAC and SDC CF extension requirements.
% v3.2 uses scale_factor and add_offset for packing geophysical variables.
% v3.3 names depth dimension as DEPTH; don not apply _SEADATANET_ string in
% TIME, POSITION and DEPTH QC variables; uses blank space as separator in
% institution_edmo_code global attribute.

% INPUT:
%         mat_tot: structure containing radial file in Codar format
%         networkData: cell array containing information about the network
%                      (metadata)
%         networkFields: field names of the cell array containing
%                       information about the network.
%         stationData: cell array containing information about the contributing
%                      stations (metadata)
%         stationFields: field names of the cell array containing
%                       information about the station.
%         timestamp: timestamp of the total file to be converted.
%         institutions: cell array containing information about all the stations.


% OUTPUT:
%         T2C_err: error flag (0 = correct, 1 = error)
%         networkData: cell array containing information about the network.
%                       It is returned in case an update of the database is
%                       needed.
%         stationData: cell array containing information about the contributing
%                       stations.
%                       It is returned in case an update of the database is
%                       needed.
%         ncFileNoPath: filename of the converted nc file, without the full path
%         ncFilesize: size of the converted nc file.


% Author: Lorenzo Corgnati
% Date: November 9, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [T2C_err,networkData,stationData,ncFileNoPath,ncFilesize] = tot2netCDF_v33(mat_tot, networkData, networkFields, stationData, stationFields,timestamp,institutions)

disp(['[' datestr(now) '] - - ' 'tot2netCDF_v33.m started.']);

T2C_err = 0;

warning('off', 'all');

%% Set scale_factor and add_offset values

scaleFactor = 0.001;
addOffset = 0;

%%

%% Retrieve site codes and coordinates

try
    % Find the index of the site_lon field
    site_lonIndexC = strfind(stationFields, 'site_lon');
    site_lonIndex = find(not(cellfun('isempty', site_lonIndexC)));
    
    % Find the index of the site_lat field
    site_latIndexC = strfind(stationFields, 'site_lat');
    site_latIndex = find(not(cellfun('isempty', site_latIndexC)));
    
    % Find the index of the station_id field
    station_idIndexC = strfind(stationFields, 'station_id');
    station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
    
    % Build the arrays for site coordinates and codes
    sitesLon = cell2mat(stationData(:,site_lonIndex));
    sitesLat = cell2mat(stationData(:,site_latIndex));
    sitesCodes = cell2mat(stationData(:,station_idIndex));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Retrieve EDMO codes and institution names

try
    % Find the EDMO_code field from network data
    NT_EDMO_codeIndex = find(not(cellfun('isempty', strfind(networkFields, 'EDMO_code'))));
    NT_EDMO_code = networkData{NT_EDMO_codeIndex};
    
    % Find the EDMO_code field from station data
    ST_EDMO_codeIndex = find(not(cellfun('isempty', strfind(stationFields, 'EDMO_code'))));
    ST_EDMO_code = cell2mat(institutions(:,ST_EDMO_codeIndex));
    ST_EDMO_code = ST_EDMO_code(ST_EDMO_code~=0);
    
    % Build the cumulative EDMO code list
    [EDMO_code,ia,ic] = unique([NT_EDMO_code; ST_EDMO_code]);
    EDMO_code = EDMO_code';
    EDMO_codeStr = sprintf('%.0d ' , EDMO_code);
%     EDMO_codeStr = EDMO_codeStr(1:end-2);% strip final comma
    
    % Find the institution_name field from network data
    NT_institution_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'institution_name'))));
    NT_institution_name = networkData{NT_institution_nameIndex};
    
    % Find the institution_name field from station data
    ST_DoIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_name'))));
    ST_institution_name = institutions(:,ST_DoIndex);
    ST_institution_name(cellfun('isempty',ST_institution_name)) = [];
    
    % Build the cumulative institution name list
    institutionList = [NT_institution_name; ST_institution_name];
    institution_names = institutionList(ia);
    institution_nameStr = strjoin(institution_names,'; ');
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Retrieve DoA, calibration types, calibration links and calibrations dates
try
    % Find the last_calibration_date field from station data
    ST_last_calibration_dateIndex = find(not(cellfun('isempty', strfind(stationFields, 'last_calibration_date'))));
    ST_last_calibration_date = datenum(stationData(:,ST_last_calibration_dateIndex));
    ST_last_calibration_date = ST_last_calibration_date(ST_last_calibration_date~=0);
    lastPatternStr = [sitesCodes(1,:) ': '];
    for lcd_idx=2:length(ST_last_calibration_date)
        lastPatternStr = [lastPatternStr datestr(ST_last_calibration_date(lcd_idx-1), 'yyyy-mm-dd') 'T' datestr(ST_last_calibration_date(lcd_idx-1), 'HH:MM:SS') 'Z; ' sitesCodes(lcd_idx,:) ': '];
    end
    lastPatternStr = [lastPatternStr datestr(ST_last_calibration_date(lcd_idx), 'yyyy-mm-dd') 'T' datestr(ST_last_calibration_date(lcd_idx), 'HH:MM:SS') 'Z'];
    
    % Find the DoA from station data
    ST_DoAIndex = find(not(cellfun('isempty', strfind(stationFields, 'DoA_estimation_method'))));
    ST_DoA = stationData(:,ST_DoAIndex);
    ST_DoA(cellfun('isempty',ST_DoA)) = [];
    %     ST_DoA = uniqueStrCell(ST_DoA);
    %     DoAStr = strjoin(ST_DoA,', ');
    DoAStr = [sitesCodes(1,:) ': '];
    for doa_idx=2:length(ST_DoA)
        DoAStr = [DoAStr ST_DoA{doa_idx} '; ' sitesCodes(doa_idx,:) ': '];
    end
    DoAStr = [DoAStr ST_DoA{doa_idx}];
    
    % Find the calibration_type from station data
    ST_calibration_typeIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_type'))));
    ST_calibration_type = stationData(:,ST_calibration_typeIndex);
    ST_calibration_type(cellfun('isempty',ST_calibration_type)) = [];
    %     ST_calibration_type = uniqueStrCell(ST_calibration_type);
    %     calibration_typeStr = strjoin(ST_calibration_type,', ');
    calibration_typeStr = [sitesCodes(1,:) ': '];
    for ct_idx=2:length(ST_calibration_type)
        calibration_typeStr = [calibration_typeStr ST_calibration_type{ct_idx-1} '; ' sitesCodes(ct_idx,:) ': '];
    end
    calibration_typeStr = [calibration_typeStr ST_calibration_type{ct_idx}];
    
    % Find the calibration_link from station data
    ST_calibration_linkIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_link'))));
    ST_calibration_link = stationData(:,ST_calibration_linkIndex);
    ST_calibration_link(cellfun('isempty',ST_calibration_link)) = [];
    %     ST_calibration_link = uniqueStrCell(ST_calibration_link);
    %     calibration_linkStr = strjoin(ST_calibration_link,', ');
    calibration_linkStr = [sitesCodes(1,:) ': '];
    for cl_idx=2:length(ST_calibration_link)
        calibration_linkStr = [calibration_linkStr ST_calibration_link{cl_idx-1} '; ' sitesCodes(cl_idx,:) ': '];
    end
    calibration_linkStr = [calibration_linkStr ST_calibration_link{ct_idx}];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Evaluate measurement vertical max and resolution

try
    transmit_central_frequencyIndex = find(not(cellfun('isempty', strfind(stationFields, 'transmit_central_frequency'))));
    txFreq = cell2mat(stationData(:,transmit_central_frequencyIndex)).*1e6; % transmit frequency in Hertz
    vertMax = (3e8)/(8*pi*min(txFreq));
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Prepare data
% Set netcdf format
ncfmt = 'netcdf4_classic';

% Set total data on a regular grid.
try
    lonGrid = unique(mat_tot.LonLat(:,1));
    latGrid = unique(mat_tot.LonLat(:,2));
    depth = 0;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

try
    % Prepare variables
    
    mat_tot.U_grid = NaN.*ones(length(lonGrid),length(latGrid),1);
    mat_tot.V_grid = NaN.*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.U_std = NaN.*ones(length(lonGrid),length(latGrid),1);
    mat_tot.V_std = NaN.*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.covariance = NaN.*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.GDOP = NaN.*ones(length(lonGrid),length(latGrid),1);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Populate variables
try
    for i=1:length(mat_tot.LonLat(:,1))
        lonGrid_idx = find(lonGrid==mat_tot.LonLat(i,1));
        latGrid_idx = find(latGrid==mat_tot.LonLat(i,2));
        % U and V components of current velocity
        if (not(isnan(mat_tot.U(i))))
            mat_tot.U_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.U(i)*0.01;
        end
        if (not(isnan(mat_tot.V(i))))
            mat_tot.V_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.V(i)*0.01;
        end
        % U and V standard errors
        if (not(isnan(mat_tot.ErrorEstimates(1,1).Uerr(i))))
            mat_tot.U_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Uerr(i))*0.01;
        end
        if (not(isnan(mat_tot.ErrorEstimates(1,1).Verr(i))))
            mat_tot.V_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Verr(i))*0.01;
        end
        % UV covariance
        if (not(isnan(mat_tot.ErrorEstimates(1,1).UVCovariance(i))))
            mat_tot.covariance(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).UVCovariance(i)*0.0001;
        end
        % GDOP
        if (not(isnan(mat_tot.ErrorEstimates(1,1).TotalErrors(i))))
            mat_tot.GDOP(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).TotalErrors(i);
        end
    end
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Set reference time
% if (T2C_err == 0)
%     timeref = datenum(1970,1,1);
%     time_units = ['seconds since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
%     [year,mon,day,hr,minutes,sec] = datevec(timeref);
% end

try
    timeref = datenum(1950,1,1);
    time_units = ['days since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
    %     [year,mon,day,hr,minutes,sec] = datevec(timeref);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Set data creation time and date and start and stop time and date of the
% data time window.
try
    creation = datestr(mat_tot.TimeStamp);
    creationTime = [creation(length(creation)-7:length(creation)) ' UTC'];
    stopTime = creationTime;
    creation = datevec(mat_tot.TimeStamp);
    creationDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
    stopDate = creationDate;
    creation = datestr(mat_tot.TimeStamp-1/24);
    if (length(creation) == 11)
        startTime = '00:00:00 UTC';
    else
        startTime = [creation(length(creation)-7:length(creation)) ' UTC'];
    end
    creation = datevec(mat_tot.TimeStamp-1/24);
    startDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C = 1;
end

% Set ADCC compliant data creation, coverage, resolution, duration times and spatial resolution.
try
    % File creation datetime
    dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
    % Data coverage period
    temporal_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'temporal_resolution'))));
    temporal_resolution = networkData{temporal_resolutionIndex};
    coverageStart = addtodate(mat_tot.TimeStamp, -temporal_resolution/2, 'minute');
    timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
    coverageEnd = addtodate(mat_tot.TimeStamp, temporal_resolution/2, 'minute');
    timeCoverageEnd = [datestr(coverageEnd, 'yyyy-mm-dd') 'T' datestr(coverageEnd, 'HH:MM:SS') 'Z'];
    % Temporal resolution and duration
    resolutionMinutes = minutes(temporal_resolution);
    [resH,resM,resS] = hms(resolutionMinutes);
    timeCoverageResolution = 'PT';
    if(resH~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resH) 'H'];
    end
    if(resM~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resM) 'M'];
    end
    if(resS~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resS) 'S'];
    end
    timeCoverageDuration = timeCoverageResolution;
    % Geospatial resolution
    latRes = mean(diff(latGrid));
    lonRes = mean(diff(lonGrid));
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C = 1;
end

% Set nc output file name
try
    ts = datevec(mat_tot.TimeStamp);
    time_str = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
    outputPathIndex = find(not(cellfun('isempty', strfind(networkFields, 'total_HFRnetCDF_folder_path'))));
    network_idIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_id'))));
    networkData{outputPathIndex} = strtrim(networkData{outputPathIndex});
    [tFB_err, ncFilePath] = totalFolderBuilder_v212(networkData{outputPathIndex}, timestamp);
    if(tFB_err == 0)
        ncfile = [ncFilePath filesep networkData{network_idIndex} '-Total_' time_str '.nc'];
        ncFileNoPath = [networkData{network_idIndex} '-Total_' time_str '.nc'];
    else
        disp(['[' datestr(now) '] - - ERROR in building the folder structure.']);
        return
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Set citation string and distribution string
try
    citation_statementIndex = find(not(cellfun('isempty', strfind(networkFields, 'citation_statement'))));
    citation_str = ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ' networkData{citation_statementIndex}];
    distribution_str = 'These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Set naming authority
try
    institution_websiteIndex = find(not(cellfun('isempty', strfind(networkFields, 'institution_website'))));
    institution_websiteStr = networkData{institution_websiteIndex};
    if(~isempty(strfind(institution_websiteStr,'http://')))
        tmpStr = strrep(institution_websiteStr,'http://','');
    elseif(~isempty(strfind(institution_websiteStr,'https://')))
        tmpStr = strrep(institution_websiteStr,'https://','');
    else
        tmpStr = institution_websiteStr;
    end
    tmpStr = strrep(tmpStr,'www.','');
    tmpStr = strrep(tmpStr,'/','');
    splitStr = strsplit(tmpStr,'.');
    naming_authorityStr = [];
    for split_idx=length(splitStr):-1:1
        naming_authorityStr = [naming_authorityStr splitStr{split_idx}];
        if(split_idx~=1)
            naming_authorityStr = [naming_authorityStr '.'];
        end
    end
    naming_authorityStr= naming_authorityStr(~isspace(naming_authorityStr));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Set collection time
try
    ts = datevec(mat_tot.TimeStamp);
    time_coll = [datestr(ts, 'yyyy-mm-dd') 'T' datestr(ts, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

try
    % Define EDIOS codes, site code, platform code, id and metadata resources
    EDIOS_Series_ID = networkData{network_idIndex};
    site_code = EDIOS_Series_ID;
    platform_code = [EDIOS_Series_ID '-Total'];
    dataID = [EDIOS_Series_ID '-Total_' datestr(mat_tot.TimeStamp, 'yyyy-mm-dd') 'T' datestr(mat_tot.TimeStamp, 'HH:MM:SS') 'Z'];
    metadata_pageIndex = find(not(cellfun('isempty', strfind(networkFields, 'metadata_page'))));
    TDS_catalog = networkData{metadata_pageIndex};
    xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Gets dimensions
try
    time_dim = size(mat_tot.TimeStamp,1);
    %         time_dim = netcdf.getConstant('unlimited');
    lat_dim = size(latGrid,1);
    lon_dim = size(lonGrid,1);
    depth_dim = 1;
    maxSite_dim = 50;
    refMax_dim = 1;
    maxInst_dim = length(EDMO_code);
    string15_dim = 15;
    string50_dim = 50;
    string250_dim = 250;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Build structures containing QC tests parameters

try
    GDOPthreshIndex = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_GDOP_threshold'))));
    Total_QC_params.GDOPThr = networkData{GDOPthreshIndex};
    var_thr_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_variance_threshold'))));
    Total_QC_params.VarThr = networkData{var_thr_Index};
    temp_der_thr_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_temporal_derivative_threshold'))));
    Total_QC_params.TempDerThr.threshold = networkData{temp_der_thr_Index};
    maxspd_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_velocity_threshold'))));
    Total_QC_params.VelThr = networkData{maxspd_Index};
    dataDens_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_data_density_threshold'))));
    Total_QC_params.DataDensityThr = networkData{dataDens_Index};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Perform QC tests

try
    % Build the names of the files of the previous two hours
    [twoHoursBefore, oneHourBefore] = twoPastHours(mat_tot.TimeStamp,temporal_resolution);
    Total_QC_params.TempDerThr.hour2 = [ncFilePath(1:length(ncFilePath)-length(twoHoursBefore.fP)) twoHoursBefore.fP filesep networkData{network_idIndex} '-Total_' twoHoursBefore.TS '.nc'];
    Total_QC_params.TempDerThr.hour1 = [ncFilePath(1:length(ncFilePath)-length(oneHourBefore.fP)) oneHourBefore.fP filesep networkData{network_idIndex} '-Total_' oneHourBefore.TS '.nc'];
    
    [overall_QCflag, varianceThreshold_QCflag, temporalDerivativeThreshold_QCflag, GDOP_QCflag, dataDensity_QCflag, velocityThreshold_QCflag] = TotalQCtests_v11(mat_tot, Total_QC_params);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end


%%

%% Set the time, position and depth quality flags
try
    % Time quality flag
    sdnTime_QCflag = 1;
    
    % Position quality flag
    sdnPosition_QCflag = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(length(lonGrid),length(latGrid),1));
    sdnPosition_QCflag(~isnan(mat_tot.U_grid)) = 1;
    
    % Depth quality flag
    sdnDepth_QCflag = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

%% Creates vars with their dimensions
try
    nccreate(ncfile,'TIME',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','double',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LATITUDE',...
        'Dimensions',{'LATITUDE',lat_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LONGITUDE',...
        'Dimensions',{'LONGITUDE',lon_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'crs',...
        'Datatype','int16',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_CRUISE',...
        'Dimensions',{'STRING50',string50_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_STATION',...
        'Dimensions',{'STRING50',string50_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_LOCAL_CDI_ID',...
        'Dimensions',{'STRING50',string50_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_EDMO_CODE',...
        'Dimensions',{'MAXINST',maxInst_dim,'TIME',time_dim},...
        'Datatype','int16',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_REFERENCES',...
        'Dimensions',{'STRING200', string250_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_XLINK',...
        'Dimensions',{'STRING200',string250_dim, 'REFMAX',refMax_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPH',...
        'Dimensions',{'DEPTH',depth_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'EWCT',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NSCT',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'EWCS',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NSCS',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CCOV',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'GDOP',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'TIME_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'POSITION_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPH_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'QCflag',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'VART_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'GDOP_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DDNS_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CSPD_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NARX',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NATX',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLTR',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLNR',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLTT',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLNT',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDR',...
        'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDT',...
        'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
        'Format',ncfmt);
    
    %% Creates attributes for the variables
    ncwriteatt(ncfile,'TIME','long_name',char('Time'));
    ncwriteatt(ncfile,'TIME','standard_name',char('time'));
    ncwriteatt(ncfile,'TIME','units',char(time_units));
    ncwriteatt(ncfile,'TIME','calendar',char('Gregorian'));
    ncwriteatt(ncfile,'TIME','axis',char('T'));
    ncwriteatt(ncfile,'TIME','sdn_parameter_name',char('Elapsed time (since 1950-01-01T00:00:00Z)'));
    ncwriteatt(ncfile,'TIME','sdn_parameter_urn',char('SDN:P01::ELTJLD01'));
    ncwriteatt(ncfile,'TIME','sdn_uom_name',char('Days'));
    ncwriteatt(ncfile,'TIME','sdn_uom_urn',char('SDN:P06::UTAA'));
    ncwriteatt(ncfile,'TIME','ancillary_variables',char('TIME_QC'));
    
    ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude of each location'));
    ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
    ncwriteatt(ncfile,'LATITUDE','units',char('degree_north'));
    ncwriteatt(ncfile,'LATITUDE','axis',char('Y'));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_urn',char('SDN:P06::DEGN'));
    ncwriteatt(ncfile,'LATITUDE','grid_mapping',char('crs'));
    ncwriteatt(ncfile,'LATITUDE','ancillary_variables',char('POSITION_QC'));
    
    ncwriteatt(ncfile,'LONGITUDE','long_name',char('Longitude of each location'));
    ncwriteatt(ncfile,'LONGITUDE','standard_name',char('longitude'));
    ncwriteatt(ncfile,'LONGITUDE','units',char('degree_east'));
    ncwriteatt(ncfile,'LONGITUDE','axis',char('X'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_uom_urn',char('SDN:P06::DEGE'));
    ncwriteatt(ncfile,'LONGITUDE','grid_mapping',char('crs'));
    ncwriteatt(ncfile,'LONGITUDE','ancillary_variables',char('POSITION_QC'));
    
    ncwriteatt(ncfile,'crs','grid_mapping_name',char('latitude_longitude'));
    ncwriteatt(ncfile,'crs','epsg_code',char('EPSG:4326'));
    ncwriteatt(ncfile,'crs','semi_major_axis',6378137.0);
    ncwriteatt(ncfile,'crs','inverse_flattening',298.257223563);
    
    ncwriteatt(ncfile,'SDN_CRUISE','long_name',char('Grid grouping label'));
    
    ncwriteatt(ncfile,'SDN_STATION','long_name',char('Grid label'));
    
    ncwriteatt(ncfile,'SDN_LOCAL_CDI_ID','long_name',char('SeaDataCloud CDI identifier'));
    ncwriteatt(ncfile,'SDN_LOCAL_CDI_ID','cf_role',char('grid_id'));
    
    ncwriteatt(ncfile,'SDN_EDMO_CODE','long_name',char('European Directory of Marine Organisations code for the CDI partner'));
    ncwriteatt(ncfile,'SDN_EDMO_CODE','units',char('1'));
    
    ncwriteatt(ncfile,'SDN_REFERENCES','long_name',char('Usage metadata reference'));
    
    ncwriteatt(ncfile,'SDN_XLINK','long_name',char('External resource linkages'));
    
    ncwriteatt(ncfile,'DEPH','long_name',char('Depth of measurement'));
    ncwriteatt(ncfile,'DEPH','standard_name',char('depth'));
    ncwriteatt(ncfile,'DEPH','units',char('m'));
    ncwriteatt(ncfile,'DEPH','axis',char('Z'));
    ncwriteatt(ncfile,'DEPH','positive',char('down'));
    ncwriteatt(ncfile,'DEPH','reference',char('sea_level'));
    ncwriteatt(ncfile,'DEPH','sdn_parameter_name',char('Depth below surface of the water body'));
    ncwriteatt(ncfile,'DEPH','sdn_parameter_urn',char('SDN:P01::ADEPZZ01'));
    ncwriteatt(ncfile,'DEPH','sdn_uom_name',char('Metres'));
    ncwriteatt(ncfile,'DEPH','sdn_uom_urn',char('SDN:P06::ULAA'));
    ncwriteatt(ncfile,'DEPH','ancillary_variables',char('DEPH_QC'));
    
    ncwriteatt(ncfile,'EWCT','long_name',char('West-east current component'));
    ncwriteatt(ncfile,'EWCT','standard_name',char('eastward_sea_water_velocity'));
    ncwriteatt(ncfile,'EWCT','units',char('m s-1'));
    ncwriteatt(ncfile,'EWCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EWCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EWCT','ioos_category',char('Currents'));
    ncwriteatt(ncfile,'EWCT','coordsys',char('geographic'));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_name',char('Eastward current velocity in the water body'));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_urn',char('SDN:P01::LCEWZZ01'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EWCT','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    %        ncwriteatt(ncfile,'EWCT','cell_methods',char('time: mean over hours time'));
%     ncwriteatt(ncfile,'EWCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'EWCT','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCT','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
    
    ncwriteatt(ncfile,'NSCT','long_name',char('South-north current component'));
    ncwriteatt(ncfile,'NSCT','standard_name',char('northward_sea_water_velocity'));
    ncwriteatt(ncfile,'NSCT','units',char('m s-1'));
    ncwriteatt(ncfile,'NSCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'NSCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'NSCT','ioos_category',char('Currents'));
    ncwriteatt(ncfile,'NSCT','coordsys',char('geographic'));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_name',char('Northward current velocity in the water body'));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_urn',char('SDN:P01::LCNSZZ01'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'NSCT','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    %        ncwriteatt(ncfile,'NSCT','cell_methods',char('time: mean over hours time'));
%     ncwriteatt(ncfile,'NSCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'NSCT','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCT','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
    
    ncwriteatt(ncfile,'EWCS','long_name',char('Standard Deviation of Surface Eastward Sea Water Velocity'));
    %        ncwriteatt(ncfile,'EWCS','standard_name',char('surface_eastward_sea_water_velocity_standard_error'));
    ncwriteatt(ncfile,'EWCS','units',char('m s-1'));
%     ncwriteatt(ncfile,'EWCS','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'EWCS','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'EWCS','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCS','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCS','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EWCS','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EWCS','sdn_parameter_name',char('Eastward current velocity standard deviation in the water body'));
    ncwriteatt(ncfile,'EWCS','sdn_parameter_urn',char('SDN:P01::SDEWZZZZ'));
    ncwriteatt(ncfile,'EWCS','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EWCS','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EWCS','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'NSCS','long_name',char('Standard Deviation of Surface Northward Sea Water Velocity'));
    %        ncwriteatt(ncfile,'NSCS','standard_name',char('surface_northward_sea_water_velocity_standard_error'));
    ncwriteatt(ncfile,'NSCS','units',char('m s-1'));
%     ncwriteatt(ncfile,'NSCS','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'NSCS','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'NSCS','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCS','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCS','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'NSCS','add_offset',double(addOffset));
    ncwriteatt(ncfile,'NSCS','sdn_parameter_name',char('Northward current velocity standard deviation in the water body'));
    ncwriteatt(ncfile,'NSCS','sdn_parameter_urn',char('SDN:P01::SDNSZZZZ'));
    ncwriteatt(ncfile,'NSCS','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'NSCS','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'NSCS','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'CCOV','long_name',char('Covariance of Surface Sea Water Velocity'));
    %         ncwriteatt(ncfile,'CCOV','standard_name',char('surface_sea_water_velocity_covariance'));
    ncwriteatt(ncfile,'CCOV','units',char('m2 s-2'));
%     ncwriteatt(ncfile,'CCOV','valid_range',int32([(-10-addOffset)./(scaleFactor^2), (10-addOffset)./(scaleFactor^2)]));
    ncwriteatt(ncfile,'CCOV','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'CCOV','valid_min',int32((-10-addOffset)./(scaleFactor^2)));
    ncwriteatt(ncfile,'CCOV','valid_max',int32((10-addOffset)./(scaleFactor^2)));
    ncwriteatt(ncfile,'CCOV','scale_factor',double(scaleFactor^2));
    ncwriteatt(ncfile,'CCOV','add_offset',double(addOffset));
    ncwriteatt(ncfile,'CCOV','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'CCOV','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'CCOV','sdn_uom_name',char('Square metres per second squared'));
    ncwriteatt(ncfile,'CCOV','sdn_uom_urn',char('SDN:P06::SQM2'));
    ncwriteatt(ncfile,'CCOV','ancillary_variables',char('QCflag'));
    
    ncwriteatt(ncfile,'GDOP','long_name',char('Geometrical Dilution Of Precision'));
    %         ncwriteatt(ncfile,'GDOP','standard_name',char('gdop'));
    ncwriteatt(ncfile,'GDOP','units',char('1'));
%     ncwriteatt(ncfile,'GDOP','valid_range',int16([(-20-addOffset)./scaleFactor, (20-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'GDOP','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'GDOP','valid_min',int16((-20-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'GDOP','valid_max',int16((20-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'GDOP','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'GDOP','add_offset',double(addOffset));
    ncwriteatt(ncfile,'GDOP','comment',char(['The Geometric Dilution of Precision (GDOP) is the coefficient of the uncertainty, which relates the uncertainties in radial and velocity vectors.' ...
        ' The GDOP is a unit-less coefficient, which characterizes the effect that radar station geometry has on the measurement and position determination errors.' ...
        ' A low GDOP corresponds to an optimal geometric configuration of radar stations, and results in accurate surface current data. Essentially, GDOP is a quantitative way to relate the radial and velocity vector uncertainties.'...
        ' Setting a threshold on GDOP for total combination avoids the combination of radials with an intersection angle below a certain value.' ...
        ' GDOP is a useful metric for filtering errant velocities due to poor geometry.']));
    ncwriteatt(ncfile,'GDOP','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'GDOP','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'GDOP','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'GDOP','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'GDOP','ancillary_variables',char('QCflag, GDOP_QC'));
    
    ncwriteatt(ncfile,'TIME_QC','long_name',char('Time Quality Flag'));
    ncwriteatt(ncfile,'TIME_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'TIME_QC','units',char('1'));
%     ncwriteatt(ncfile,'TIME_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'TIME_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'TIME_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'TIME_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'TIME_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'TIME_QC','comment',char('OceanSITES quality flagging for temporal coordinate.'));
    ncwriteatt(ncfile,'TIME_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'TIME_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'POSITION_QC','long_name',char('Position Quality Flags'));
    ncwriteatt(ncfile,'POSITION_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'POSITION_QC','units',char('1'));
%     ncwriteatt(ncfile,'POSITION_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'POSITION_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'POSITION_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'POSITION_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'POSITION_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'POSITION_QC','comment',char('OceanSITES quality flagging for position coordinates.'));
    ncwriteatt(ncfile,'POSITION_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'POSITION_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'DEPH_QC','long_name',char('Depth Quality Flag'));
    ncwriteatt(ncfile,'DEPH_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'DEPH_QC','units',char('1'));
%     ncwriteatt(ncfile,'DEPH_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'DEPH_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'DEPH_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'DEPH_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'DEPH_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'DEPH_QC','comment',char('OceanSITES quality flagging for depth coordinate.'));
    ncwriteatt(ncfile,'DEPH_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'DEPH_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'QCflag','long_name',char('Overall Quality Flags'));
    ncwriteatt(ncfile,'QCflag','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'QCflag','units',char('1'));
%     ncwriteatt(ncfile,'QCflag','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'QCflag','valid_min',int8(0));
    ncwriteatt(ncfile,'QCflag','valid_max',int8(9));
    ncwriteatt(ncfile,'QCflag','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'QCflag','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'QCflag','comment',char('OceanSITES quality flagging for all QC tests.'));
    ncwriteatt(ncfile,'QCflag','scale_factor',int8(1));
    ncwriteatt(ncfile,'QCflag','add_offset',int8(0));
    
    ncwriteatt(ncfile,'VART_QC','long_name',char('Variance Threshold Quality Flags'));
    ncwriteatt(ncfile,'VART_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'VART_QC','units',char('1'));
%     ncwriteatt(ncfile,'VART_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'VART_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'VART_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'VART_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'VART_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'VART_QC','comment',char(['OceanSITES quality flagging for variance threshold QC test. ' ...
        'Test not applicable to Direction Finding systems. The Temporal Derivative test is applied.' ...
        'Threshold set to ' num2str(Total_QC_params.TempDerThr.threshold) ' m/s. ']));
    ncwriteatt(ncfile,'VART_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'VART_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'GDOP_QC','long_name',char('GDOP Threshold Quality Flags'));
    ncwriteatt(ncfile,'GDOP_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'GDOP_QC','units',char('1'));
%     ncwriteatt(ncfile,'GDOP_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'GDOP_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'GDOP_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'GDOP_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'GDOP_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'GDOP_QC','comment',char(['OceanSITES quality flagging for GDOP threshold QC test. ' ...
        'Threshold set to ' num2str(Total_QC_params.GDOPThr) '.']));
    ncwriteatt(ncfile,'GDOP_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'GDOP_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'DDNS_QC','long_name',char('Data Density Threshold Quality Flags'));
    ncwriteatt(ncfile,'DDNS_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'DDNS_QC','units',char('1'));
%     ncwriteatt(ncfile,'DDNS_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'DDNS_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'DDNS_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'DDNS_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'DDNS_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'DDNS_QC','comment',char(['OceanSITES quality flagging for Data density threshold QC test. ' ...
        'Threshold set to ' num2str(Total_QC_params.DataDensityThr) ' radials.']));
    ncwriteatt(ncfile,'DDNS_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'DDNS_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity Threshold Quality Flags'));
    ncwriteatt(ncfile,'CSPD_QC','conventions',char('Copernicus Marine in situ reference table 2'));
    ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
%     ncwriteatt(ncfile,'CSPD_QC','valid_range',int8([0 9]));
    ncwriteatt(ncfile,'CSPD_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'CSPD_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'CSPD_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'CSPD_QC','comment',char(['OceanSITES quality flagging for Velocity threshold QC test. ' ...
        'Threshold set to ' num2str(Total_QC_params.VelThr) ' m/s.']));
    ncwriteatt(ncfile,'CSPD_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'CSPD_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'NARX','long_name',char('Number of Receive Antennas'));
    ncwriteatt(ncfile,'NARX','units',char('1'));
%     ncwriteatt(ncfile,'NARX','valid_range',int8([0 maxSite_dim]));
    ncwriteatt(ncfile,'NARX','valid_min',int8(0));
    ncwriteatt(ncfile,'NARX','valid_max',int8(maxSite_dim));
    ncwriteatt(ncfile,'NARX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NARX','add_offset',int8(0));
    ncwriteatt(ncfile,'NARX','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'NARX','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'NARX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NARX','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'NATX','long_name',char('Number of Transmit Antennas'));
    ncwriteatt(ncfile,'NATX','units',char('1'));
%     ncwriteatt(ncfile,'NATX','valid_range',int8([0 maxSite_dim]));
    ncwriteatt(ncfile,'NATX','valid_min',int8(0));
    ncwriteatt(ncfile,'NATX','valid_max',int8(maxSite_dim));
    ncwriteatt(ncfile,'NATX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NATX','add_offset',int8(0));
    ncwriteatt(ncfile,'NATX','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'NATX','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'NATX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NATX','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'SLTR','long_name',char('Receive Antenna Latitudes'));
    ncwriteatt(ncfile,'SLTR','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTR','units','degree_north');
%     ncwriteatt(ncfile,'SLTR','valid_range',int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLTR','valid_min',int32((-90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTR','valid_max',int32((90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_urn',char('SDN:P06::DEGN'));
    
    ncwriteatt(ncfile,'SLNR','long_name',char('Receive Antenna Longitudes'));
    ncwriteatt(ncfile,'SLNR','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNR','units','degree_east');
%     ncwriteatt(ncfile,'SLNR','valid_range',int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLNR','valid_min',int32((-180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNR','valid_max',int32((180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_urn',char('SDN:P06::DEGE'));
    
    ncwriteatt(ncfile,'SLTT','long_name',char('Transmit Antenna Latitudes'));
    ncwriteatt(ncfile,'SLTT','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTT','units','degree_north');
%     ncwriteatt(ncfile,'SLTT','valid_range',int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLTT','valid_min',int32((-90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTT','valid_max',int32((90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_urn',char('SDN:P06::DEGN'));
    
    ncwriteatt(ncfile,'SLNT','long_name',char('Transmit Antenna Longitudes'));
    ncwriteatt(ncfile,'SLNT','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNT','units','degree_east');
%     ncwriteatt(ncfile,'SLNT','valid_range',int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLNT','valid_min',int32((-180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNT','valid_max',int32((180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_urn',char('SDN:P06::DEGE'));
    
    ncwriteatt(ncfile,'SCDR','long_name',char('Receive Antenna Codes'));
    ncwriteatt(ncfile,'SCDR','units',char('1'));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDR','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'SCDT','long_name',char('Transmit Antenna Codes'));
    ncwriteatt(ncfile,'SCDT','units',char('1'));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDT','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    %% Writes values in variables
    %         ncwrite(ncfile,'TIME',int32((mat_tot.TimeStamp-timeref)*86400));
    ncwrite(ncfile,'TIME',mat_tot.TimeStamp-timeref);
    ncwrite(ncfile,'LATITUDE',latGrid);
    ncwrite(ncfile,'LONGITUDE',lonGrid);
    ncwrite(ncfile,'crs',0);
    ncwrite(ncfile,'SDN_CRUISE',site_code');
    ncwrite(ncfile,'SDN_STATION',platform_code');
    ncwrite(ncfile,'SDN_LOCAL_CDI_ID',dataID');
    ncwrite(ncfile,'SDN_EDMO_CODE',EDMO_code');
    ncwrite(ncfile,'SDN_REFERENCES',TDS_catalog');
    ncwrite(ncfile,'SDN_XLINK',xlink');
    ncwrite(ncfile,'DEPH',depth);
    ncwrite(ncfile,'EWCT',mat_tot.U_grid);
    ncwrite(ncfile,'NSCT',mat_tot.V_grid);
    ncwrite(ncfile,'EWCS',mat_tot.U_std);
    ncwrite(ncfile,'NSCS',mat_tot.V_std);
    ncwrite(ncfile,'CCOV',mat_tot.covariance);
    ncwrite(ncfile,'GDOP',mat_tot.GDOP);
    ncwrite(ncfile,'NARX',length(sitesLat));
    ncwrite(ncfile,'NATX',length(sitesLat));
    ncwrite(ncfile,'SLTR',sitesLat);
    ncwrite(ncfile,'SLNR',sitesLon);
    ncwrite(ncfile,'SLTT',sitesLat);
    ncwrite(ncfile,'SLNT',sitesLon);
    ncwrite(ncfile,'SCDR',sitesCodes');
    ncwrite(ncfile,'SCDT',sitesCodes');
    ncwrite(ncfile,'TIME_QC',int8(sdnTime_QCflag));
    ncwrite(ncfile,'POSITION_QC',sdnPosition_QCflag);
    ncwrite(ncfile,'DEPH_QC',int8(sdnDepth_QCflag));
    ncwrite(ncfile,'QCflag',overall_QCflag);
    ncwrite(ncfile,'VART_QC',temporalDerivativeThreshold_QCflag);
    ncwrite(ncfile,'GDOP_QC',GDOP_QCflag);
    ncwrite(ncfile,'DDNS_QC',dataDensity_QCflag);
    ncwrite(ncfile,'CSPD_QC',velocityThreshold_QCflag);
    
    %% Define global attributes
    
    % MANDATORY ATTRIBUTES
    % Discovery and Identification
    ncwriteatt(ncfile,'/','site_code',char(site_code));
    ncwriteatt(ncfile,'/','platform_code',char(platform_code));
    ncwriteatt(ncfile,'/','data_mode',char('R'));
    ncwriteatt(ncfile,'/','DoA_estimation_method',char(DoAStr));
    ncwriteatt(ncfile,'/','calibration_type',char(calibration_typeStr));
    ncwriteatt(ncfile,'/','last_calibration_date',char(lastPatternStr));
    ncwriteatt(ncfile,'/','calibration_link',char(calibration_linkStr));
    titleIndex = find(not(cellfun('isempty', strfind(networkFields, 'title'))));
    ncwriteatt(ncfile,'/','title',char(networkData{titleIndex}));
    summaryIndex = find(not(cellfun('isempty', strfind(networkFields, 'summary'))));
    ncwriteatt(ncfile,'/','summary',char(networkData{summaryIndex}));
    ncwriteatt(ncfile,'/','source',char('coastal structure'));
    ncwriteatt(ncfile,'/','source_platform_category_code',char('17'));
    ncwriteatt(ncfile,'/','institution',char(institution_nameStr));
    ncwriteatt(ncfile,'/','institution_edmo_code',char(EDMO_codeStr));
    ncwriteatt(ncfile,'/','data_assembly_center',char('European HFR Node'));
    ncwriteatt(ncfile,'/','id',char(dataID));
    % Geo-spatial-temporal
    ncwriteatt(ncfile,'/','data_type', char('HF radar total data'));
    ncwriteatt(ncfile,'/','feature_type',char('surface'));
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    ncwriteatt(ncfile,'/','geospatial_lat_min',char(num2str(networkData{geospatial_lat_minIndex})));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    ncwriteatt(ncfile,'/','geospatial_lat_max',char(num2str(networkData{geospatial_lat_maxIndex})));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    ncwriteatt(ncfile,'/','geospatial_lon_min',char(num2str(networkData{geospatial_lon_minIndex})));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    ncwriteatt(ncfile,'/','geospatial_lon_max',char(num2str(networkData{geospatial_lon_maxIndex})));
    ncwriteatt(ncfile,'/','geospatial_vertical_min', char('0'));
    ncwriteatt(ncfile,'/','geospatial_vertical_max', char(num2str(vertMax)));
    ncwriteatt(ncfile, '/','time_coverage_start',char(timeCoverageStart));
    ncwriteatt(ncfile, '/','time_coverage_end',char(timeCoverageEnd));
    % Conventions used
    ncwriteatt(ncfile,'/','format_version',char('v2.1.2'));
    ncwriteatt(ncfile,'/','Conventions',char('CF-1.6 Jerico-Next-deliverable-D5.14 Copernicus-InSituTAC-Manual-1.0 Copernicus-InSituTAC-SRD-1.4 Copernicus-InSituTAC-ParametersList-3.1.0'));
    % Publication information
    ncwriteatt(ncfile,'/','update_interval',char('void'));
    ncwriteatt(ncfile,'/','citation',char(citation_str));
    ncwriteatt(ncfile,'/','distribution_statement',char(distribution_str));
    ncwriteatt(ncfile,'/','publisher_name',char('European HFR Node'));
    ncwriteatt(ncfile,'/','publisher_url',char('http://eurogoos.eu/'));
    ncwriteatt(ncfile,'/','publisher_email',char('euhfrnode@azti.es'));
    licenseIndex = find(not(cellfun('isempty', strfind(networkFields, 'license'))));
    ncwriteatt(ncfile,'/','license',char(networkData{licenseIndex}));
    acknowledgmentIndex = find(not(cellfun('isempty', strfind(networkFields, 'acknowledgment'))));
    ncwriteatt(ncfile,'/','acknowledgment',char(networkData{acknowledgmentIndex}));
    % Provenance
    ncwriteatt(ncfile,'/','date_created',char(dateCreated));
    ncwriteatt(ncfile,'/','history',char([time_coll ' data collected. ' dateCreated ' netCDF file created and sent to European HFR Node']));
    ncwriteatt(ncfile,'/','date_modified',char(dateCreated));
    ncwriteatt(ncfile,'/','date_update',char(dateCreated));
    ncwriteatt(ncfile,'/','processing_level',char('3B'));
    contributor_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_name'))));
    ncwriteatt(ncfile,'/','contributor_name',char(networkData{contributor_nameIndex}));
    contributor_roleIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_role'))));
    ncwriteatt(ncfile,'/','contributor_role',char(networkData{contributor_roleIndex}));
    contributor_emailIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_email'))));
    ncwriteatt(ncfile,'/','contributor_email',char(networkData{contributor_emailIndex}));
    
    % RECOMMENDED ATTRIBUTES
    % Discovery and Identification
    projectIndex = find(not(cellfun('isempty', strfind(networkFields, 'project'))));
    ncwriteatt(ncfile,'/','project',char(networkData{projectIndex}));
    ncwriteatt(ncfile,'/','naming_authority',char(naming_authorityStr));
    ncwriteatt(ncfile,'/','keywords',char('OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF'));
    ncwriteatt(ncfile,'/','keywords_vocabulary',char('GCMD Science Keywords'));
    commentIndex = find(not(cellfun('isempty', strfind(networkFields, 'comment'))));
    ncwriteatt(ncfile,'/','comment',char(networkData{commentIndex}));
    ncwriteatt(ncfile,'/','data_language',char('eng'));
    ncwriteatt(ncfile,'/','data_character_set',char('utf8'));
    ncwriteatt(ncfile,'/','metadata_language',char('eng'));
    ncwriteatt(ncfile,'/','metadata_character_set',char('utf8'));
    ncwriteatt(ncfile,'/','topic_category',char('oceans'));
    network_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_name'))));
    ncwriteatt(ncfile,'/','network',char(networkData{network_nameIndex}));
    % Geo-spatial-temporal
    areaIndex = find(not(cellfun('isempty', strfind(networkFields, 'area'))));
    ncwriteatt(ncfile,'/','area',char(networkData{areaIndex}));
    ncwriteatt(ncfile,'/','geospatial_lat_units',char('degree_north'));
    ncwriteatt(ncfile,'/','geospatial_lon_units',char('degree_east'));
    ncwriteatt(ncfile,'/','geospatial_lat_resolution',char(num2str(latRes)));
    ncwriteatt(ncfile,'/','geospatial_lon_resolution',char(num2str(lonRes)));
    ncwriteatt(ncfile,'/','geospatial_vertical_resolution', char(num2str(vertMax)));
    ncwriteatt(ncfile,'/','geospatial_vertical_units', char('m'));
    ncwriteatt(ncfile,'/','geospatial_vertical_positive', char('down'));
    ncwriteatt(ncfile, '/','time_coverage_duration',char(timeCoverageDuration));
    ncwriteatt(ncfile, '/','time_coverage_resolution',char(timeCoverageResolution));
    ncwriteatt(ncfile,'/','reference_system',char('EPSG:4806'));
    grid_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'grid_resolution'))));
    ncwriteatt(ncfile,'/','grid_resolution',char(num2str(networkData{grid_resolutionIndex})));
    ncwriteatt(ncfile,'/','cdm_data_type',char('Grid'));
    % Conventions used
    ncwriteatt(ncfile,'/','netcdf_version',char(netcdf.inqLibVers));
    ncwriteatt(ncfile,'/','netcdf_format',char(ncfmt));
    
    % OTHER ATTRIBUTES
    ncwriteatt(ncfile,'/','metadata_contact',char('lorenzo.corgnati@sp.ismar.cnr.it'));
    ncwriteatt(ncfile,'/','metadata_date_stamp',char(dateCreated));
    ncwriteatt(ncfile,'/','standard_name_vocabulary',char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));
    ncwriteatt(ncfile,'/','sensor',char('CODAR SeaSonde'));
    ncwriteatt(ncfile,'/','institution_reference',char(institution_websiteStr));
    ncwriteatt(ncfile,'/','date_issued',char(dateCreated));
    ncwriteatt(ncfile,'/','software_name',char('HFR_Combiner'));
    ncwriteatt(ncfile,'/','software_version',char('v3.3'));
    ncwriteatt(ncfile,'/','references',char('High Frequency Radar European common data and metadata model Reference Card: all you need to know about High Frequency Radar (HFR) data harmonization at a glance. http://www.marineinsitu.eu/wp-content/uploads/2018/02/HFR_Data_Model_Reference_Card_v1.pdf'));
    
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%% Retrieve information about the nc file
try
    ncfileInfo = dir(ncfile);
    ncFilesize = ncfileInfo.bytes/1024;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

%%

if(T2C_err==0)
    disp(['[' datestr(now) '] - - ' 'tot2netCDF_v33.m successfully executed.']);
end

return

