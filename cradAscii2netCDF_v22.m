%% cradAscii2netCDF_v22.m
% This function converts hourly radial files from native WERA ASCII format
% .crad_ascii in netCDF-4 format according to the European common data and metadata
% model integrating CMEMS-INSTAC and SDC CF extension requirements.
% v2.2 complies with Copernicus-InSituTAC-FormatManual-1.4, 
% Copernicus-InSituTAC-SRD-1.41 and Copernicus-InSituTAC-ParametersList-3.2.0

% INPUT:
%         radFilename: filename of the radial file to be converted (including full path)
%         networkData: cell array containing information about the network
%                      (metadata)
%         networkFields: field names of the cell array containing
%                       information about the network.
%         stationData: cell array containing information about the station
%                      (metadata)
%         stationFields: field names of the cell array containing
%                       information about the station.
%         timestamp: timestamp of the radial file to be converted

% OUTPUT:
%         cA2C_err: error flag (0 = correct, 1 = error)
%         networkData: cell array containing information about the network.
%                       It is returned in case an update of the database is
%                       needed.
%         ncFileNoPath: filename of the converted nc file, without the full path
%         ncFilesize: size of the converted nc file.


% Author: Lorenzo Corgnati
% Date: September 28, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [cA2C_err,networkData,ncFileNoPath,ncFilesize] = cradAscii2netCDF_v22(radFilename,networkData,networkFields,stationData,stationFields,timestamp)

disp(['[' datestr(now) '] - - ' 'cradAscii2netCDF_v22.m started.']);

cA2C_err = 0;

warning('off', 'all');

%% Set scale_factor and add_offset values

scaleFactor = 0.001;
addOffset = 0;

%%

%% Retrieve the file header, the data table and the column names of the data table

try
    % Retrieve top-left point of the first gridcell, cell size and number of lon and lat gridcells
    [grb,grb,grb,topLeftLat,grb,topLeftLon,grb,cellSize]=textread(radFilename,'%s %f %s %f %s %f %s %f',1,'headerlines', 7);
    assert(topLeftLat>=-90.0 & topLeftLat<=90.0,'Grid top-left latitude out of range');
    assert(topLeftLon>=-180.0 & topLeftLon<=180.0,'Grid top-left latitude out of range');
    
    [grb,grb,grb,lonCells,latCells]=textread(radFilename,'%f %f %f %u %u',1,'headerlines', 8);
    
    % Read the data table
    [latC,lonC,kur,snV,snV2,sn,pwr] = textread(radFilename,'%f %f %u %f %f %f %f','headerlines', 9);
    ascTable=[lonC,latC,kur,snV,snV2,sn,pwr];
    tableFields={'LonC'    'LatC'    'KUR'    'SNV'    'SNS'    'SNR'    'PWR'}; % TO BE CANCELLED AFTER JAN'S CONFIRMATION
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Retrieve site code and coordinates

try
    station_idIndex = find(not(cellfun('isempty', strfind(stationFields, 'station_id'))));
    siteCode = stationData{station_idIndex};
    site_lonIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lon'))));
    siteLon = stationData{site_lonIndex};
    site_latIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lat'))));
    siteLat = stationData{site_latIndex};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%


%% Create the regular grid

try
    % Retrieve lat lon points present in the file
    unqLon = (sort(unique(lonC(lonC~=-999)))).*(180/pi);
    unqLat = (sort(unique(latC(latC~=-999)))).*(180/pi);
    unqLat = flipud(unqLat);
    
    % Retrieve the steps in longitude and latitude
    stepLon = mean(diff(unqLon));
    stepLat = mean(diff(unqLat));
    
    % Extend the sequence to create the grid
    % Start of the sequence
    seqLon = topLeftLon:stepLon:unqLon(1);
    if(~isempty(seqLon))
        if((abs(seqLon(end)-unqLon(1)))<(abs(stepLon)/2))
            seqLon = seqLon(1:end-1);
        end
    end
    seqLon = [seqLon unqLon'];
    
    seqLat = topLeftLat:stepLat:unqLat(1);
    if(~isempty(seqLat))
        if((abs(seqLat(end)-unqLat(1)))<(abs(stepLat)/2))
            seqLat = seqLat(1:end-1);
        end
    end
    seqLat = [seqLat unqLat'];
    
    % End of the sequence
    for lon_idx=length(seqLon)+1:lonCells
        seqLon = [seqLon seqLon(end)+stepLon];
    end
    
    for lat_idx=length(seqLat)+1:latCells
        seqLat = [seqLat seqLat(end)+stepLat];
    end
    
%     % Shift coordinates from top-left corner to center of the grid cells
%     for lat_idx=1:latCells
%         for lon_idx=1:lonCells
%             [gridLon(lon_idx) gridLat(lat_idx)] = km2lonlat(seqLon(lon_idx),seqLat(lat_idx),cellSize/2,-(cellSize/2));
%         end
%     end

    % Change seqLon and seqLat variable names to gridLon and gridLat
    % because they are the center of the grid cells and not the corners
    gridLon = seqLon;
    gridLat = seqLat;
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Create the TUV structure and fill it with the total data

try
    % Fill the TUV structure with the total data
    [cA2C_err,mat_rad] = cradAsciiTable2RAD(ascTable,tableFields,timestamp,gridLon,gridLat,siteLon,siteLat);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Retrieve DoA, calibration types, calibration links and calibrations dates
try
    % Find the last_calibration_date field from station data
    ST_last_calibration_dateIndex = find(not(cellfun('isempty', strfind(stationFields, 'last_calibration_date'))));
    ST_last_calibration_date = datenum(stationData(ST_last_calibration_dateIndex));
    ST_last_calibration_date = ST_last_calibration_date(ST_last_calibration_date~=0);
    lastPatternStr = [datestr(ST_last_calibration_date, 'yyyy-mm-dd') 'T' datestr(ST_last_calibration_date, 'HH:MM:SS') 'Z'];
    
    % Find the DoA from station data
    ST_DoAIndex = find(not(cellfun('isempty', strfind(stationFields, 'DoA_estimation_method'))));
    DoAStr = stationData{ST_DoAIndex};
    
    % Find the calibration_type from station data
    ST_calibration_typeIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_type'))));
    calibration_typeStr = stationData{ST_calibration_typeIndex};
    
    % Find the calibration_link from station data
    ST_calibration_linkIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_link'))));
    calibration_linkStr = stationData{ST_calibration_linkIndex};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Evaluate measurement vertical max and resolution

try
    transmit_central_frequencyIndex = find(not(cellfun('isempty', strfind(stationFields, 'transmit_central_frequency'))));
    txFreq = stationData{transmit_central_frequencyIndex}*1e6; % transmit frequency in Hertz
    vertMax = (3e8)/(8*pi*txFreq);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Prepare data
% Set netcdf format
ncfmt = 'netcdf4_classic';

try
    timeref = datenum(1950,1,1);
    time_units = ['days since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
    %     [year,mon,day,hr,minutes,sec] = datevec(timeref);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set data creation time and date and start and stop time and date of the
% data time window.
try
    creation = datestr(mat_rad.TimeStamp);
    creationTime = [creation(length(creation)-7:length(creation)) ' UTC'];
    stopTime = creationTime;
    creation = datevec(mat_rad.TimeStamp);
    creationDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
    stopDate = creationDate;
    creation = datestr(mat_rad.TimeStamp-1/24);
    if (length(creation) == 11)
        startTime = '00:00:00 UTC';
    else
        startTime = [creation(length(creation)-7:length(creation)) ' UTC'];
    end
    creation = datevec(mat_rad.TimeStamp-1/24);
    startDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set ADCC compliant data creation, coverage times and spatial resolution.
try
    % File creation datetime
    dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
    % Data coverage period
    temporal_resolutionIndex = find(not(cellfun('isempty', strfind(stationFields, 'temporal_resolution'))));
    temporal_resolution = stationData{temporal_resolutionIndex};
    coverageStart = addtodate(mat_rad.TimeStamp, -temporal_resolution/2, 'minute');
    timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
    coverageEnd = addtodate(mat_rad.TimeStamp, temporal_resolution/2, 'minute');
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
    latRes = abs(stepLat);
    lonRes = abs(stepLon);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set nc output file name
try
    ts = datevec(mat_rad.TimeStamp);
    fileTime = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
    outputPathIndex = find(not(cellfun('isempty', strfind(stationFields, 'radial_HFRnetCDF_folder_path'))));
    network_idIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_id'))));
    stationData{outputPathIndex} = strtrim(stationData{outputPathIndex});
    %v2.2
    [rFB_err, ncFilePath] = radialFolderBuilder_v22(stationData{outputPathIndex},siteCode,timestamp);
    if(rFB_err == 0)
        ncfile = [ncFilePath filesep networkData{network_idIndex} '-' siteCode '_' fileTime '.nc'];
        ncFileNoPath = [networkData{network_idIndex} '-' siteCode '_' fileTime '.nc'];
    else
        disp(['[' datestr(now) '] - - ERROR in building the folder structure.']);
        return
    end
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set citation string and distribution string
try
    citation_statementIndex = find(not(cellfun('isempty', strfind(networkFields, 'citation_statement'))));
    citation_str = ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ' networkData{citation_statementIndex}];
    distribution_str = 'These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set naming authority
try
    institution_websiteIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_website'))));
    institution_websiteStr = stationData{institution_websiteIndex};
    %     if(~isempty(strfind(institution_websiteStr,'http://')))
    %         tmpStr = strrep(institution_websiteStr,'http://','');
    %     elseif(~isempty(strfind(institution_websiteStr,'https://')))
    %         tmpStr = strrep(institution_websiteStr,'https://','');
    %     else
    %         tmpStr = institution_websiteStr;
    %     end
    %     tmpStr = strrep(tmpStr,'www.','');
    %     tmpStr = strrep(tmpStr,'/','');
    %     splitStr = strsplit(tmpStr,'.');
    %     naming_authorityStr = [];
    %     for split_idx=length(splitStr):-1:1
    %         naming_authorityStr = [naming_authorityStr splitStr{split_idx}];
    %         if(split_idx~=1)
    %             naming_authorityStr = [naming_authorityStr '.'];
    %         end
    %     end
    %     naming_authorityStr= naming_authorityStr(~isspace(naming_authorityStr));
%     naming_authorityStr = 'eu.eurogoos';
    naming_authorityStr = 'Copernicus Marine In Situ';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Set collection time
try
    ts = datevec(mat_rad.TimeStamp);
    time_coll = [datestr(ts, 'yyyy-mm-dd') 'T' datestr(ts, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

try
    % Define EDIOS codes, site code, platform code, id and metadata resources
    EDIOS_Series_ID = networkData{network_idIndex};
    EDIOS_Platform_ID = siteCode;
    EDMO_codeIndex = find(not(cellfun('isempty', strfind(stationFields, 'EDMO_code'))));
    EDMO_code = stationData{EDMO_codeIndex};
    site_code = EDIOS_Series_ID;
    platform_code = [EDIOS_Series_ID '-' EDIOS_Platform_ID];
    dataID = [EDIOS_Series_ID '-' EDIOS_Platform_ID '_' datestr(mat_rad.TimeStamp, 'yyyy-mm-dd') 'T' datestr(mat_rad.TimeStamp, 'HH:MM:SS') 'Z'];
    metadata_pageIndex = find(not(cellfun('isempty', strfind(networkFields, 'metadata_page'))));
    TDS_catalog = networkData{metadata_pageIndex};
    xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Gets dimensions
try
    time_dim = size(mat_rad.TimeStamp,1);
    %         time_dim = netcdf.getConstant('unlimited');
    lat_dim = length(gridLat);
    lon_dim = length(gridLon);
    depth_dim = 1;
    maxSite_dim = size(siteCode,1);
    maxInst_dim = length(EDMO_code);
    refMax_dim = 1;
    string15_dim = 15;
    string50_dim = 50;
    string250_dim = 250;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Build structures containing QC tests parameters

try
    maxspd_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_velocity_threshold'))));
    Radial_QC_params.VelThr = stationData{maxspd_Index};
    var_thr_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_variance_threshold'))));
    Radial_QC_params.VarThr = stationData{var_thr_Index};
    temp_der_thr_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_temporal_derivative_threshold'))));
    Radial_QC_params.TempDerThr.threshold = stationData{temp_der_thr_Index};
    med_filt_RCLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_RCLim'))));
    med_filt_AngLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_AngLim'))));
    med_filt_CurLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_CurLim'))));
    Radial_QC_params.MedFilt = [stationData{med_filt_RCLim_Index},stationData{med_filt_AngLim_Index},stationData{med_filt_CurLim_Index}];
    avg_rad_bear_min_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_average_radial_bearing_min'))));
    avg_rad_bear_max_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_average_radial_bearing_max'))));
    Radial_QC_params.AvgRadBear = [stationData{avg_rad_bear_min_Index},stationData{avg_rad_bear_max_Index}];
    rad_cnt_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_radial_count_threshold'))));
    Radial_QC_params.RadCnt = stationData{rad_cnt_Index};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Perform QC tests

try
    % Build the names of the files of the previous two hours
    [twoHoursBefore, oneHourBefore] = twoPastHours(mat_rad.TimeStamp,temporal_resolution);
    Radial_QC_params.TempDerThr.hour2 = [ncFilePath(1:length(ncFilePath)-length(twoHoursBefore.fP)) twoHoursBefore.fP filesep networkData{network_idIndex} '-' stationData{station_idIndex} '_' twoHoursBefore.TS '.nc'];
    Radial_QC_params.TempDerThr.hour1 = [ncFilePath(1:length(ncFilePath)-length(oneHourBefore.fP)) oneHourBefore.fP filesep networkData{network_idIndex} '-' stationData{station_idIndex} '_' oneHourBefore.TS '.nc'];
    
    [overall_QCflag, overWater_QCflag, varianceThreshold_QCflag, temporalDerivativeThreshold_QCflag, velocityThreshold_QCflag, medianFilter_QCflag, averageRadialBearing_QC_flag, radialVelocityMedianFiltered, radialCount_QC_flag] = cradAsciiRadialQCtests_v11(mat_rad, Radial_QC_params);
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

%% Set the time, position and depth quality flags
try
    % Time quality flag
    sdnTime_QCflag = 1;
    % Position quality flag
    sdnPosition_QCflag = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(length(gridLat),length(gridLon),1));
    sdnPosition_QCflag(~isnan(mat_rad.rdva)) = 1;
    
    % Depth quality flag
    sdnDepth_QCflag = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
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
        'Dimensions',{'STRING250', string250_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_XLINK',...
        'Dimensions',{'STRING250',string250_dim, 'REFMAX',refMax_dim, 'TIME',time_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPH',...
        'Dimensions',{'DEPTH',depth_dim},...
        'Datatype','single',...
        'FillValue',netcdf.getConstant('NC_FILL_FLOAT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'RDVA',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DRVA',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
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
    
    nccreate(ncfile,'HCSS',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'EACC',...
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
    
    nccreate(ncfile,'OWTR_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'MDFL_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'VART_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CSPD_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'AVRB_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'RDCT_QC',...
        'Dimensions',{'TIME',time_dim},...
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
        'FillValue',char(' '),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDT',...
        'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',char(' '),...
        'Format',ncfmt);
    
    %% Creates attributes for the variables
    ncwriteatt(ncfile,'TIME','long_name',char('Time'));
    ncwriteatt(ncfile,'TIME','standard_name',char('time'));
    ncwriteatt(ncfile,'TIME','units',char(time_units));
    ncwriteatt(ncfile,'TIME','calendar',char('standard'));
    ncwriteatt(ncfile,'TIME','axis',char('T'));
    ncwriteatt(ncfile,'TIME','valid_min',double(-90000));
    ncwriteatt(ncfile,'TIME','valid_max',double(90000));
    ncwriteatt(ncfile,'TIME','uncertainty',char(' '));
    ncwriteatt(ncfile,'TIME','sdn_parameter_name',char('Elapsed time (since 1950-01-01T00:00:00Z)'));
    ncwriteatt(ncfile,'TIME','sdn_parameter_urn',char('SDN:P01::ELTJLD01'));
    ncwriteatt(ncfile,'TIME','sdn_uom_name',char('Days'));
    ncwriteatt(ncfile,'TIME','sdn_uom_urn',char('SDN:P06::UTAA'));
    ncwriteatt(ncfile,'TIME','ancillary_variables',char('TIME_QC'));
    
    ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude of each location'));
    ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
    ncwriteatt(ncfile,'LATITUDE','units',char('degree_north'));
    ncwriteatt(ncfile,'LATITUDE','axis',char('Y'));
    ncwriteatt(ncfile,'LATITUDE','valid_min',single(-90));
    ncwriteatt(ncfile,'LATITUDE','valid_max',single(90));
    ncwriteatt(ncfile,'LATITUDE','uncertainty',char(' '));
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
    ncwriteatt(ncfile,'LONGITUDE','valid_min',single(-180));
    ncwriteatt(ncfile,'LONGITUDE','valid_max',single(180));
    ncwriteatt(ncfile,'LONGITUDE','uncertainty',char(' '));
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
    
    ncwriteatt(ncfile,'DEPH','long_name',char('Depth'));
    ncwriteatt(ncfile,'DEPH','standard_name',char('depth'));
    ncwriteatt(ncfile,'DEPH','units',char('m'));
    ncwriteatt(ncfile,'DEPH','axis',char('Z'));
    ncwriteatt(ncfile,'DEPH','valid_min',single(-12000));
    ncwriteatt(ncfile,'DEPH','valid_max',single(12000));
    ncwriteatt(ncfile,'DEPH','uncertainty',char(' '));
    ncwriteatt(ncfile,'DEPH','positive',char('down'));
    ncwriteatt(ncfile,'DEPH','reference',char('sea_level'));
    ncwriteatt(ncfile,'DEPH','sdn_parameter_name',char('Depth below surface of the water body'));
    ncwriteatt(ncfile,'DEPH','sdn_parameter_urn',char('SDN:P01::ADEPZZ01'));
    ncwriteatt(ncfile,'DEPH','sdn_uom_name',char('Metres'));
    ncwriteatt(ncfile,'DEPH','sdn_uom_urn',char('SDN:P06::ULAA'));
    ncwriteatt(ncfile,'DEPH','ancillary_variables',char('DEPH_QC'));
    ncwriteatt(ncfile,'DEPH','data_mode',char('R'));
    
    ncwriteatt(ncfile,'RDVA','long_name',char('Radial sea water velocity away from instrument'));
    ncwriteatt(ncfile,'RDVA','standard_name',char('radial_sea_water_velocity_away_from_instrument'));
    ncwriteatt(ncfile,'RDVA','units',char('m s-1'));
    ncwriteatt(ncfile,'RDVA','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'RDVA','add_offset',double(addOffset));
    ncwriteatt(ncfile,'RDVA','sdn_parameter_name',char('Speed of water current (Eulerian measurement) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'RDVA','sdn_parameter_urn',char('SDN:P01::LCSAWVRD'));
    ncwriteatt(ncfile,'RDVA','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'RDVA','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'RDVA','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'RDVA','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'RDVA','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'RDVA','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, RDCT_QC'));
    ncwriteatt(ncfile,'RDVA','data_mode',char('R'));
    
    ncwriteatt(ncfile,'DRVA','long_name',char('Direction of radial vector away from instrument'));
    ncwriteatt(ncfile,'DRVA','standard_name',char('direction_of_radial_vector_away_from_instrument'));
    ncwriteatt(ncfile,'DRVA','units',char('degree_true'));
    ncwriteatt(ncfile,'DRVA','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'DRVA','add_offset',double(addOffset));
    ncwriteatt(ncfile,'DRVA','sdn_parameter_name',char('Direction (towards) of water current (Eulerian measurement) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'DRVA','sdn_parameter_urn',char('SDN:P01::LCDAWVRD'));
    ncwriteatt(ncfile,'DRVA','sdn_uom_name',char('Degrees True'));
    ncwriteatt(ncfile,'DRVA','sdn_uom_urn',char('SDN:P06::UABB'));
    ncwriteatt(ncfile,'DRVA','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'DRVA','valid_min',int32(0));
    ncwriteatt(ncfile,'DRVA','valid_max',int32((360-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'DRVA','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, AVRB_QC, RDCT_QC'));
    ncwriteatt(ncfile,'DRVA','data_mode',char('R'));
    
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
    ncwriteatt(ncfile,'EWCT','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCT','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EWCT','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC'));
    ncwriteatt(ncfile,'EWCT','data_mode',char('R'));
    
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
    ncwriteatt(ncfile,'NSCT','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCT','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'NSCT','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC'));
    ncwriteatt(ncfile,'NSCT','data_mode',char('R'));
    
    ncwriteatt(ncfile,'HCSS','long_name',char('Radial variance of current velocity over coverage period'));
    ncwriteatt(ncfile,'HCSS','standard_name',char(' '));
    ncwriteatt(ncfile,'HCSS','units',char('m2 s-2'));
    ncwriteatt(ncfile,'HCSS','valid_min',int32((-10-addOffset)./(scaleFactor^2)));
    ncwriteatt(ncfile,'HCSS','valid_max',int32((10-addOffset)./(scaleFactor^2)));
    ncwriteatt(ncfile,'HCSS','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'HCSS','scale_factor',double(scaleFactor^2));
    ncwriteatt(ncfile,'HCSS','add_offset',double(addOffset));
    ncwriteatt(ncfile,'HCSS','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'HCSS','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'HCSS','sdn_uom_name',char('Square metres per second squared'));
    ncwriteatt(ncfile,'HCSS','sdn_uom_urn',char('SDN:P06::SQM2'));
    ncwriteatt(ncfile,'HCSS','ancillary_variables',char('QCflag, VART_QC'));
    ncwriteatt(ncfile,'HCSS','data_mode',char('R'));
    
    ncwriteatt(ncfile,'EACC','long_name',char('Radial accuracy of current velocity over coverage period'));
    ncwriteatt(ncfile,'EACC','standard_name',char(' '));
    ncwriteatt(ncfile,'EACC','units',char('m s-1'));
    ncwriteatt(ncfile,'EACC','valid_min',int16((-10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EACC','valid_max',int16((10-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'EACC','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'EACC','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EACC','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EACC','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'EACC','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'EACC','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EACC','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EACC','ancillary_variables',char('QCflag, VART_QC'));
    ncwriteatt(ncfile,'EACC','data_mode',char('R'));
    
    ncwriteatt(ncfile,'TIME_QC','long_name',char('Time quality flag'));
    ncwriteatt(ncfile,'TIME_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'TIME_QC','units',char('1'));
    ncwriteatt(ncfile,'TIME_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'TIME_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'TIME_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'TIME_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'TIME_QC','comment',char('OceanSITES quality flagging for temporal coordinate.'));
    ncwriteatt(ncfile,'TIME_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'TIME_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'POSITION_QC','long_name',char('Position quality flag'));
    ncwriteatt(ncfile,'POSITION_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'POSITION_QC','units',char('1'));
    ncwriteatt(ncfile,'POSITION_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'POSITION_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'POSITION_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'POSITION_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'POSITION_QC','comment',char('OceanSITES quality flagging for position coordinates.'));
    ncwriteatt(ncfile,'POSITION_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'POSITION_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'DEPH_QC','long_name',char('Depth quality flag'));
    ncwriteatt(ncfile,'DEPH_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'DEPH_QC','units',char('1'));
    ncwriteatt(ncfile,'DEPH_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'DEPH_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'DEPH_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'DEPH_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'DEPH_QC','comment',char('OceanSITES quality flagging for depth coordinate.'));
    ncwriteatt(ncfile,'DEPH_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'DEPH_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'QCflag','long_name',char('Overall quality flag'));
    ncwriteatt(ncfile,'QCflag','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'QCflag','units',char('1'));
    ncwriteatt(ncfile,'QCflag','valid_min',int8(0));
    ncwriteatt(ncfile,'QCflag','valid_max',int8(9));
    ncwriteatt(ncfile,'QCflag','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'QCflag','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'QCflag','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'QCflag','comment',char('OceanSITES quality flagging for all QC tests.'));
    ncwriteatt(ncfile,'QCflag','scale_factor',int8(1));
    ncwriteatt(ncfile,'QCflag','add_offset',int8(0));
    
    ncwriteatt(ncfile,'OWTR_QC','long_name',char('Over-water quality flag'));
    ncwriteatt(ncfile,'OWTR_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'OWTR_QC','units',char('1'));
    ncwriteatt(ncfile,'OWTR_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'OWTR_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'OWTR_QC','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'OWTR_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'OWTR_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'OWTR_QC','comment',char('OceanSITES quality flagging for Over-water QC test.'));
    ncwriteatt(ncfile,'OWTR_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'OWTR_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'MDFL_QC','long_name',char('Median filter quality flag'));
    ncwriteatt(ncfile,'MDFL_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'MDFL_QC','units',char('1'));
    ncwriteatt(ncfile,'MDFL_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'MDFL_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'MDFL_QC','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'MDFL_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'MDFL_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'MDFL_QC','comment',char(['OceanSITES quality flagging for Median Filter QC test. ' ...
        'Threshold set to ' num2str(Radial_QC_params.MedFilt(1)) ' km, ' num2str(Radial_QC_params.MedFilt(2)) ' deg, ' ...
        num2str(Radial_QC_params.MedFilt(3)) ' m/s, ']));
    ncwriteatt(ncfile,'MDFL_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'MDFL_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'VART_QC','long_name',char('Variance threshold quality flag'));
    ncwriteatt(ncfile,'VART_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'VART_QC','units',char('1'));
    ncwriteatt(ncfile,'VART_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'VART_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'VART_QC','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'VART_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'VART_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'VART_QC','comment',char(['OceanSITES quality flagging for variance threshold QC test. ' ...
        'Threshold set to ' num2str(Radial_QC_params.VarThr) ' m2/s2. ']));
    ncwriteatt(ncfile,'VART_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'VART_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity threshold quality flag'));
    ncwriteatt(ncfile,'CSPD_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
    ncwriteatt(ncfile,'CSPD_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'CSPD_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'CSPD_QC','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'CSPD_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'CSPD_QC','comment',char(['OceanSITES quality flagging for Velocity threshold QC test. ' ...
        'Threshold set to ' num2str(Radial_QC_params.VelThr) ' m/s.']));
    ncwriteatt(ncfile,'CSPD_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'CSPD_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'AVRB_QC','long_name',char('Average radial bearing quality flag'));
    ncwriteatt(ncfile,'AVRB_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'AVRB_QC','units',char('1'));
    ncwriteatt(ncfile,'AVRB_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'AVRB_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'AVRB_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'AVRB_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'AVRB_QC','comment',char(['OceanSITES quality flagging for Average Radial Bearing QC test. Thresholds set to [' num2str(Radial_QC_params.AvgRadBear(1)) '-' num2str(Radial_QC_params.AvgRadBear(2)) '] deg.']));
    ncwriteatt(ncfile,'AVRB_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'AVRB_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'RDCT_QC','long_name',char('Radial count quality flag'));
    ncwriteatt(ncfile,'RDCT_QC','conventions',char('Copernicus Marine In Situ reference table 2'));
    ncwriteatt(ncfile,'RDCT_QC','units',char('1'));
    ncwriteatt(ncfile,'RDCT_QC','valid_min',int8(0));
    ncwriteatt(ncfile,'RDCT_QC','valid_max',int8(9));
    ncwriteatt(ncfile,'RDCT_QC','flag_values',int8([0 1 2 3 4 5 6 7 8 9]));
    ncwriteatt(ncfile,'RDCT_QC','flag_meanings',char('no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value'));
    ncwriteatt(ncfile,'RDCT_QC','comment',char(['OceanSITES quality flagging for Radial Count QC test. Thresholds set to ' num2str(Radial_QC_params.RadCnt) ' vectors.']));
    ncwriteatt(ncfile,'RDCT_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'RDCT_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'NARX','long_name',char('Number of receive antennas'));
    ncwriteatt(ncfile,'NARX','standard_name',char(' '));
    ncwriteatt(ncfile,'NARX','units',char('1'));
    ncwriteatt(ncfile,'NARX','valid_min',int8(0));
    ncwriteatt(ncfile,'NARX','valid_max',int8(127));
    ncwriteatt(ncfile,'NARX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NARX','add_offset',int8(0));
    ncwriteatt(ncfile,'NARX','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'NARX','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'NARX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NARX','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'NARX','data_mode',char('R'));
    
    ncwriteatt(ncfile,'NATX','long_name',char('Number of transmit antennas'));
    ncwriteatt(ncfile,'NATX','standard_name',char(' '));
    ncwriteatt(ncfile,'NATX','units',char('1'));
    ncwriteatt(ncfile,'NATX','valid_min',int8(0));
    ncwriteatt(ncfile,'NATX','valid_max',int8(127));
    ncwriteatt(ncfile,'NATX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NATX','add_offset',int8(0));
    ncwriteatt(ncfile,'NATX','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'NATX','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'NATX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NATX','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'NATX','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SLTR','long_name',char('Receive antenna latitudes'));
    ncwriteatt(ncfile,'SLTR','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTR','units','degree_north');
    ncwriteatt(ncfile,'SLTR','valid_min',int32((-90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTR','valid_max',int32((90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_urn',char('SDN:P06::DEGN'));
    ncwriteatt(ncfile,'SLTR','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SLNR','long_name',char('Receive antenna longitudes'));
    ncwriteatt(ncfile,'SLNR','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNR','units','degree_east');
    ncwriteatt(ncfile,'SLNR','valid_min',int32((-180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNR','valid_max',int32((180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_urn',char('SDN:P06::DEGE'));
    ncwriteatt(ncfile,'SLNR','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SLTT','long_name',char('Transmit antenna latitudes'));
    ncwriteatt(ncfile,'SLTT','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTT','units','degree_north');
    ncwriteatt(ncfile,'SLTT','valid_min',int32((-90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTT','valid_max',int32((90-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLTT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_urn',char('SDN:P06::DEGN'));
    ncwriteatt(ncfile,'SLTT','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SLNT','long_name',char('Transmit antenna longitudes'));
    ncwriteatt(ncfile,'SLNT','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNT','units','degree_east');
    ncwriteatt(ncfile,'SLNT','valid_min',int32((-180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNT','valid_max',int32((180-addOffset)./scaleFactor));
    ncwriteatt(ncfile,'SLNT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_urn',char('SDN:P06::DEGE'));
    ncwriteatt(ncfile,'SLNT','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SCDR','long_name',char('Receive antenna codes'));
    ncwriteatt(ncfile,'SCDR','standard_name',char(' '));
    ncwriteatt(ncfile,'SCDR','units',char('1'));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'SCDR','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDR','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'SCDR','data_mode',char('R'));
    
    ncwriteatt(ncfile,'SCDT','long_name',char('Transmit antenna codes'));
    ncwriteatt(ncfile,'SCDT','standard_name',char(' '));
    ncwriteatt(ncfile,'SCDT','units',char('1'));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_name',char(' '));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_urn',char(' '));
    ncwriteatt(ncfile,'SCDT','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDT','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'SCDT','data_mode',char('R'));
    
    %% Writes values in variables
    %         ncwrite(ncfile,'TIME',int32((mat_tot.TimeStamp-timeref)*86400));
    ncwrite(ncfile,'TIME',mat_rad.TimeStamp-timeref);
    ncwrite(ncfile,'LATITUDE',gridLat);
    ncwrite(ncfile,'LONGITUDE',gridLon);
    ncwrite(ncfile,'crs',0);
    ncwrite(ncfile,'SDN_CRUISE',site_code');
    ncwrite(ncfile,'SDN_STATION',platform_code');
    ncwrite(ncfile,'SDN_LOCAL_CDI_ID',dataID');
    ncwrite(ncfile,'SDN_EDMO_CODE',EDMO_code');
    ncwrite(ncfile,'SDN_REFERENCES',TDS_catalog');
    ncwrite(ncfile,'SDN_XLINK',xlink');
    ncwrite(ncfile,'DEPH',mat_rad.depth);
    ncwrite(ncfile,'RDVA',mat_rad.rdva');
    ncwrite(ncfile,'DRVA',mat_rad.drva');
    ncwrite(ncfile,'EWCT',mat_rad.ewct');
    ncwrite(ncfile,'NSCT',mat_rad.nsct');    
    ncwrite(ncfile,'EACC',mat_rad.eacc');
    ncwrite(ncfile,'HCSS',mat_rad.hcss');
    ncwrite(ncfile,'NARX',length(siteLat));
    ncwrite(ncfile,'NATX',length(siteLat));
    ncwrite(ncfile,'SLTR',siteLat');
    ncwrite(ncfile,'SLNR',siteLon');
    ncwrite(ncfile,'SLTT',siteLat');
    ncwrite(ncfile,'SLNT',siteLon');
    ncwrite(ncfile,'SCDR',siteCode');
    ncwrite(ncfile,'SCDT',siteCode');
    ncwrite(ncfile,'TIME_QC',sdnTime_QCflag);
    ncwrite(ncfile,'POSITION_QC',sdnPosition_QCflag');
    ncwrite(ncfile,'DEPH_QC',sdnDepth_QCflag);
    ncwrite(ncfile,'QCflag',overall_QCflag');
    ncwrite(ncfile,'OWTR_QC',overWater_QCflag');
    ncwrite(ncfile,'VART_QC',varianceThreshold_QCflag');
    ncwrite(ncfile,'CSPD_QC',velocityThreshold_QCflag');
    ncwrite(ncfile,'MDFL_QC',medianFilter_QCflag');
    ncwrite(ncfile,'AVRB_QC',averageRadialBearing_QC_flag);
    ncwrite(ncfile,'RDCT_QC',radialCount_QC_flag);
        
    %% Define global attributes
    
    % MANDATORY ATTRIBUTES
    % Discovery and Identification
    ncwriteatt(ncfile,'/','site_code',char(site_code));
    ncwriteatt(ncfile,'/','platform_code',char(platform_code));
    ncwriteatt(ncfile,'/','platform_name',char(platform_code));
    ncwriteatt(ncfile,'/','wmo_platform_code',char(''));
    ncwriteatt(ncfile,'/','ices_platform_code',char(''));
    ncwriteatt(ncfile,'/','data_mode',char('R'));
    ncwriteatt(ncfile,'/','DoA_estimation_method',char(DoAStr));
    ncwriteatt(ncfile,'/','calibration_type',char(calibration_typeStr));
    ncwriteatt(ncfile,'/','last_calibration_date',char(lastPatternStr));
    ncwriteatt(ncfile,'/','calibration_link',char(calibration_linkStr));
    titleIndex = find(not(cellfun('isempty', strfind(networkFields, 'title'))));
    ncwriteatt(ncfile,'/','title',char(networkData{titleIndex}));
    summaryIndex = find(not(cellfun('isempty', strfind(stationFields, 'summary'))));
    ncwriteatt(ncfile,'/','summary',char(stationData{summaryIndex}));
    ncwriteatt(ncfile,'/','source',char('coastal structure'));
    ncwriteatt(ncfile,'/','source_platform_category_code',char('17'));
    institution_nameIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_name'))));
    ncwriteatt(ncfile,'/','institution',char(stationData{institution_nameIndex}));
    ncwriteatt(ncfile,'/','institution_edmo_code',char(num2str(EDMO_code)));
    ncwriteatt(ncfile,'/','data_assembly_center',char('European HFR Node'));
    ncwriteatt(ncfile,'/','id',char(dataID));
    % Geo-spatial-temporal
    ncwriteatt(ncfile,'/','data_type', char('HF radar radial data'));
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
    ncwriteatt(ncfile,'/','time_coverage_start',char(timeCoverageStart));
    ncwriteatt(ncfile,'/','time_coverage_end',char(timeCoverageEnd));
    ncwriteatt(ncfile,'/','bottom_depth', char(''));
    % Conventions used
    ncwriteatt(ncfile,'/','format_version',char('1.4'));
    ncwriteatt(ncfile,'/','Conventions',char('CF-1.6 Copernicus-InSituTAC-FormatManual-1.42 Copernicus-InSituTAC-SRD-1.5 Copernicus-InSituTAC-ParametersList-3.2.1'));
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
    ncwriteatt(ncfile,'/','processing_level',char('2B'));
    contributor_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_name'))));
    ncwriteatt(ncfile,'/','contributor_name',char(networkData{contributor_nameIndex}));
    contributor_roleIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_role'))));
    ncwriteatt(ncfile,'/','contributor_role',char(networkData{contributor_roleIndex}));
    contributor_emailIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_email'))));
    ncwriteatt(ncfile,'/','contributor_email',char(networkData{contributor_emailIndex}));
    ncwriteatt(ncfile,'/','contact',char('euhfrnode@azti.es cmems-service@ifremer.fr'));
    
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
    ncwriteatt(ncfile,'/','geospatial_lat_units',char('degrees_north'));
    ncwriteatt(ncfile,'/','geospatial_lon_units',char('degrees_east'));
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
    ncwriteatt(ncfile,'/','cdm_data_type',char('grid'));
    % Conventions used
%     ncwriteatt(ncfile,'/','netcdf_version',char(netcdf.inqLibVers));
    ncwriteatt(ncfile,'/','netcdf_format',char(ncfmt));
    ncwriteatt(ncfile,'/','netcdf_version',char('netCDF-4 classic model'));
    
    % OTHER ATTRIBUTES
    ncwriteatt(ncfile,'/','metadata_contact',char('lorenzo.corgnati@sp.ismar.cnr.it'));
    ncwriteatt(ncfile,'/','metadata_date_stamp',char(dateCreated));
    ncwriteatt(ncfile,'/','standard_name_vocabulary',char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));
    ncwriteatt(ncfile,'/','sensor',char('WERA'));
    ncwriteatt(ncfile,'/','institution_references',char(institution_websiteStr));
    ncwriteatt(ncfile,'/','date_issued',char(dateCreated));
    ncwriteatt(ncfile,'/','software_name',char('HFR_Node__Centralized_Processing'));
    ncwriteatt(ncfile,'/','software_version',char('v2.2'));
    ncwriteatt(ncfile,'/','references',char('http://marine.copernicus.eu http://www.marineinsitu.eu http://www.marineinsitu.eu/wp-content/uploads/2018/02/HFR_Data_Model_Reference_Card_v1.pdf'));
    ncwriteatt(ncfile,'/','doi',char(''));
    ncwriteatt(ncfile,'/','pi_name',char(''));
    ncwriteatt(ncfile,'/','qc_manual',char('Recommendation Report 2 on improved common procedures for HFR QC analysis http://dx.doi.org/10.25607/OBP-944'));
    ncwriteatt(ncfile,'/','wmo_inst_type',char(''));
    
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%% Retrieve information about the nc file
try
    ncfileInfo = dir(ncfile);
    ncFilesize = ncfileInfo.bytes/1024;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cA2C_err = 1;
end

%%

if(cA2C_err==0)
    disp(['[' datestr(now) '] - - ' 'cradAscii2netCDF_v22.m successfully executed.']);
end

return

