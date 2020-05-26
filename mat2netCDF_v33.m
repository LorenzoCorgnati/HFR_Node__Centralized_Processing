%% mat2netCDF_v33.m
% This function converts hourly radial files from MetNo remppaed Codar format .ruv
% in netCDF-4 format according to the European common data and metadata
% model integrating CMEMS-INSTAC and SDC CF extension requirements.
% v3.2 uses scale_factor and add_offset for packing geophysical variables.
% v3.3 names depth dimension as DEPTH; don not apply _SEADATANET_ string in
% TIME, POSITION and DEPTH QC variables; uses blank space as separator in
% institution_edmo_code global attribute.

% INPUT:
%         HFRP_RUV: structure containing the radial data
%         networkData: cell array containing information about the network
%                      (metadata)
%         networkFields: field names of the cell array containing
%                       information about the network.
%         stationData: cell array containing information about the station
%                      (metadata)
%         stationFields: field names of the cell array containing
%                       information about the station.
%         timestamp: timestamp of the total file to be converted

% OUTPUT:
%         R2C_err: error flag (0 = correct, 1 = error)
%         networkData: cell array containing information about the network.
%                       It is returned in case an update of the database is
%                       needed.
%         stationData: cell array containing information about the station.
%                       It is returned in case an update of the database is
%                       needed.
%         ncFileNoPath: filename of the converted nc file, without the full path
%         ncFilesize: size of the converted nc file.
%         station_tbUpdateFlag: flag for updating the station_tb table of the database.



% Author: Lorenzo Corgnati
% Date: November 9, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [R2C_err,networkData,stationData,ncFileNoPath,ncFilesize,station_tbUpdateFlag] = mat2netCDF_v33(HFRP_RUV,networkData,networkFields,stationData,stationFields,timestamp)

disp(['[' datestr(now) '] - - ' 'mat2netCDF_v33.m started.']);

R2C_err = 0;

PatternDate = '0000 00 00  00 00 00';

warning('off', 'all');

%% Set scale_factor and add_offset values

scaleFactor = 0.001;
addOffset = 0;

%%

%% Set the flag for updating the station_tb table of the database

station_tbUpdateFlag = 0;

%%

%% Retrieve the LLUVSpec version

try
    % Read the file header
    radHeader = HFRP_RUV.OtherMetadata.Header;
    
    % Retrieve the LLUVSpec version
    [LSC_err,LLUVSpec] = LLUVSpecChecker(radHeader);
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Set the origin coordinates

try
    site_lonIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lon'))));
    site_latIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lat'))));
    Origin = [stationData{site_latIndex},stationData{site_lonIndex}];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end


%%

%% Define the maximum number of allowed sites and of xlink references

maxsite = 1;
refmax = 1;

%%

%% Build the output file name and the citation and distribution statements

try
    switch HFRP_RUV.Type
        case 'RDLMeasured'
            fileType = 'm';
        case 'RDLIdeal'
            fileType = 'i';
    end
    
    station_idIndex = find(not(cellfun('isempty', strfind(stationFields, 'station_id'))));
    siteCode = stationData{station_idIndex};
    
    ts = datevec(HFRP_RUV.TimeStamp);
    fileTime = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
    
    outputPathIndex = find(not(cellfun('isempty', strfind(stationFields, 'radial_HFRnetCDF_folder_path'))));
    network_idIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_id'))));
    stationData{outputPathIndex} = strtrim(stationData{outputPathIndex});
    %v2.1.2
    [rFB_err, ncFilePath] = radialFolderBuilder_v212(stationData{outputPathIndex},siteCode,timestamp);
    if(rFB_err == 0)
        ncfile = [ncFilePath filesep networkData{network_idIndex} '-' siteCode '_' fileTime '.nc'];
        ncFileNoPath = [networkData{network_idIndex} '-' siteCode '_' fileTime '.nc'];
    else
        disp(['[' datestr(now) '] - - ERROR in building the folder structure.']);
        return
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

% Set citation string and distribution string
try
    citation_statementIndex = find(not(cellfun('isempty', strfind(networkFields, 'citation_statement'))));
    citation_str = ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ' networkData{citation_statementIndex}];
    distribution_str = 'These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Set naming authority

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
    naming_authorityStr = 'eu.eurogoos';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Define EDIOS and EDMO codes, site code, platform code, id and metadata resources

try
    EDIOS_Series_ID = networkData{network_idIndex};
    EDIOS_Platform_ID = siteCode;
    EDMO_codeIndex = find(not(cellfun('isempty', strfind(stationFields, 'EDMO_code'))));
    EDMO_code = stationData{EDMO_codeIndex};
    site_code = EDIOS_Series_ID;
    platform_code = [EDIOS_Series_ID '-' EDIOS_Platform_ID];
    id = [EDIOS_Series_ID '-' EDIOS_Platform_ID '_' datestr(HFRP_RUV.TimeStamp, 'yyyy-mm-dd') 'T' datestr(HFRP_RUV.TimeStamp, 'HH:MM:SS') 'Z'];
    metadata_pageIndex = find(not(cellfun('isempty', strfind(networkFields, 'metadata_page'))));
    TDS_catalog = networkData{metadata_pageIndex};
    xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Load radial data
try
    R = struct( ...
        'lond', HFRP_RUV.OtherMetadata.RawData(:, 1), ...
        'latd', HFRP_RUV.OtherMetadata.RawData(:, 2), ...
        'velu', HFRP_RUV.OtherMetadata.RawData(:, 3), ...
        'velv', HFRP_RUV.OtherMetadata.RawData(:, 4), ...
        'vflg', HFRP_RUV.OtherMetadata.RawData(:, 5), ...
        'espc', HFRP_RUV.OtherMetadata.RawData(:, 6), ...
        'etmp', HFRP_RUV.OtherMetadata.RawData(:, 7), ...
        'maxv', HFRP_RUV.OtherMetadata.RawData(:, 8), ...
        'minv', HFRP_RUV.OtherMetadata.RawData(:, 9), ...
        'ersc', HFRP_RUV.OtherMetadata.RawData(:,10), ...
        'ertc', HFRP_RUV.OtherMetadata.RawData(:,11), ...
        'xdst', HFRP_RUV.OtherMetadata.RawData(:,12), ...
        'ydst', HFRP_RUV.OtherMetadata.RawData(:,13), ...
        'rnge', HFRP_RUV.OtherMetadata.RawData(:,14), ...
        'bear', HFRP_RUV.OtherMetadata.RawData(:,15), ...
        'velo', HFRP_RUV.OtherMetadata.RawData(:,16), ...
        'head', HFRP_RUV.OtherMetadata.RawData(:,17), ...
        'sprc', HFRP_RUV.OtherMetadata.RawData(:,18) ...
        );
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
    return
end

%%

%% Retrieve information from original CODAR fields
try
    % Retrieve information from header
    for head_idx=1:length(radHeader)
        splitLine = regexp(radHeader{head_idx}, ' ', 'split');
        
        if (strcmp(splitLine{1}, '%RangeResolutionKMeters:'))
            RangeResolutionKMeters = str2double(splitLine{2});
        end
        if (strcmp(splitLine{1}, '%RangeResolutionMeters:'))
            RangeResolutionMeters = str2double(splitLine{2});
        end
        if(strcmp(splitLine{1}, '%AngularResolution:'))
            AngularResolution = str2double(splitLine{2});
        end
        if(strcmp(splitLine{1}, '%UUID:'))
            UUID = strrep(radHeader{head_idx}(length('%UUID:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%Manufacturer:'))
            manufacturer = strrep(radHeader{head_idx}(length('%Manufacturer:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RangeStart:'))
            RangeStart = strrep(radHeader{head_idx}(length('%RangeStart:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RangeEnd:'))
            RangeEnd = strrep(radHeader{head_idx}(length('%RangeEnd:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RangeResolutionKMeters:'))
            RangeResolutionKMeters_str = strrep(radHeader{head_idx}(length('%RangeResolutionKMeters:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RangeResolutionMeters:'))
            RangeResolutionMeters_str = strrep(radHeader{head_idx}(length('%RangeResolutionMeters:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%AntennaBearing:'))
            AntennaBearing = strrep(radHeader{head_idx}(length('%AntennaBearing:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%ReferenceBearing:'))
            ReferenceBearing = strrep(radHeader{head_idx}(length('%ReferenceBearing:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%AngularResolution:'))
            AngularResolution_str = strrep(radHeader{head_idx}(length('%AngularResolution:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%SpatialResolution:'))
            SpatialResolution = strrep(radHeader{head_idx}(length('%SpatialResolution:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternDate:'))
            PatternDate = strrep(radHeader{head_idx}(length('%PatternDate:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternResolution:'))
            PatternResolution = strrep(radHeader{head_idx}(length('%PatternResolution:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%TransmitCenterFreqMHz:'))
            TransmitCenterFreqMHz = strrep(radHeader{head_idx}(length('%TransmitCenterFreqMHz:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%DopplerResolutionHzPerBin:'))
            DopplerResolutionHzPerBin = strrep(radHeader{head_idx}(length('%DopplerResolutionHzPerBin:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%FirstOrderMethod:'))
            FirstOrderMethod = strrep(radHeader{head_idx}(length('%FirstOrderMethod:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%BraggSmoothingPoints:'))
            BraggSmoothingPoints = strrep(radHeader{head_idx}(length('%BraggSmoothingPoints:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%BraggHasSecondOrder:'))
            BraggHasSecondOrder = strrep(radHeader{head_idx}(length('%BraggHasSecondOrder:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RadialBraggPeakDropOff:'))
            RadialBraggPeakDropOff = strrep(radHeader{head_idx}(length('%RadialBraggPeakDropOff:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RadialBraggPeakNull:'))
            RadialBraggPeakNull = strrep(radHeader{head_idx}(length('%RadialBraggPeakNull:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RadialBraggNoiseThreshold:'))
            RadialBraggNoiseThreshold = strrep(radHeader{head_idx}(length('%RadialBraggNoiseThreshold:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternAmplitudeCorrections:'))
            PatternAmplitudeCorrections = strrep(radHeader{head_idx}(length('%PatternAmplitudeCorrections:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternPhaseCorrections:'))
            PatternPhaseCorrections = strrep(radHeader{head_idx}(length('%PatternPhaseCorrections:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternAmplitudeCalculations:'))
            PatternAmplitudeCalculations = strrep(radHeader{head_idx}(length('%PatternAmplitudeCalculations:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternPhaseCalculations:'))
            PatternPhaseCalculations = strrep(radHeader{head_idx}(length('%PatternPhaseCalculations:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RadialMusicParameters:'))
            RadialMusicParameters = strrep(radHeader{head_idx}(length('%RadialMusicParameters:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%MergedCount:'))
            MergedCount = strrep(radHeader{head_idx}(length('%MergedCount:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%RadialMinimumMergePoints:'))
            RadialMinimumMergePoints = strrep(radHeader{head_idx}(length('%RadialMinimumMergePoints:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%FirstOrderCalc:'))
            FirstOrderCalc = strrep(radHeader{head_idx}(length('%FirstOrderCalc:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%MergeMethod:'))
            MergeMethod = strrep(radHeader{head_idx}(length('%MergeMethod:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%PatternMethod:'))
            PatternMethod = strrep(radHeader{head_idx}(length('%PatternMethod:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%TransmitSweepRateHz:'))
            TransmitSweepRateHz = strrep(radHeader{head_idx}(length('%TransmitSweepRateHz:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%TransmitBandwidthKHz:'))
            TransmitBandwidthKHz = strrep(radHeader{head_idx}(length('%TransmitBandwidthKHz:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%SpectraRangeCells:'))
            SpectraRangeCells = strrep(radHeader{head_idx}(length('%SpectraRangeCells:')+2:length(radHeader{head_idx})), '"', '');
        end
        if(strcmp(splitLine{1}, '%SpectraDopplerCells:'))
            SpectraDopplerCells = strrep(radHeader{head_idx}(length('%SpectraDopplerCells:')+2:length(radHeader{head_idx})), '"', '');
        end
    end
    
    % Evaluate last pattern date
    patternTS = datenum([str2double(PatternDate(1:5)) str2double(PatternDate(6:8)) str2double(PatternDate(9:11))]);
    if(patternTS~=0)
        lastPatternVec = datevec(patternTS);
        lastPatternStr = [datestr(lastPatternVec, 'yyyy-mm-dd') 'T' datestr(lastPatternVec, 'HH:MM:SS') 'Z'];
    else
        lastPatternStr = '';
    end
    
    % Check if the last calibration date stored in the database is updated
    last_calibration_dateIndex = find(not(cellfun('isempty', strfind(stationFields, 'last_calibration_date'))));
    if(patternTS~=0)
        if(size(stationData{last_calibration_dateIndex},2)==10)
            if(datenum(stationData{last_calibration_dateIndex})<patternTS)
                station_tbUpdateFlag = 1;
                stationData{last_calibration_dateIndex} = datestr(patternTS, 'yyyy-mm-dd');
            end
        else
            station_tbUpdateFlag = 1;
            stationData{last_calibration_dateIndex} = datestr(patternTS, 'yyyy-mm-dd');
        end
    end
    
    % Time Stamp
    R.time = HFRP_RUV.TimeStamp;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
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

%% Retrieve time coverage period, resolution, duration and metadata date stamp
try
    temporal_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'temporal_resolution'))));
    temporal_resolution = networkData{temporal_resolutionIndex};
    coverageStart = addtodate(HFRP_RUV.TimeStamp, -temporal_resolution/2, 'minute');
    timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
    coverageEnd = addtodate(HFRP_RUV.TimeStamp, temporal_resolution/2, 'minute');
    timeCoverageEnd = [datestr(coverageEnd, 'yyyy-mm-dd') 'T' datestr(coverageEnd, 'HH:MM:SS') 'Z'];
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
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Retrieve the file creation datetime (i.e. metadata date stamp)
try
    dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Set reference time
try
    timeref = datenum(1950,1,1);
    time_units = [datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%


%% Set collection time
try
    time_coll = [datestr(HFRP_RUV.TimeStamp, 'yyyy-mm-dd') 'T' datestr(HFRP_RUV.TimeStamp, 'HH:MM:SS') 'Z'];
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Generate grid domain
try
    
    %     range_dim = min(R.rnge):RangeResolutionKMeters:max(R.rnge);
    
    % Force the range variable to have the same number of values, in order to have
    % all radial files with the same range dimension. It's crucial for the
    % THREDDS time aggregation.
    number_of_range_cellsIndex = find(not(cellfun('isempty', strfind(stationFields, 'number_of_range_cells'))));
    if(exist('RangeResolutionKMeters','var') ~= 0)
        range_dim = 0:RangeResolutionKMeters:stationData{number_of_range_cellsIndex};
    elseif(exist('RangeResolutionMeters','var') ~= 0)
        range_dim = 0:RangeResolutionMeters*0.001:stationData{number_of_range_cellsIndex};
    end
    
    if(~exist('AngularResolution','var') ~= 0)
        AngularResolution = 5;
    end
    bearing_dim_1 = min(R.bear):-AngularResolution:0;
    bearing_dim_1 = sort(bearing_dim_1,'asc');
    bearing_dim_2 = min(R.bear):AngularResolution:max(R.bear);
    bearing_dim_3 = max(R.bear):AngularResolution:360;
    bearing_dim = unique([bearing_dim_1 bearing_dim_2 bearing_dim_3]);
    
    [bearing, range] = meshgrid(bearing_dim, range_dim);
    
    % % CHECK - see that our grid captures the sample space right
    % figure
    % plot( range, bearing, 'go' )
    % hold on
    % plot( R.rnge, R.bear, 'b.')
    
    % Generate lat/lon from bearing and range
    % This is the most accurate method of defining lat/lon values for given
    % range & bearing. Other methods of interpolating and extrapolating using
    % griddata or otherwise from range, bearing, lon and lat are inconsistent
    % and generally not reliable from file to file. The reason we need to do
    % this at all is because there cannot be missing values in auxillary
    % coordinate variables.
    [latd, lond] = reckon(Origin(1), Origin(2), km2deg(range), bearing);
    
    % Set geographical resolutions
    
    % Latitude and longitude resolution (degrees)
    lonArr = unique(lond);
    lonArr = sort(lonArr);
    latArr = unique(latd);
    latArr = sort(latArr);
    latDiff = diff(latArr);
    lonDiff = diff(lonArr);
    latRes = abs(mean(latDiff));
    lonRes = abs(mean(lonDiff));
    
    % % CHECK
    % figure
    % plot(lond, latd, 'ko')
    % hold on
    % plot(R.lond, R.latd, 'r.')
    % daspect([1.2 1 1])
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Determine the mapping index from data to full grid
try
    range_map_idx = nan(1, numel(R.rnge));
    bearing_map_idx = nan(1, numel(R.bear));
    for I = 1 : numel(R.rnge)
        [~, range_map_idx(I)] = min(abs(range_dim - R.rnge(I)));
        [~, bearing_map_idx(I)] = min(abs(bearing_dim - R.bear(I)));
    end
    clear I;
    rb_map_idx = sub2ind(size(range), range_map_idx, bearing_map_idx);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Populate gridded matrices
try
    % Eastward Velocity
    velu = nan(size(range));
    velu(rb_map_idx) = R.velu./100;
    %     velu(rb_map_idx) = 0;
    
    % Northward Velocity
    velv = nan(size(range));
    velv(rb_map_idx) = R.velv./100;
    %    velv(rb_map_idx) = 0;
    
    % Current Speed - flip sign so positive velocity is away from radar
    % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
    % CODAR reports positive speeds as toward the radar.
    velo = nan(size(range));
    velo(rb_map_idx) = -R.velo./100;
    
    % Current Heading
    head = nan(size(range));
    head(rb_map_idx) = R.head;
    
    % Current Bearing
    bear = nan(size(range));
    bear(rb_map_idx) = R.bear;
    
    % Vector On Water Flag
    owtr = nan(size(range));
    owtr(rb_map_idx) = R.vflg;
    
    % Spatial Quality
    espc = nan(size(range));
    espc(rb_map_idx) = R.espc./100;
    espc( espc==9.99 ) = NaN; % native bad-value
    
    % Temporal Quality
    etmp = nan(size(range));
    etmp(rb_map_idx) = R.etmp./100;
    etmp( etmp==9.99 ) = NaN; % native bad-value
    
    % Velocity Maximum - flip sign so positive velocity is away from radar
    % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
    % CODAR reports positive speeds as toward the radar.
    maxv = nan(size(range));
    maxv(rb_map_idx) = -R.maxv/100;
    maxv( maxv==-9.99 ) = NaN; % native bad-value
    
    % Velocity Minimum - flip sign so positive velocity is away from radar
    % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
    % CODAR reports positive speeds as toward the radar.
    minv = nan(size(range));
    minv(rb_map_idx) = -R.minv/100;
    minv( minv==-9.99 ) = NaN; % native bad-value
    
    % Spatial Count
    ersc = nan(size(range));
    ersc(rb_map_idx) = R.ersc;
    ersc(ersc>127) = 127;
    ersc( ersc==999 ) = NaN; % native bad-value
    
    % Temporal Count
    ertc = nan(size(range));
    ertc(rb_map_idx) = R.ertc;
    ertc(ertc>127) = 127;
    ertc( ertc==999 ) = NaN; % native bad-value
    
    % X-Distance
    xdst = nan(size(range));
    xdst(rb_map_idx) = R.xdst;
    
    % Y-Distance
    ydst = nan(size(range));
    ydst(rb_map_idx) = R.ydst;
    
    % Spectra Range Cell
    sprc = nan(size(range));
    sprc(rb_map_idx) = R.sprc;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
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
    [twoHoursBefore, oneHourBefore] = twoPastHours(R.time,temporal_resolution);
    Radial_QC_params.TempDerThr.hour2 = [ncFilePath(1:length(ncFilePath)-length(twoHoursBefore.fP)) twoHoursBefore.fP filesep networkData{network_idIndex} '-' stationData{station_idIndex} '_' twoHoursBefore.TS '.nc'];
    Radial_QC_params.TempDerThr.hour1 = [ncFilePath(1:length(ncFilePath)-length(oneHourBefore.fP)) oneHourBefore.fP filesep networkData{network_idIndex} '-' stationData{station_idIndex} '_' oneHourBefore.TS '.nc'];
    
    [overall_QCflag, overWater_QCflag, varianceThreshold_QCflag, temporalDerivativeThreshold_QCflag, velocityThreshold_QCflag, medianFilter_QCflag, averageRadialBearing_QC_flag, radialVelocityMedianFiltered, radialCount_QC_flag] = ruvRadialQCtests_v11(bear, lond, latd, owtr, etmp, head, velo, Radial_QC_params);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Scaling for data packing

% Bearing and heading are only reported to 10ths and can be
% reported as short unsigned integers when scaled. However, Bearing is a
% coordinate variable and cannot be scaled by CF metadata convention.
% head = round(head.*10);

velo = round((velo-addOffset)./scaleFactor);
head = round((head-addOffset)./scaleFactor);
velu = round((velu-addOffset)./scaleFactor);
velv = round((velv-addOffset)./scaleFactor);
espc = round((espc-addOffset)./scaleFactor);
etmp = round((etmp-addOffset)./scaleFactor);
maxv = round((maxv-addOffset)./scaleFactor);
minv = round((minv-addOffset)./scaleFactor);
xdst = round((xdst-addOffset)./scaleFactor);
ydst = round((ydst-addOffset)./scaleFactor);

%%

%% Fill Values
try
    % Non-scaled variables
    ersc(isnan(ersc)) = netcdf.getConstant('NC_FILL_SHORT');
    ertc(isnan(ertc)) = netcdf.getConstant('NC_FILL_SHORT');
    sprc(isnan(sprc)) = netcdf.getConstant('NC_FILL_SHORT');
    
    % Scaled variables
    velo(isnan(velo)) = netcdf.getConstant('NC_FILL_SHORT');
    head(isnan(head)) = netcdf.getConstant('NC_FILL_INT');
    radialVelocityMedianFiltered(isnan(radialVelocityMedianFiltered)) = netcdf.getConstant('NC_FILL_SHORT');
    velu(isnan(velu)) = netcdf.getConstant('NC_FILL_SHORT');
    velv(isnan(velv)) = netcdf.getConstant('NC_FILL_SHORT');
    espc(isnan(espc)) = netcdf.getConstant('NC_FILL_SHORT');
    etmp(isnan(etmp)) = netcdf.getConstant('NC_FILL_SHORT');
    maxv(isnan(maxv)) = netcdf.getConstant('NC_FILL_SHORT');
    minv(isnan(minv)) = netcdf.getConstant('NC_FILL_SHORT');
    xdst(isnan(xdst)) = netcdf.getConstant('NC_FILL_INT');
    ydst(isnan(ydst)) = netcdf.getConstant('NC_FILL_INT');
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Set the time, position and depth quality flags

try
    % Time quality flag
    sdnTime_QCflag = 1;
    
    % Position quality flag
    sdnPosition_QCflag = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(velu,1),size(velu,2),1));
    sdnPosition_QCflag(velu~=netcdf.getConstant('NC_FILL_SHORT')) = 1;
    
    % Depth quality flag
    sdnDepth_QCflag = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Populate site latitude, site longitude and site code variables

try
    siteLat = Origin(1);
    siteLon = Origin(2);
    numSites = 1;
    siteLat((length(siteLat)+1):maxsite) = netcdf.getConstant('NC_FILL_SHORT');
    siteLon((length(siteLon)+1):maxsite) = netcdf.getConstant('NC_FILL_SHORT');
    %     siteCodeArr(1,:) = siteCode;
    siteCodeArr = siteCode;
    for sC_idx=size(siteCode,1)+1:maxsite
        siteCodeArr(sC_idx,:) = blanks(size(siteCode,2));
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%


%% Define dimensions

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

try
    % Keep the classic data model to increase compatibility since we aren't using
    % features specific to the NetCDF-4.0 model.
    mode = netcdf.getConstant('NETCDF4');
    mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
    ncid = netcdf.create(ncfile, mode);
    hist_create = [time_coll ' data collected. ' dateCreated ' netCDF file created and sent to European HFR Node'];
    
    % Add dimensions (in order of T, Z, Y, X for CF/COARDS compliance)
    
    %         dimid_t = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('unlimited'));
    dimid_t = netcdf.defDim(ncid, 'TIME', 1);
    dimid_bearing = netcdf.defDim(ncid, 'BEAR', numel( bearing_dim));
    dimid_range = netcdf.defDim(ncid, 'RNGE', numel( range_dim));
    dimid_depth = netcdf.defDim(ncid, 'DEPTH', 1);
    dimid_maxsite = netcdf.defDim(ncid, 'MAXSITE', maxsite);
    dimid_maxinst = netcdf.defDim(ncid, 'MAXINST', length(EDMO_code));
    dimid_refmax = netcdf.defDim(ncid, 'REFMAX', refmax);
    dimid_string_siteCodeArr = netcdf.defDim(ncid, ['STRING' num2str(length(siteCodeArr))], length(siteCodeArr));
    dimid_string_site_code = netcdf.defDim(ncid, ['STRING' num2str(length(site_code))], length(site_code));
    dimid_string_platform_code = netcdf.defDim(ncid, ['STRING' num2str(length(platform_code))], length(platform_code));
    dimid_string_sdn_local_cdi_id = netcdf.defDim(ncid, ['STRING' num2str(length(id))], length(id));
    dimid_string_sdn_references = netcdf.defDim(ncid, ['STRING' num2str(length(TDS_catalog))], length(TDS_catalog));
    dimid_string_sdn_xlink = netcdf.defDim(ncid, ['STRING' num2str(length(xlink))], length(xlink));
    
    %dimid_lat = netcdf.defDim(ncid, 'lat', length(LatU));
    %dimid_lon = netcdf.defDim(ncid, 'lon', length(LonU));
    
    %%
    
    %% Add coordinate variables and attributes
    
    % For time coordinate variable, when using the standard or gregorian
    % calendar, units should be in seconds and shifted forward beyond
    % 1582/10/15 to avoid Julian -> Gregorian crossover and potential
    % problems with the udunits package, if used. Year length should be
    % exactly 365.2425 days long for the Gregorian/standard and proleptic
    % Gregorian calendar.
    %
    % MATLAB's datenum uses the proleptic gregorian calendar as defined
    % by the CF standard. When using standard_name attribute for time
    % coordinate be sure to include the calendar and units attributes.
    %
    % Use int for time data type because need 8 significant figures to
    % represent number of seconds in 1 year (31556952 seconds), float
    % isn't enough (7 sig. figures). Would need to use double beyond
    % int. Besides, don't need sub-second accuracy. Unix time will
    % overflow int in 2038. New epoch will be needed then to keep int
    % representation. Unix time epoch is: 1970-01-01 00:00:00Z = 0 seconds, 86,400 sec/day
    varid_t = netcdf.defVar(ncid, 'TIME', 'double', dimid_t);
    netcdf.putAtt(ncid, varid_t, 'long_name', 'Time');
    netcdf.putAtt(ncid, varid_t, 'standard_name', 'time');
    netcdf.putAtt(ncid, varid_t, 'units', ['days since ' time_units]);
    netcdf.putAtt(ncid, varid_t, 'calendar', 'Gregorian');
    netcdf.putAtt(ncid, varid_t, 'axis', 'T');
    netcdf.putAtt(ncid, varid_t, 'sdn_parameter_name', 'Elapsed time (since 1950-01-01T00:00:00Z)');
    netcdf.putAtt(ncid, varid_t, 'sdn_parameter_urn', 'SDN:P01::ELTJLD01');
    netcdf.putAtt(ncid, varid_t, 'sdn_uom_name', 'Days');
    netcdf.putAtt(ncid, varid_t, 'sdn_uom_urn', 'SDN:P06::UTAA');
    netcdf.putAtt(ncid, varid_t, 'ancillary_variables', 'TIME_QC');
    
    % Bearing (arbitrary 'y' dimension)
    varid_bearing = netcdf.defVar(ncid, 'BEAR', 'float', dimid_bearing);
    netcdf.putAtt(ncid, varid_bearing, 'axis', 'X');
    netcdf.putAtt(ncid, varid_bearing, 'long_name', 'Bearing away from instrument');
    netcdf.putAtt(ncid, varid_bearing, 'units', 'degree_true');
    netcdf.putAtt(ncid, varid_bearing, 'sdn_parameter_name', 'Bearing');
    netcdf.putAtt(ncid, varid_bearing, 'sdn_parameter_urn', 'SDN:P01::BEARRFTR');
    netcdf.putAtt(ncid, varid_bearing, 'sdn_uom_name', 'Degrees true');
    netcdf.putAtt(ncid, varid_bearing, 'sdn_uom_urn', 'SDN:P06::UABB');
    netcdf.putAtt(ncid, varid_bearing, 'ancillary_variables', 'POSITION_QC');
    
    % Range (arbitrary 'x' dimension)
    varid_range = netcdf.defVar(ncid, 'RNGE', 'float', dimid_range);
    netcdf.putAtt(ncid, varid_range, 'axis', 'Y');
    netcdf.putAtt(ncid, varid_range, 'long_name', 'Range away from instrument');
    netcdf.putAtt(ncid, varid_range, 'units', 'km');
    netcdf.putAtt(ncid, varid_range, 'sdn_parameter_name', 'Range (from fixed reference point) by unspecified GPS system');
    netcdf.putAtt(ncid, varid_range, 'sdn_parameter_urn', 'SDN:P01::RIFNAX01');
    netcdf.putAtt(ncid, varid_range, 'sdn_uom_name', 'Kilometres');
    netcdf.putAtt(ncid, varid_range, 'sdn_uom_urn', 'SDN:P06::ULKM');
    netcdf.putAtt(ncid, varid_range, 'ancillary_variables', 'POSITION_QC');
    
    % Depth (arbitrary 'z' dimension)
    varid_depth = netcdf.defVar(ncid, 'DEPH', 'float', dimid_depth);
    netcdf.putAtt(ncid, varid_depth, 'long_name', 'Depth');
    netcdf.putAtt(ncid, varid_depth, 'standard_name', 'depth');
    netcdf.putAtt(ncid, varid_depth, 'units', 'm');
    netcdf.putAtt(ncid, varid_depth, 'axis', 'Z');
    netcdf.putAtt(ncid, varid_depth, 'positive', 'down');
    netcdf.putAtt(ncid, varid_depth, 'reference', 'sea_level');
    netcdf.putAtt(ncid, varid_depth, 'sdn_parameter_name', 'Depth below surface of the water body');
    netcdf.putAtt(ncid, varid_depth, 'sdn_parameter_urn', 'SDN:P01::ADEPZZ01');
    netcdf.putAtt(ncid, varid_depth, 'sdn_uom_name', 'Metres');
    netcdf.putAtt(ncid, varid_depth, 'sdn_uom_urn', 'SDN:P06::ULAA');
    netcdf.putAtt(ncid, varid_depth, 'ancillary_variables', 'DEPH_QC');
    
    
    %% Add auxillary coordinate variables to provide mapping from range and bearing to lat, lon.
    
    % A minimum of 4 significant figures to the right of the decimal
    % place is needed to keep resolution below 10's of meters. Using
    % float data type for lat yields at least 5 significant digits to
    % the right of the decimal giving ~1/2m resolution in latitude.
    % Nine significant figures (at least 7 to the right of the
    % decimal) could be achieved using int data type but need to
    % introduce a scale factor.
    
    % Latitude
    varid_lat = netcdf.defVar( ncid, 'LATITUDE', 'float', [dimid_range dimid_bearing] );
    netcdf.defVarDeflate(ncid, varid_lat, true, true, 6);
    netcdf.putAtt( ncid, varid_lat, 'standard_name', 'latitude' );
    netcdf.putAtt( ncid, varid_lat, 'long_name', 'Latitude of each location' );
    netcdf.putAtt( ncid, varid_lat, 'units', 'degree_north' );
    %     netcdf.putAtt(ncid, varid_lat, 'valid_range', single( [-90 90]));
    netcdf.putAtt(ncid, varid_lat, 'valid_min', single(-90));
    netcdf.putAtt(ncid, varid_lat, 'valid_max', single(90));
    %         netcdf.putAtt( ncid, varid_lat, 'coordinates', 'BEAR RNGE' );
    netcdf.putAtt(ncid, varid_lat, '_FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
    %     netcdf.putAtt(ncid, varid_lat, 'axis', 'Y');
    netcdf.putAtt(ncid, varid_lat, 'sdn_parameter_name', 'Latitude north');
    netcdf.putAtt(ncid, varid_lat, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
    netcdf.putAtt(ncid, varid_lat, 'sdn_uom_name', 'Degrees north');
    netcdf.putAtt(ncid, varid_lat, 'sdn_uom_urn', 'SDN:P06::DEGN');
    netcdf.putAtt(ncid, varid_lat, 'grid_mapping', 'crs');
    netcdf.putAtt(ncid, varid_lat, 'ancillary_variables', 'POSITION_QC');
    
    % Longitude
    varid_lon = netcdf.defVar( ncid, 'LONGITUDE', 'float', [dimid_range dimid_bearing] );
    netcdf.defVarDeflate(ncid, varid_lon, true, true, 6);
    netcdf.putAtt( ncid, varid_lon, 'standard_name', 'longitude' );
    netcdf.putAtt( ncid, varid_lon, 'long_name', 'Longitude of each location' );
    netcdf.putAtt( ncid, varid_lon, 'units', 'degree_east' );
    %     netcdf.putAtt(ncid, varid_lon, 'valid_range', single( [-180 180]));
    netcdf.putAtt(ncid, varid_lon, 'valid_min', single(-180));
    netcdf.putAtt(ncid, varid_lon, 'valid_max', single(180));
    %         netcdf.putAtt( ncid, varid_lon, 'coordinates', 'BEAR RNGE' );
    netcdf.putAtt(ncid, varid_lon, '_FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
    netcdf.putAtt(ncid, varid_lon, 'sdn_parameter_name', 'Longitude east');
    netcdf.putAtt(ncid, varid_lon, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
    netcdf.putAtt(ncid, varid_lon, 'sdn_uom_name', 'Degrees east');
    netcdf.putAtt(ncid, varid_lon, 'sdn_uom_urn', 'SDN:P06::DEGE');
    netcdf.putAtt(ncid, varid_lon, 'grid_mapping', 'crs');
    netcdf.putAtt(ncid, varid_lon, 'ancillary_variables', 'POSITION_QC');
    
    % crs
    %     varid_crs = netcdf.defVar( ncid, 'crs', 'short', dimid_t);
    varid_crs = netcdf.defVar( ncid, 'crs', 'short', []);
    netcdf.defVarDeflate(ncid, varid_crs, true, true, 6);
    netcdf.putAtt( ncid, varid_crs, 'grid_mapping_name', 'latitude_longitude' );
    netcdf.putAtt( ncid, varid_crs, 'epsg_code', 'EPSG:4326' );
    netcdf.putAtt( ncid, varid_crs, 'semi_major_axis', double(6378137.0) );
    netcdf.putAtt(ncid, varid_crs, 'inverse_flattening', double(298.257223563));
    
    %%
    
    %% Add SDN namespace variables
    
    % To enforce homogeneity in the codes and interoperability with SDC, the site_code has to be set equal
    % to the EDIOS Series id of the HFR network and platform_code must include the EDIOS Platform id of the
    % HFR site, i.e.:
    %           SDN_CRUISE=site_code=EDIOS-Series-id
    %           SDN_STATION=platform_code=EDIOS-Series-id_Total (for total current data files)
    %           SDN_STATION=platform_code=EDIOS-Series-id_ EDIOS-Platform-id (for radial current data files)
    
    % SDN_CRUISE
    varid_sdncruise = netcdf.defVar( ncid, 'SDN_CRUISE', 'char', [dimid_string_site_code dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdncruise, true, true, 6);
    netcdf.putAtt( ncid, varid_sdncruise, 'long_name', 'Grid grouping label');
    
    % SDN_STATION
    varid_sdnstation = netcdf.defVar( ncid, 'SDN_STATION', 'char', [dimid_string_platform_code dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdnstation, true, true, 6);
    netcdf.putAtt( ncid, varid_sdnstation, 'long_name', 'Grid label');
    
    % SDN_LOCAL_CDI_ID
    varid_sdnlocalcdiid = netcdf.defVar( ncid, 'SDN_LOCAL_CDI_ID', 'char', [dimid_string_sdn_local_cdi_id dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdnlocalcdiid, true, true, 6);
    netcdf.putAtt( ncid, varid_sdnlocalcdiid, 'long_name', 'SeaDataCloud CDI identifier');
    netcdf.putAtt( ncid, varid_sdnlocalcdiid, 'cf_role', 'grid_id');
    
    % SDN_EDMO_CODE
    varid_sdnedmocode = netcdf.defVar( ncid, 'SDN_EDMO_CODE', 'short', [dimid_maxinst dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdnedmocode, true, true, 6);
    netcdf.putAtt( ncid, varid_sdnedmocode, 'long_name', 'European Directory of Marine Organisations code for the CDI partner');
    netcdf.putAtt( ncid, varid_sdnedmocode, 'units', char('1'));
    
    % SDN_REFERENCES
    varid_sdnreferences = netcdf.defVar( ncid, 'SDN_REFERENCES', 'char', [dimid_string_sdn_references dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdnreferences, true, true, 6);
    netcdf.putAtt( ncid, varid_sdnreferences, 'long_name', 'Usage metadata reference');
    
    % SDN_XLINK
    varid_sdnxlink = netcdf.defVar(ncid, 'SDN_XLINK', 'char', [dimid_string_sdn_xlink dimid_refmax dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sdnxlink, true, true, 6);
    netcdf.putAtt( ncid, varid_sdnxlink, 'long_name', 'External resource linkages');
    
    %%
    
    %% Add data variables
    
    % radial_sea_water_velocity_away_from_instrument:
    % A velocity is a vector quantity. Radial velocity away from instrument
    % means the component of the velocity along the line of sight of the
    % instrument where positive implies movement away from the instrument (i.e.
    % outward). The "instrument" (examples are radar and lidar) is the device
    % used to make an observation.
    varid_speed = netcdf.defVar(ncid, 'RDVA', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_speed, true, true, 6);
    %     netcdf.putAtt(ncid, varid_speed, 'valid_range', int16( [(-10-addOffset)./scaleFactor (10-addOffset)./scaleFactor]));
    netcdf.putAtt(ncid, varid_speed, 'valid_min', int16( (-10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_speed, 'valid_max', int16( (10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_speed, 'standard_name', 'radial_sea_water_velocity_away_from_instrument');
    netcdf.putAtt(ncid, varid_speed, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_speed, 'standard_name', 'radial_sea_water_velocity_away_from_instrument');
    netcdf.putAtt(ncid, varid_speed, 'long_name', 'Radial sea water velocity away from instrument');
    netcdf.putAtt(ncid, varid_speed, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_speed, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_speed, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_speed, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
    netcdf.putAtt(ncid, varid_speed, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
    netcdf.putAtt(ncid, varid_speed, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_speed, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt( ncid, varid_speed, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_speed, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, RDCT_QC');
    
    % (radial current) direction
    % direction_of_radial_vector_away_from_instrument:
    % The direction_of_radial_vector_away_from_instrument is the direction in
    % which the instrument itself is pointing. The direction is measured
    % positive clockwise from due north. The "instrument" (examples are radar
    % and lidar) is the device used to make an observation. "direction_of_X"
    % means direction of a vector, a bearing.
    varid_direction = netcdf.defVar(ncid, 'DRVA', 'int', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_direction, true, true, 6);
    %     netcdf.putAtt(ncid, varid_direction, 'valid_range', int32( [0 (360-addOffset)./scaleFactor]));
    netcdf.putAtt(ncid, varid_direction, 'valid_min', int32(0));
    netcdf.putAtt(ncid, varid_direction, 'valid_max', int32((360-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_direction, 'standard_name', 'direction_of_radial_vector_away_from_instrument');
    netcdf.putAtt(ncid, varid_direction, 'long_name', 'Direction of radial vector away from instrument');
    netcdf.putAtt(ncid, varid_direction, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_direction, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_direction, 'units', 'degree_true');
    netcdf.putAtt(ncid, varid_direction, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_direction, 'sdn_parameter_name', 'Current direction (Eulerian) in the water body by directional range-gated radar');
    netcdf.putAtt(ncid, varid_direction, 'sdn_parameter_urn', 'SDN:P01::LCDAWVRD');
    netcdf.putAtt(ncid, varid_direction, 'sdn_uom_name', 'Degrees True');
    netcdf.putAtt(ncid, varid_direction, 'sdn_uom_urn', 'SDN:P06::UABB');
    netcdf.putAtt( ncid, varid_direction, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_direction, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, AVRB_QC, RDCT_QC');
    
    % u
    varid_u = netcdf.defVar(ncid, 'EWCT', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_u, true, true, 6);
    %     netcdf.putAtt(ncid, varid_u, 'valid_range', int16( [(-10-addOffset)./scaleFactor (10-addOffset)./scaleFactor]));
    netcdf.putAtt(ncid, varid_u, 'valid_min', int16( (-10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_u, 'valid_max', int16( (10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_u, 'standard_name', 'eastward_sea_water_velocity');
    netcdf.putAtt(ncid, varid_u, 'long_name', 'West-east current component');
    netcdf.putAtt(ncid, varid_u, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_u, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_u, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_u, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_u, 'sdn_parameter_name', 'Eastward current velocity in the water body');
    netcdf.putAtt(ncid, varid_u, 'sdn_parameter_urn', 'SDN:P01::LCEWZZ01');
    netcdf.putAtt(ncid, varid_u, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_u, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_u, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_u, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC');
    
    % v
    varid_v = netcdf.defVar(ncid, 'NSCT', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_v, true, true, 6);
    %     netcdf.putAtt(ncid, varid_v, 'valid_range', int16( [(-10-addOffset)./scaleFactor (10-addOffset)./scaleFactor]));
    netcdf.putAtt(ncid, varid_v, 'valid_min', int16( (-10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_v, 'valid_max', int16( (10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_v, 'standard_name', 'northward_sea_water_velocity');
    netcdf.putAtt(ncid, varid_v, 'long_name', 'South-north current component');
    netcdf.putAtt(ncid, varid_v, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_v, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_v, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_v, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_v, 'sdn_parameter_name', 'Northward current velocity in the water body');
    netcdf.putAtt(ncid, varid_v, 'sdn_parameter_urn', 'SDN:P01::LCNSZZ01');
    netcdf.putAtt(ncid, varid_v, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_v, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_v, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_v, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC');
    
    % Spatial Quality
    varid_espc = netcdf.defVar(ncid, 'ESPC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_espc, true, true, 6);
    %     netcdf.putAtt(ncid, varid_espc, 'valid_range', int16( [(-32-addOffset)./scaleFactor (32-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_espc, 'valid_min', int16( (-32-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_espc, 'valid_max', int16( (32-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_espc, 'long_name', 'Radial standard deviation of current velocity over the scatter patch');
    netcdf.putAtt(ncid, varid_espc, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_espc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_espc, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_espc, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_espc, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_espc, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_espc, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_espc, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_espc, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_espc, 'ancillary_variables', 'QCflag, VART_QC');
    
    % Temporal Quality
    varid_etmp = netcdf.defVar(ncid, 'ETMP', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_etmp, true, true, 6);
    %     netcdf.putAtt(ncid, varid_etmp, 'valid_range', int16( [(-32-addOffset)./scaleFactor (32-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_etmp, 'valid_min', int16( (-32-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_etmp, 'valid_max', int16( (32-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_etmp, 'long_name', 'Radial standard deviation of current velocity over coverage period');
    netcdf.putAtt(ncid, varid_etmp, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_etmp, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_etmp, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_etmp, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_etmp, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_etmp, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_etmp, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_etmp, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_etmp, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_etmp, 'ancillary_variables', 'QCflag, VART_QC');
    
    % Velocity Maximum
    varid_maxv = netcdf.defVar(ncid, 'MAXV', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_maxv, true, true, 6);
    netcdf.putAtt(ncid, varid_maxv, 'long_name', 'Radial sea water velocity away from instrument maximum');
    %     netcdf.putAtt(ncid, varid_maxv, 'valid_range', int16( [(-10-addOffset)./scaleFactor (10-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_maxv, 'valid_min', int16( (-10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_maxv, 'valid_max', int16( (10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_maxv, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_maxv, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_maxv, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_maxv, 'units', 'm s-1' );
    netcdf.putAtt(ncid, varid_maxv, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
    netcdf.putAtt(ncid, varid_maxv, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
    netcdf.putAtt(ncid, varid_maxv, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_maxv, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_maxv, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_maxv, 'ancillary_variables', 'QCflag, MDFL_QC, CSPD_QC, VART_QC');
    
    % Velocity Minimum
    varid_minv = netcdf.defVar(ncid, 'MINV', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_minv, true, true, 6);
    netcdf.putAtt( ncid, varid_minv, 'long_name','Radial sea water velocity away from instrument minimum');
    %     netcdf.putAtt(ncid, varid_minv, 'valid_range', int16( [(-10-addOffset)./scaleFactor (10-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_minv, 'valid_min', int16( (-10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_minv, 'valid_max', int16( (10-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_minv, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_minv, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_minv, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_minv, 'units', 'm s-1');
    netcdf.putAtt(ncid, varid_minv, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
    netcdf.putAtt(ncid, varid_minv, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
    netcdf.putAtt(ncid, varid_minv, 'sdn_uom_name', 'Metres per second');
    netcdf.putAtt(ncid, varid_minv, 'sdn_uom_urn', 'SDN:P06::UVAA');
    netcdf.putAtt(ncid, varid_minv, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_minv, 'ancillary_variables', 'QCflag, MDFL_QC, CSPD_QC, VART_QC');
    
    % Spatial Count
    varid_ersc = netcdf.defVar(ncid, 'ERSC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_ersc, true, true, 6);
    netcdf.putAtt(ncid, varid_ersc, 'long_name', 'Radial sea water velocity spatial quality count' );
    %     netcdf.putAtt(ncid, varid_ersc, 'valid_range', int16( [0 127] ));
    netcdf.putAtt(ncid, varid_ersc, 'valid_min', int16(0));
    netcdf.putAtt(ncid, varid_ersc, 'valid_max', int16(127));
    netcdf.putAtt(ncid, varid_ersc, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_ersc, 'scale_factor', int16(1));
    netcdf.putAtt(ncid, varid_ersc, 'add_offset', int16(0));
    netcdf.putAtt(ncid, varid_ersc, 'units', '1');
    netcdf.putAtt(ncid, varid_ersc, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_ersc, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_ersc, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_ersc, 'sdn_uom_urn', 'SDN:P06::UUUU');
    netcdf.putAtt(ncid, varid_ersc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_ersc, 'ancillary_variables', 'QCflag');
    
    % Temporal Count
    varid_ertc = netcdf.defVar(ncid, 'ERTC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_ertc, true, true, 6);
    netcdf.putAtt(ncid, varid_ertc, 'long_name', 'Radial sea water velocity temporal quality count');
    %     netcdf.putAtt(ncid, varid_ertc, 'valid_range', int16( [0 127] ));
    netcdf.putAtt(ncid, varid_ertc, 'valid_min', int16(0));
    netcdf.putAtt(ncid, varid_ertc, 'valid_max', int16(127));
    netcdf.putAtt(ncid, varid_ertc, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_ertc, 'scale_factor', int16(1));
    netcdf.putAtt(ncid, varid_ertc, 'add_offset', int16(0));
    netcdf.putAtt(ncid, varid_ertc, 'units', '1');
    netcdf.putAtt(ncid, varid_ertc, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_ertc, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_ertc, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_ertc, 'sdn_uom_urn', 'SDN:P06::UUUU');
    netcdf.putAtt(ncid, varid_ertc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_ertc, 'ancillary_variables', 'QCflag');
    
    % X-Distance
    varid_xdst = netcdf.defVar(ncid, 'XDST', 'int', [dimid_range dimid_bearing]);
    netcdf.defVarDeflate(ncid, varid_xdst, true, true, 6);
    netcdf.putAtt(ncid, varid_xdst, 'long_name', 'Eastward distance from instrument');
    %     netcdf.putAtt(ncid, varid_xdst, 'valid_range', int32( [0 (1e3-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_xdst, 'valid_min', int32(0));
    netcdf.putAtt(ncid, varid_xdst, 'valid_max', int32((1e3-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_xdst, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_xdst, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_xdst, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_xdst, 'units', 'km');
    netcdf.putAtt(ncid, varid_xdst, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_xdst, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_xdst, 'sdn_uom_name', 'Kilometres');
    netcdf.putAtt(ncid, varid_xdst, 'sdn_uom_urn', 'SDN:P06::ULKM');
    netcdf.putAtt(ncid, varid_xdst, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_xdst, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
    
    % Y-Distance
    varid_ydst = netcdf.defVar(ncid, 'YDST', 'int', [dimid_range dimid_bearing]);
    netcdf.defVarDeflate(ncid, varid_ydst, true, true, 6);
    netcdf.putAtt(ncid, varid_ydst, 'long_name', 'Northward distance from instrument');
    %     netcdf.putAtt(ncid, varid_ydst, 'valid_range', int32( [0 (1e3-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_ydst, 'valid_min', int32(0));
    netcdf.putAtt(ncid, varid_ydst, 'valid_max', int32((1e3-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_ydst, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_ydst, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_ydst, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_ydst, 'units', 'km');
    netcdf.putAtt(ncid, varid_ydst, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_ydst, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_ydst, 'sdn_uom_name', 'Kilometres');
    netcdf.putAtt(ncid, varid_ydst, 'sdn_uom_urn', 'SDN:P06::ULKM');
    netcdf.putAtt(ncid, varid_ydst, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_ydst, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
    
    % Spectra Range Cell
    varid_sprc = netcdf.defVar(ncid, 'SPRC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sprc, true, true, 6);
    netcdf.putAtt(ncid, varid_sprc, 'long_name', 'Radial sea water velocity cross spectra range cell');
    %     netcdf.putAtt(ncid, varid_sprc, 'valid_range', int16( [0 127] ));
    netcdf.putAtt(ncid, varid_sprc, 'valid_min', int16(0));
    netcdf.putAtt(ncid, varid_sprc, 'valid_max', int16(127));
    netcdf.putAtt(ncid, varid_sprc, '_FillValue', netcdf.getConstant('NC_FILL_SHORT'));
    netcdf.putAtt(ncid, varid_sprc, 'scale_factor', int16(1));
    netcdf.putAtt(ncid, varid_sprc, 'add_offset', int16(0));
    netcdf.putAtt(ncid, varid_sprc, 'units', '1');
    netcdf.putAtt(ncid, varid_sprc, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_sprc, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_sprc, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_sprc, 'sdn_uom_urn', 'SDN:P06::UUUU');
    netcdf.putAtt(ncid, varid_sprc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    netcdf.putAtt(ncid, varid_sprc, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
    
    % Number of receive antennas
    varid_narx = netcdf.defVar(ncid, 'NARX', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_narx, true, true, 6);
    netcdf.putAtt(ncid, varid_narx, 'long_name', 'Number of receive antennas');
    %     netcdf.putAtt(ncid, varid_narx, 'valid_range', int8([0 maxsite] ));
    netcdf.putAtt(ncid, varid_narx, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_narx, 'valid_max', int8(maxsite));
    netcdf.putAtt(ncid, varid_narx, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_narx, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_narx, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_narx, 'units', char('1'));
    netcdf.putAtt(ncid, varid_narx, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_narx, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_narx, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_narx, 'sdn_uom_urn', 'SDN:P06::UUUU');
    %         netcdf.putAtt(ncid, varid_sprc, 'coordinates', 'TIME' );
    
    % Number of transmit antennas
    varid_natx = netcdf.defVar(ncid, 'NATX', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_natx, true, true, 6);
    netcdf.putAtt(ncid, varid_natx, 'long_name', 'Number of transmit antennas');
    %     netcdf.putAtt(ncid, varid_natx, 'valid_range', int8([0 maxsite] ));
    netcdf.putAtt(ncid, varid_natx, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_natx, 'valid_max', int8(maxsite));
    netcdf.putAtt(ncid, varid_natx, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_natx, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_natx, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_natx, 'units', char('1'));
    netcdf.putAtt(ncid, varid_natx, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_natx, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_natx, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_natx, 'sdn_uom_urn', 'SDN:P06::UUUU');
    %         netcdf.putAtt(ncid, varid_natx, 'coordinates', 'TIME' );
    
    % Receive antenna latitudes
    varid_sltr = netcdf.defVar(ncid, 'SLTR', 'int', [dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sltr, true, true, 6);
    netcdf.putAtt(ncid, varid_sltr, 'long_name', 'Receive antenna latitudes');
    netcdf.putAtt( ncid, varid_sltr, 'standard_name', 'latitude' );
    %     netcdf.putAtt(ncid, varid_sltr, 'valid_range', int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_sltr, 'valid_min', int32((-90-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_sltr, 'valid_max', int32((90-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_sltr, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_sltr, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_sltr, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_sltr, 'units', 'degree_north');
    netcdf.putAtt(ncid, varid_sltr, 'sdn_parameter_name', 'Latitude north');
    netcdf.putAtt(ncid, varid_sltr, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
    netcdf.putAtt(ncid, varid_sltr, 'sdn_uom_name', 'Degrees north');
    netcdf.putAtt(ncid, varid_sltr, 'sdn_uom_urn', 'SDN:P06::DEGN');
    netcdf.putAtt(ncid, varid_sltr, 'coordinates', 'TIME MAXSITE' );
    
    % Receive antenna longitudes
    varid_slnr = netcdf.defVar(ncid, 'SLNR', 'int', [dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_slnr, true, true, 6);
    netcdf.putAtt(ncid, varid_slnr, 'long_name', 'Receive antenna longitudes');
    netcdf.putAtt( ncid, varid_slnr, 'standard_name', 'longitude' );
    %     netcdf.putAtt(ncid, varid_slnr, 'valid_range', int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_slnr, 'valid_min', int32((-180-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_slnr, 'valid_max', int32((180-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_slnr, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_slnr, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_slnr, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_slnr, 'units', 'degree_east');
    netcdf.putAtt(ncid, varid_slnr, 'sdn_parameter_name', 'Longitude east');
    netcdf.putAtt(ncid, varid_slnr, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
    netcdf.putAtt(ncid, varid_slnr, 'sdn_uom_name', 'Degrees east');
    netcdf.putAtt(ncid, varid_slnr, 'sdn_uom_urn', 'SDN:P06::DEGE');
    netcdf.putAtt(ncid, varid_slnr, 'coordinates', 'TIME MAXSITE' );
    
    % Transmit antenna latitudes
    varid_sltt = netcdf.defVar(ncid, 'SLTT', 'int', [dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_sltt, true, true, 6);
    netcdf.putAtt(ncid, varid_sltt, 'long_name', 'Transmit antenna latitudes');
    netcdf.putAtt( ncid, varid_sltt, 'standard_name', 'latitude' );
    %     netcdf.putAtt(ncid, varid_sltt, 'valid_range', int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_sltt, 'valid_min', int32((-90-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_sltt, 'valid_max', int32((90-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_sltt, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_sltt, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_sltt, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_sltt, 'units', 'degree_north');
    netcdf.putAtt(ncid, varid_sltt, 'sdn_parameter_name', 'Latitude north');
    netcdf.putAtt(ncid, varid_sltt, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
    netcdf.putAtt(ncid, varid_sltt, 'sdn_uom_name', 'Degrees north');
    netcdf.putAtt(ncid, varid_sltt, 'sdn_uom_urn', 'SDN:P06::DEGN');
    netcdf.putAtt(ncid, varid_sltt, 'coordinates', 'TIME MAXSITE' );
    
    % Transmit antenna longitudes
    varid_slnt = netcdf.defVar(ncid, 'SLNT', 'int', [dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_slnt, true, true, 6);
    netcdf.putAtt(ncid, varid_slnt, 'long_name', 'Transmit antenna longitudes');
    netcdf.putAtt( ncid, varid_slnt, 'standard_name', 'longitude' );
    %     netcdf.putAtt(ncid, varid_slnt, 'valid_range', int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    netcdf.putAtt(ncid, varid_slnt, 'valid_min', int32((-180-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_slnt, 'valid_max', int32((180-addOffset)./scaleFactor));
    netcdf.putAtt(ncid, varid_slnt, '_FillValue', netcdf.getConstant('NC_FILL_INT'));
    netcdf.putAtt(ncid, varid_slnt, 'scale_factor', single(scaleFactor));
    netcdf.putAtt(ncid, varid_slnt, 'add_offset', single(addOffset));
    netcdf.putAtt(ncid, varid_slnt, 'units', 'degree_east');
    netcdf.putAtt(ncid, varid_slnt, 'sdn_parameter_name', 'Longitude east');
    netcdf.putAtt(ncid, varid_slnt, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
    netcdf.putAtt(ncid, varid_slnt, 'sdn_uom_name', 'Degrees east');
    netcdf.putAtt(ncid, varid_slnt, 'sdn_uom_urn', 'SDN:P06::DEGE');
    netcdf.putAtt(ncid, varid_slnt, 'coordinates', 'TIME MAXSITE' );
    
    % Receive antenna codes
    varid_scdr = netcdf.defVar(ncid, 'SCDR', 'char', [dimid_string_siteCodeArr dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_scdr, true, true, 6);
    netcdf.putAtt(ncid, varid_scdr, 'long_name', 'Receive antenna codes');
    netcdf.putAtt(ncid, varid_scdr, 'units', '1');
    netcdf.putAtt(ncid, varid_scdr, '_FillValue', netcdf.getConstant('NC_FILL_CHAR'));
    netcdf.putAtt(ncid, varid_scdr, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_scdr, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_scdr, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_scdr, 'sdn_uom_urn', 'SDN:P06::UUUU');
    
    % Transmit antenna codes
    varid_scdt = netcdf.defVar(ncid, 'SCDT', 'char', [dimid_string_siteCodeArr dimid_maxsite dimid_t]);
    netcdf.defVarDeflate(ncid, varid_scdt, true, true, 6);
    netcdf.putAtt(ncid, varid_scdt, 'long_name', 'Transmit antenna codes');
    netcdf.putAtt(ncid, varid_scdt, 'units', '1');
    netcdf.putAtt(ncid, varid_scdt, '_FillValue', netcdf.getConstant('NC_FILL_CHAR'));
    netcdf.putAtt(ncid, varid_scdt, 'sdn_parameter_name', '');
    netcdf.putAtt(ncid, varid_scdt, 'sdn_parameter_urn', '');
    netcdf.putAtt(ncid, varid_scdt, 'sdn_uom_name', 'Dimensionless');
    netcdf.putAtt(ncid, varid_scdt, 'sdn_uom_urn', 'SDN:P06::UUUU');
    
    %%
    
    %% Add QC variables
    
    % Time QC Flag
    varid_tqc = netcdf.defVar(ncid, 'TIME_QC', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_tqc, true, true, 6);
    netcdf.putAtt(ncid, varid_tqc, 'long_name', 'Time quality flag');
    netcdf.putAtt(ncid, varid_tqc, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_tqc, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_tqc, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_tqc, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_tqc, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_tqc, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_tqc, 'comment', 'OceanSITES quality flagging for temporal coordinate.');
    netcdf.putAtt(ncid, varid_tqc, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_tqc, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_tqc, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_tqc, 'units', '1');
    
    % Position QC Flag
    varid_posqc = netcdf.defVar(ncid, 'POSITION_QC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_posqc, true, true, 6);
    netcdf.putAtt(ncid, varid_posqc, 'long_name', 'Position quality flags');
    netcdf.putAtt(ncid, varid_posqc, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_posqc, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_posqc, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_posqc, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_posqc, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_posqc, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_posqc, 'comment', 'OceanSITES quality flagging for position coordinates');
    netcdf.putAtt(ncid, varid_posqc, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_posqc, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_posqc, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_posqc, 'units', '1');
    netcdf.putAtt(ncid, varid_posqc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Depth QC Flag
    varid_dqc = netcdf.defVar(ncid, 'DEPH_QC', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_dqc, true, true, 6);
    netcdf.putAtt(ncid, varid_dqc, 'long_name', 'Depth quality flag');
    netcdf.putAtt(ncid, varid_dqc, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_dqc, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_dqc, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_dqc, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_dqc, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_dqc, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_dqc, 'comment', 'OceanSITES quality flagging for depth coordinate.');
    netcdf.putAtt(ncid, varid_dqc, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_dqc, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_dqc, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_dqc, 'units', '1');
    
    % Overall QC Flag
    varid_ovqc = netcdf.defVar(ncid, 'QCflag', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_ovqc, true, true, 6);
    netcdf.putAtt(ncid, varid_ovqc, 'long_name', 'Overall quality flags');
    netcdf.putAtt(ncid, varid_ovqc, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_ovqc, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_ovqc, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_ovqc, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_ovqc, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_ovqc, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_ovqc, 'comment', 'OceanSITES quality flagging for all QC tests.');
    netcdf.putAtt(ncid, varid_ovqc, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_ovqc, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_ovqc, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_ovqc, 'units', '1');
    netcdf.putAtt(ncid, varid_ovqc, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Over-water QC Flag
    varid_owtr = netcdf.defVar(ncid, 'OWTR_QC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_owtr, true, true, 6);
    netcdf.putAtt(ncid, varid_owtr, 'long_name', 'Over-water quality flags');
    netcdf.putAtt(ncid, varid_owtr, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_owtr, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_owtr, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_owtr, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_owtr, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_owtr, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_owtr, 'comment', 'OceanSITES quality flagging for Over-water QC test.');
    netcdf.putAtt(ncid, varid_owtr, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_owtr, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_owtr, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_owtr, 'units', '1');
    netcdf.putAtt(ncid, varid_owtr, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Median Filter QC Flag
    varid_mdfl = netcdf.defVar(ncid, 'MDFL_QC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_mdfl, true, true, 6);
    netcdf.putAtt(ncid, varid_mdfl, 'long_name', 'Median filter quality flags');
    netcdf.putAtt(ncid, varid_mdfl, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_mdfl, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_mdfl, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_mdfl, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_mdfl, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_mdfl, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_mdfl, 'comment', ['OceanSITES quality flagging for Median Filter QC test. Threshold set to ' num2str(Radial_QC_params.MedFilt(1)) ' km, ' num2str(Radial_QC_params.MedFilt(2)) ' deg, ' num2str(Radial_QC_params.MedFilt(3)) ' m/s, ']);
    netcdf.putAtt(ncid, varid_mdfl, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_mdfl, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_mdfl, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_mdfl, 'units', '1');
    netcdf.putAtt(ncid, varid_mdfl, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Variance Threshold QC Flag
    varid_vart = netcdf.defVar(ncid, 'VART_QC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_vart, true, true, 6);
    netcdf.putAtt(ncid, varid_vart, 'long_name', 'Variance threshold quality flags');
    netcdf.putAtt(ncid, varid_vart, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_vart, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_vart, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_vart, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_vart, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_vart, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_vart, 'comment', ['OceanSITES quality flagging for Variance Threshold QC test. Test not applicable to Direction Finding systems. The Temporal Derivative test is applied. Threshold set to ' num2str(Radial_QC_params.TempDerThr.threshold) ' m/s.']);
    netcdf.putAtt(ncid, varid_vart, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_vart, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_vart, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_vart, 'units', '1');
    netcdf.putAtt(ncid, varid_vart, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Velocity Threshold QC Flag
    varid_cspd = netcdf.defVar(ncid, 'CSPD_QC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
    netcdf.defVarDeflate(ncid, varid_cspd, true, true, 6);
    netcdf.putAtt(ncid, varid_cspd, 'long_name', 'Velocity threshold quality flags');
    netcdf.putAtt(ncid, varid_cspd, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_cspd, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_cspd, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_cspd, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_cspd, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_cspd, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_cspd, 'comment', ['OceanSITES quality flagging for Velocity Threshold QC test. Threshold set to ' num2str(Radial_QC_params.VelThr) ' m/s.']);
    netcdf.putAtt(ncid, varid_cspd, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_cspd, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_cspd, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_cspd, 'units', '1');
    netcdf.putAtt(ncid, varid_cspd, 'coordinates', 'TIME DEPH LATITUDE LONGITUDE' );
    
    % Average Radial Bearing QC Flag
    varid_avrb = netcdf.defVar(ncid, 'AVRB_QC', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_avrb, true, true, 6);
    netcdf.putAtt(ncid, varid_avrb, 'long_name', 'Average radial bearing quality flag');
    netcdf.putAtt(ncid, varid_avrb, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_avrb, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_avrb, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_avrb, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_avrb, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_avrb, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_avrb, 'comment', ['OceanSITES quality flagging for Average Radial Bearing QC test. Thresholds set to [' num2str(Radial_QC_params.AvgRadBear(1)) '-' num2str(Radial_QC_params.AvgRadBear(2)) '] deg.']);
    netcdf.putAtt(ncid, varid_avrb, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_avrb, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_avrb, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_avrb, 'units', '1');
    
    % Radial Count QC Flag
    varid_rdct = netcdf.defVar(ncid, 'RDCT_QC', 'byte', dimid_t);
    netcdf.defVarDeflate(ncid, varid_rdct, true, true, 6);
    netcdf.putAtt(ncid, varid_rdct, 'long_name', 'Radial count quality flag');
    netcdf.putAtt(ncid, varid_rdct, 'conventions', 'Copernicus Marine in situ reference table 2');
    %     netcdf.putAtt(ncid, varid_rdct, 'valid_range', int8( [0 9]));
    netcdf.putAtt(ncid, varid_rdct, 'valid_min', int8(0));
    netcdf.putAtt(ncid, varid_rdct, 'valid_max', int8(9));
    netcdf.putAtt(ncid, varid_rdct, 'flag_values', int8( [0 1 2 3 4 5 6 7 8 9]));
    netcdf.putAtt(ncid, varid_rdct, 'flag_meanings', 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
    netcdf.putAtt(ncid, varid_rdct, 'comment', ['OceanSITES quality flagging for Radial Count QC test. Thresholds set to ' num2str(Radial_QC_params.RadCnt) ' vectors.']);
    netcdf.putAtt(ncid, varid_rdct, '_FillValue', netcdf.getConstant('NC_FILL_BYTE'));
    netcdf.putAtt(ncid, varid_rdct, 'scale_factor', int8(1));
    netcdf.putAtt(ncid, varid_rdct, 'add_offset', int8(0));
    netcdf.putAtt(ncid, varid_rdct, 'units', '1');
    
    %%
    
    %% Add global attributes
    
    varid_global = netcdf.getConstant('global');
    
    % MANDATORY ATTRIBUTES
    
    % Discovery and Identification
    netcdf.putAtt(ncid, varid_global, 'site_code', site_code);
    netcdf.putAtt(ncid, varid_global, 'platform_code', platform_code);
    netcdf.putAtt(ncid, varid_global, 'data_mode', 'R');
    DoAIndex = find(not(cellfun('isempty', strfind(stationFields, 'DoA_estimation_method'))));
    netcdf.putAtt(ncid, varid_global, 'DoA_estimation_method', stationData{DoAIndex});
    calibration_typeIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_type'))));
    netcdf.putAtt(ncid, varid_global, 'calibration_type', stationData{calibration_typeIndex});
    netcdf.putAtt(ncid, varid_global, 'last_calibration_date', lastPatternStr);
    calibration_linkIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_link'))));
    netcdf.putAtt(ncid, varid_global, 'calibration_link', stationData{calibration_linkIndex});
    titleIndex = find(not(cellfun('isempty', strfind(networkFields, 'title'))));
    netcdf.putAtt(ncid, varid_global, 'title', networkData{titleIndex});
    summaryIndex = find(not(cellfun('isempty', strfind(stationFields, 'summary'))));
    netcdf.putAtt(ncid, varid_global, 'summary', stationData{summaryIndex});
    netcdf.putAtt(ncid, varid_global, 'source', 'coastal structure');
    netcdf.putAtt(ncid, varid_global, 'source_platform_category_code', '17');
    institution_nameIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_name'))));
    netcdf.putAtt(ncid, varid_global, 'institution', stationData{institution_nameIndex});
    netcdf.putAtt(ncid, varid_global, 'institution_edmo_code', num2str(EDMO_code));
    netcdf.putAtt(ncid, varid_global, 'data_assembly_center', 'European HFR Node');
    netcdf.putAtt(ncid, varid_global, 'id', id);
    % Geo-spatial-temporal
    netcdf.putAtt(ncid, varid_global, 'data_type', 'HF radar radial data');
    netcdf.putAtt(ncid, varid_global, 'feature_type', 'surface');
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    netcdf.putAtt(ncid, varid_global, 'geospatial_lat_min', num2str(networkData{geospatial_lat_minIndex}));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    netcdf.putAtt(ncid, varid_global, 'geospatial_lat_max', num2str(networkData{geospatial_lat_maxIndex}));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    netcdf.putAtt(ncid, varid_global, 'geospatial_lon_min', num2str(networkData{geospatial_lon_minIndex}));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    netcdf.putAtt(ncid, varid_global, 'geospatial_lon_max', num2str(networkData{geospatial_lon_maxIndex}));
    netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_max', num2str(vertMax));
    netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_min', '0');
    netcdf.putAtt(ncid, varid_global, 'time_coverage_start', timeCoverageStart);
    netcdf.putAtt(ncid, varid_global, 'time_coverage_end', timeCoverageEnd);
    % Conventions used
    netcdf.putAtt(ncid, varid_global, 'format_version', 'v2.1.2');
    netcdf.putAtt(ncid, varid_global, 'Conventions', 'CF-1.6 Jerico-Next-deliverable-D5.14 Copernicus-InSituTAC-Manual-1.0 Copernicus-InSituTAC-SRD-1.4 Copernicus-InSituTAC-ParametersList-3.1.0');
    % Publication information
    netcdf.putAtt(ncid, varid_global, 'update_interval', 'void');
    netcdf.putAtt(ncid, varid_global, 'citation', citation_str);
    netcdf.putAtt(ncid, varid_global, 'distribution_statement', distribution_str);
    netcdf.putAtt(ncid, varid_global, 'publisher_name', 'European HFR Node');
    netcdf.putAtt(ncid, varid_global, 'publisher_url', 'http://eurogoos.eu/');
    netcdf.putAtt(ncid, varid_global, 'publisher_email', 'euhfrnode@azti.es');
    licenseIndex = find(not(cellfun('isempty', strfind(networkFields, 'license'))));
    netcdf.putAtt(ncid, varid_global, 'license', networkData{licenseIndex});
    acknowledgmentIndex = find(not(cellfun('isempty', strfind(networkFields, 'acknowledgment'))));
    netcdf.putAtt(ncid, varid_global, 'acknowledgment', networkData{acknowledgmentIndex});
    % Provenance
    netcdf.putAtt(ncid, varid_global, 'date_created', dateCreated);
    netcdf.putAtt(ncid, varid_global, 'history', hist_create);
    netcdf.putAtt(ncid, varid_global, 'date_modified', dateCreated);
    netcdf.putAtt(ncid, varid_global, 'date_update', dateCreated);
    netcdf.putAtt(ncid, varid_global, 'processing_level', '2B');
    contributor_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_name'))));
    netcdf.putAtt(ncid, varid_global, 'contributor_name', networkData{contributor_nameIndex});
    contributor_roleIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_role'))));
    netcdf.putAtt(ncid, varid_global, 'contributor_role', networkData{contributor_roleIndex});
    contributor_emailIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_email'))));
    netcdf.putAtt(ncid, varid_global, 'contributor_email', networkData{contributor_emailIndex});
    
    % RECOMMENDED ATTRIBUTES
    
    % Discovery and Identification
    projectIndex = find(not(cellfun('isempty', strfind(networkFields, 'project'))));
    netcdf.putAtt(ncid, varid_global, 'project', networkData{projectIndex});
    netcdf.putAtt(ncid, varid_global, 'naming_authority', naming_authorityStr);
    netcdf.putAtt(ncid, varid_global, 'keywords', 'OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF');
    netcdf.putAtt(ncid, varid_global, 'keywords_vocabulary', 'GCMD Science Keywords');
    commentIndex = find(not(cellfun('isempty', strfind(networkFields, 'comment'))));
    netcdf.putAtt(ncid, varid_global, 'comment', networkData{commentIndex});
    netcdf.putAtt(ncid, varid_global, 'data_language', 'eng');
    netcdf.putAtt(ncid, varid_global, 'data_character_set', 'utf8');
    netcdf.putAtt(ncid, varid_global, 'metadata_language', 'eng');
    netcdf.putAtt(ncid, varid_global, 'metadata_character_set', 'utf8');
    netcdf.putAtt(ncid, varid_global, 'topic_category', 'oceans');
    network_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_name'))));
    netcdf.putAtt(ncid, varid_global, 'network', networkData{network_nameIndex});
    % Geo-spatial-temporal
    areaIndex = find(not(cellfun('isempty', strfind(networkFields, 'area'))));
    netcdf.putAtt(ncid, varid_global, 'area', networkData{areaIndex});
    netcdf.putAtt(ncid, varid_global, 'geospatial_lat_units', 'degree_north');
    netcdf.putAtt(ncid, varid_global, 'geospatial_lon_units', 'degree_east');
    netcdf.putAtt(ncid, varid_global, 'geospatial_lat_resolution', num2str(latRes));
    netcdf.putAtt(ncid, varid_global, 'geospatial_lon_resolution', num2str(lonRes));
    netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_resolution', num2str(vertMax));
    netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_units', 'm');
    netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_positive', 'down');
    netcdf.putAtt(ncid, varid_global, 'time_coverage_duration', timeCoverageDuration);
    netcdf.putAtt(ncid, varid_global, 'time_coverage_resolution', timeCoverageResolution);
    netcdf.putAtt(ncid, varid_global, 'reference_system', 'EPSG:4806');
    netcdf.putAtt(ncid, varid_global, 'cdm_data_type', 'Grid');
    % Conventions used
    netcdf.putAtt(ncid, varid_global, 'netcdf_version', netcdf.inqLibVers);
    netcdf.putAtt(ncid, varid_global, 'netcdf_format', 'netcdf4_classic');
    
    % OTHER ATTRIBUTES
    netcdf.putAtt(ncid, varid_global, 'metadata_contact', 'lorenzo.corgnati@sp.ismar.cnr.it');
    netcdf.putAtt(ncid, varid_global, 'metadata_date_stamp', dateCreated);
    netcdf.putAtt(ncid, varid_global, 'standard_name_vocabulary', 'NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6');
    netcdf.putAtt(ncid, varid_global, 'sensor', 'CODAR SeaSonde');
    netcdf.putAtt(ncid, varid_global, 'institution_reference', institution_websiteStr);
    netcdf.putAtt(ncid, varid_global, 'references', 'High Frequency Radar European common data and metadata model Reference Card: all you need to know about High Frequency Radar (HFR) data harmonization at a glance. http://www.marineinsitu.eu/wp-content/uploads/2018/02/HFR_Data_Model_Reference_Card_v1.pdf');
    netcdf.putAtt(ncid, varid_global, 'software_name', 'HFR_Combiner');
    netcdf.putAtt(ncid, varid_global, 'software_version', 'v3.3');
    netcdf.putAtt(ncid, varid_global, 'date_issued', dateCreated);
    
    % Globals sourced from radial file metadata
    if (exist('UUID', 'var') == 0)
        UUID = '';
    end
    netcdf.putAtt(ncid, varid_global, 'UUID', UUID);
    
    if (exist('manufacturer', 'var') == 0)
        manufacturer = '';
    end
    netcdf.putAtt(ncid, varid_global, 'manufacturer', manufacturer);
    
    if (exist('RangeStart', 'var') == 0)
        RangeStart = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RangeStart', RangeStart);
    
    if (exist('RangeEnd', 'var') == 0)
        RangeEnd = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RangeEnd', RangeEnd);
    
    if (exist('RangeResolutionKMeters_str', 'var') == 0)
        RangeResolutionKMeters_str = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RangeResolutionKMeters', RangeResolutionKMeters_str);
    
    if (exist('RangeResolutionMeters_str', 'var') == 0)
        RangeResolutionMeters_str = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RangeResolutionMeters', RangeResolutionMeters_str);
    
    if (exist('AntennaBearing', 'var') == 0)
        AntennaBearing = '';
    end
    netcdf.putAtt(ncid, varid_global, 'AntennaBearing', AntennaBearing);
    
    if (exist('ReferenceBearing', 'var') == 0)
        ReferenceBearing = '';
    end
    netcdf.putAtt(ncid, varid_global, 'ReferenceBearing', ReferenceBearing);
    
    if (exist('AngularResolution_str', 'var') == 0)
        AngularResolution_str = '';
    end
    netcdf.putAtt(ncid, varid_global, 'AngularResolution', AngularResolution_str);
    
    if (exist('SpatialResolution', 'var') == 0)
        SpatialResolution = '';
    end
    netcdf.putAtt(ncid, varid_global, 'SpatialResolution', SpatialResolution);
    
    if (exist('PatternResolution', 'var') == 0)
        PatternResolution = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternResolution', PatternResolution);
    
    if (exist('TransmitCenterFreqMHz', 'var') == 0)
        TransmitCenterFreqMHz = '';
    end
    netcdf.putAtt(ncid, varid_global, 'TransmitCenterFreqMHz', TransmitCenterFreqMHz);
    
    if (exist('DopplerResolutionHzPerBin', 'var') == 0)
        DopplerResolutionHzPerBin = '';
    end
    netcdf.putAtt(ncid, varid_global, 'DopplerResolutionHzPerBin', DopplerResolutionHzPerBin);
    
    if (exist('FirstOrderMethod', 'var') == 0)
        FirstOrderMethod = '';
    end
    netcdf.putAtt(ncid, varid_global, 'FirstOrderMethod', FirstOrderMethod);
    
    if (exist('BraggSmoothingPoints', 'var') == 0)
        BraggSmoothingPoints = '';
    end
    netcdf.putAtt(ncid, varid_global, 'BraggSmoothingPoints', BraggSmoothingPoints);
    
    if (exist('BraggHasSecondOrder', 'var') == 0)
        BraggHasSecondOrder = '';
    end
    netcdf.putAtt(ncid, varid_global, 'BraggHasSecondOrder', BraggHasSecondOrder);
    
    if (exist('RadialBraggPeakDropOff', 'var') == 0)
        RadialBraggPeakDropOff = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RadialBraggPeakDropOff', RadialBraggPeakDropOff);
    
    if (exist('RadialBraggPeakNull', 'var') == 0)
        RadialBraggPeakNull = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RadialBraggPeakNull', RadialBraggPeakNull);
    
    if (exist('RadialBraggNoiseThreshold', 'var') == 0)
        RadialBraggNoiseThreshold = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RadialBraggNoiseThreshold', RadialBraggNoiseThreshold);
    
    if (exist('PatternAmplitudeCorrections', 'var') == 0)
        PatternAmplitudeCorrections = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternAmplitudeCorrections', PatternAmplitudeCorrections);
    
    if (exist('PatternPhaseCorrections', 'var') == 0)
        PatternPhaseCorrections = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternPhaseCorrections', PatternPhaseCorrections);
    
    if (exist('PatternAmplitudeCalculations', 'var') == 0)
        PatternAmplitudeCalculations = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternAmplitudeCalculations', PatternAmplitudeCalculations);
    
    if (exist('PatternPhaseCalculations', 'var') == 0)
        PatternPhaseCalculations = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternPhaseCalculations', PatternPhaseCalculations);
    
    if (exist('RadialMusicParameters', 'var') == 0)
        RadialMusicParameters = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RadialMusicParameters', RadialMusicParameters);
    
    if (exist('MergedCount', 'var') == 0)
        MergedCount = '';
    end
    netcdf.putAtt(ncid, varid_global, 'MergedCount', MergedCount);
    
    if (exist('RadialMinimumMergePoints', 'var') == 0)
        RadialMinimumMergePoints = '';
    end
    netcdf.putAtt(ncid, varid_global, 'RadialMinimumMergePoints', RadialMinimumMergePoints);
    
    if (exist('FirstOrderCalc', 'var') == 0)
        FirstOrderCalc = '';
    end
    netcdf.putAtt(ncid, varid_global, 'FirstOrderCalc', FirstOrderCalc);
    
    if (exist('MergeMethod', 'var') == 0)
        MergeMethod = '';
    end
    netcdf.putAtt(ncid, varid_global, 'MergeMethod', MergeMethod);
    
    if (exist('PatternMethod', 'var') == 0)
        PatternMethod = '';
    end
    netcdf.putAtt(ncid, varid_global, 'PatternMethod', PatternMethod);
    
    if (exist('TransmitSweepRateHz', 'var') == 0)
        TransmitSweepRateHz = '';
    end
    netcdf.putAtt(ncid, varid_global, 'TransmitSweepRateHz', TransmitSweepRateHz);
    
    if (exist('TransmitBandwidthKHz', 'var') == 0)
        TransmitBandwidthKHz = '';
    end
    netcdf.putAtt(ncid, varid_global, 'TransmitBandwidthKHz', TransmitBandwidthKHz);
    
    if (exist('SpectraRangeCells', 'var') == 0)
        SpectraRangeCells = '';
    end
    netcdf.putAtt(ncid, varid_global, 'SpectraRangeCells', SpectraRangeCells);
    
    if (exist('SpectraDopplerCells', 'var') == 0)
        SpectraDopplerCells = '';
    end
    netcdf.putAtt(ncid, varid_global, 'SpectraDopplerCells', SpectraDopplerCells);
    
    %%
    
    
    %% Exit definition mode
    netcdf.endDef(ncid);
    
    %%
    
    %% Write variable data
    netcdf.putVar(ncid, varid_bearing, bearing_dim);
    netcdf.putVar(ncid, varid_range, range_dim);
    netcdf.putVar(ncid, varid_t, 0, R.time - datenum(1950,1,1));
    netcdf.putVar(ncid, varid_depth, 0);
    netcdf.putVar(ncid, varid_lat, latd);
    netcdf.putVar(ncid, varid_lon, lond);
    netcdf.putVar(ncid, varid_crs, 0);
    netcdf.putVar(ncid, varid_sdncruise, site_code);
    netcdf.putVar(ncid, varid_sdnstation, platform_code);
    netcdf.putVar(ncid, varid_sdnedmocode, EDMO_code);
    netcdf.putVar(ncid, varid_sdnlocalcdiid, id);
    netcdf.putVar(ncid, varid_sdnreferences, TDS_catalog);
    netcdf.putVar(ncid, varid_sdnxlink, xlink');
    netcdf.putVar(ncid, varid_speed, velo);
    netcdf.putVar(ncid, varid_direction, head);
    netcdf.putVar(ncid, varid_u, velu);
    netcdf.putVar(ncid, varid_v, velv);
    netcdf.putVar(ncid, varid_espc, espc);
    netcdf.putVar(ncid, varid_etmp, etmp);
    netcdf.putVar(ncid, varid_maxv, maxv);
    netcdf.putVar(ncid, varid_minv, minv);
    netcdf.putVar(ncid, varid_ersc, ersc);
    netcdf.putVar(ncid, varid_ertc, ertc);
    netcdf.putVar(ncid, varid_xdst, xdst);
    netcdf.putVar(ncid, varid_ydst, ydst);
    netcdf.putVar(ncid, varid_sprc, sprc);
    netcdf.putVar(ncid, varid_narx, numSites);
    netcdf.putVar(ncid, varid_natx, numSites);
    netcdf.putVar(ncid, varid_sltr, round((siteLat-addOffset)./scaleFactor));
    netcdf.putVar(ncid, varid_slnr, round((siteLon-addOffset)./scaleFactor));
    netcdf.putVar(ncid, varid_sltt, round((siteLat-addOffset)./scaleFactor));
    netcdf.putVar(ncid, varid_slnt, round((siteLon-addOffset)./scaleFactor));
    netcdf.putVar(ncid, varid_scdr, siteCodeArr');
    netcdf.putVar(ncid, varid_scdt, siteCodeArr');
    netcdf.putVar(ncid, varid_tqc, int8(sdnTime_QCflag));
    netcdf.putVar(ncid, varid_posqc, sdnPosition_QCflag);
    netcdf.putVar(ncid, varid_dqc, int8(sdnDepth_QCflag));
    netcdf.putVar(ncid, varid_ovqc, overall_QCflag);
    netcdf.putVar(ncid, varid_owtr, overWater_QCflag);
    netcdf.putVar(ncid, varid_mdfl, medianFilter_QCflag);
    netcdf.putVar(ncid, varid_vart, temporalDerivativeThreshold_QCflag);
    netcdf.putVar(ncid, varid_cspd, velocityThreshold_QCflag);
    netcdf.putVar(ncid, varid_avrb, averageRadialBearing_QC_flag);
    netcdf.putVar(ncid, varid_rdct, radialCount_QC_flag);
    
    %%
    
    %% Close file
    netcdf.close(ncid);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end
%%

%% Retrieve information about the nc file
try
    ncfileInfo = dir(ncfile);
    ncFilesize = ncfileInfo.bytes/1024;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

if(R2C_err==0)
    disp(['[' datestr(now) '] - - ' 'mat2netCDF_v33.m successfully executed.']);
end

return
