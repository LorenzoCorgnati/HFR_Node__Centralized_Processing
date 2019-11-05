%% remapMetNoNCradial2ruv.m
% This function remaps the MetNo radial netCDF files to the ruv structure.

% INPUT:
%         ncRad: OpenDAP url to radial netCDF file.
%         time_idx: time index of the hourly radial to be remapped within 
%                   the aggregated netCDF file.

% OUTPUT:
%         rMNn2r_err: error flag (0 = correct, 1 = error)
%         REMAP: radial in ruv structure.


% Author: Lorenzo Corgnati
% Date: November 4, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [rMNn2r_err,REMAP] = remapMetNoNCradial2ruv(ncRad,time_idx)

disp(['[' datestr(now) '] - - ' 'remapMetNoNCradial2ruv.m started.']);

rMNn2r_err = 0;

warning('off', 'all');

%% Load netCDF files

try
    % Coordinate variables
    NC.time = ncread(ncRad,'time',time_idx,1);
    NC.bearing = ncread(ncRad,'bearing');
    NC.range = ncread(ncRad,'range');
    
    NC.Lat = ncread(ncRad,'lat');
    NC.Lon = ncread(ncRad,'lon');
    
    % Data variables
    NC.direction = ncread(ncRad,'direction',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.velo = ncread(ncRad,'velo',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.U = ncread(ncRad,'u',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.V = ncread(ncRad,'v',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.espc = ncread(ncRad,'espc',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.etmp = ncread(ncRad,'etmp',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.maxv = ncread(ncRad,'maxv',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.minv = ncread(ncRad,'minv',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.ersc = ncread(ncRad,'ersc',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.ertc = ncread(ncRad,'ertc',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.xdst = ncread(ncRad,'xdst',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.ydst = ncread(ncRad,'ydst',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.sprc = ncread(ncRad,'sprc',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
    NC.x = ncread(ncRad,'x');
    NC.y = ncread(ncRad,'y');
    NC.vflg = ncread(ncRad,'vflg',[1,1,time_idx],[length(NC.range),length(NC.bearing),1]);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rMNn2r_err = 1;
    return
end

%%

%% Remap gridded data to arrays according to HFR_Progs data structure

try
    % Initialize REMAP structure
    REMAP = RADIALstruct(1);
    REMAP.Type = ['RDL' ncreadatt(ncRad,'/','PatternType')];
    REMAP.SiteName = ncreadatt(ncRad,'/','site');
    if((strcmp(REMAP.SiteName,'BERL')) || (strcmp(REMAP.SiteName,'JOMF')))
        REMAP.SiteCode = 1;
    else
        REMAP.SiteCode = 2;
    end
    originCell = strsplit(ncreadatt(ncRad,'/','Origin'));
    REMAP.SiteOrigin = [str2double(originCell{2}) str2double(originCell{1})];
    REMAP.TimeStamp = (double(NC.time)/86400) + datenum('01-01-1970');
    REMAP.TimeZone = 'GMT+0.000';
    
    % Retrieve ruv header information from netCDF global attributes
    REMAP.OtherMetadata.Header{1,1} = ['%CTF: ' ncreadatt(ncRad,'/','CTF')];
    REMAP.OtherMetadata.Header{2,1} = ['%FileType: ' ncreadatt(ncRad,'/','filetype')];
    REMAP.OtherMetadata.Header{3,1} = ['%LLUVSpec: ' ncreadatt(ncRad,'/','LLUVSpec')];
    REMAP.OtherMetadata.Header{4,1} = ['%UUID: ' ncreadatt(ncRad,'/','UUID')];
    REMAP.OtherMetadata.Header{5,1} = ['Manufacturer: ' ncreadatt(ncRad,'/','manufacturer')];
    REMAP.OtherMetadata.Header{6,1} = ['%Site: ' ncreadatt(ncRad,'/','site')];
    tsVec = datevec(REMAP.TimeStamp);
    REMAP.OtherMetadata.Header{7,1} = ['%TimeStamp: ' num2str(tsVec(1),'%04d') ' ' num2str(tsVec(2),'%02d') ' ' num2str(tsVec(3),'%02d') '  ' num2str(tsVec(4),'%02d') ' ' num2str(tsVec(5),'%02d') ' ' num2str(tsVec(6),'%02d')];
    REMAP.OtherMetadata.Header{8,1} = ['%TimeZone: ' ncreadatt(ncRad,'/','TimeZone')];
    REMAP.OtherMetadata.Header{9,1} = ['%TimeCoverage: ' ncreadatt(ncRad,'/','TimeCoverage')];
    REMAP.OtherMetadata.Header{10,1} = ['%Origin:  ' ncreadatt(ncRad,'/','Origin')];
    REMAP.OtherMetadata.Header{11,1} = ['%GreatCircle: ' ncreadatt(ncRad,'/','GreatCircle')];
    REMAP.OtherMetadata.Header{12,1} = ['%GeodVersion: ' ncreadatt(ncRad,'/','GeodVersion')];
    REMAP.OtherMetadata.Header{13,1} = ['%LLUVTrustData: ' ncreadatt(ncRad,'/','LLUVTrustData')];
    REMAP.OtherMetadata.Header{14,1} = ['%RangeStart: ' ncreadatt(ncRad,'/','RangeStart')];
    REMAP.OtherMetadata.Header{15,1} = ['%RangeEnd: ' ncreadatt(ncRad,'/','RangeEnd')];
    REMAP.OtherMetadata.Header{16,1} = ['%RangeResolutionKMeters: ' ncreadatt(ncRad,'/','RangeResolutionKMeters')];
    REMAP.OtherMetadata.Header{17,1} = ['%AntennaBearing: ' ncreadatt(ncRad,'/','AntennaBearing')];
    REMAP.OtherMetadata.Header{18,1} = ['%ReferenceBearing: ' ncreadatt(ncRad,'/','ReferenceBearing')];
    REMAP.OtherMetadata.Header{19,1} = ['%AngularResolution: ' ncreadatt(ncRad,'/','AngularResolution')];
    REMAP.OtherMetadata.Header{20,1} = ['%SpatialResolution: ' ncreadatt(ncRad,'/','SpatialResolution')];
    REMAP.OtherMetadata.Header{21,1} = ['%PatternType: ' ncreadatt(ncRad,'/','PatternType')];
    REMAP.OtherMetadata.Header{22,1} = ['%PatternDate: ' ncreadatt(ncRad,'/','PatternDate')];
    REMAP.OtherMetadata.Header{23,1} = ['%PatternResolution: ' ncreadatt(ncRad,'/','PatternResolution')];
    REMAP.OtherMetadata.Header{24,1} = ['%TransmitCenterFreqMHz: ' ncreadatt(ncRad,'/','TransmitCenterFreqMHz')];
    REMAP.OtherMetadata.Header{25,1} = ['%DopplerResolutionHzPerBin: ' ncreadatt(ncRad,'/','DopplerResolutionHzPerBin')];
    REMAP.OtherMetadata.Header{26,1} = ['%FirstOrderMethod: ' ncreadatt(ncRad,'/','FirstOrderMethod')];
    REMAP.OtherMetadata.Header{27,1} = ['%BraggSmoothingPoints: ' ncreadatt(ncRad,'/','BraggSmoothingPoints')];
    REMAP.OtherMetadata.Header{28,1} = ['%BraggHasSecondOrder: ' ncreadatt(ncRad,'/','BraggHasSecondOrder')];
    REMAP.OtherMetadata.Header{29,1} = ['%RadialBraggPeakDropOff: ' ncreadatt(ncRad,'/','RadialBraggPeakDropOff')];
    REMAP.OtherMetadata.Header{30,1} = ['%RadialBraggPeakNull: ' ncreadatt(ncRad,'/','RadialBraggPeakNull')];
    REMAP.OtherMetadata.Header{31,1} = ['%RadialBraggNoiseThreshold: ' ncreadatt(ncRad,'/','RadialBraggNoiseThreshold')];
    REMAP.OtherMetadata.Header{32,1} = ['%PatternAmplitudeCorrections: ' ncreadatt(ncRad,'/','PatternAmplitudeCorrections')];
    REMAP.OtherMetadata.Header{33,1} = ['%PatternPhaseCorrections: ' ncreadatt(ncRad,'/','PatternPhaseCorrections')];
    REMAP.OtherMetadata.Header{34,1} = ['%PatternAmplitudeCalculations: ' ncreadatt(ncRad,'/','PatternAmplitudeCalculations')];
    REMAP.OtherMetadata.Header{35,1} = ['%PatternPhaseCalculations: ' ncreadatt(ncRad,'/','PatternPhaseCalculations')];
    REMAP.OtherMetadata.Header{36,1} = ['%RadialMusicParameters: ' ncreadatt(ncRad,'/','RadialMusicParameters')];
    REMAP.OtherMetadata.Header{37,1} = ['%MergedCount: ' ncreadatt(ncRad,'/','MergedCount')];
    REMAP.OtherMetadata.Header{38,1} = ['%RadialMinimumMergePoints: ' ncreadatt(ncRad,'/','RadialMinimumMergePoint')];
    REMAP.OtherMetadata.Header{39,1} = ['%FirstOrderCalc: ' ncreadatt(ncRad,'/','FirstOrderCalc')];
    REMAP.OtherMetadata.Header{40,1} = ['%MergeMethod: ' ncreadatt(ncRad,'/','MergeMethod')];
    REMAP.OtherMetadata.Header{41,1} = ['%PatternMethod: ' ncreadatt(ncRad,'/','PatternMethod')];
    REMAP.OtherMetadata.Header{42,1} = ['%TransmitSweepRateHz: ' ncreadatt(ncRad,'/','TransmitSweepRateHz')];
    REMAP.OtherMetadata.Header{43,1} = ['%TransmitBandwidthKHz: ' ncreadatt(ncRad,'/','TransmitBandwidthKHz')];
    REMAP.OtherMetadata.Header{44,1} = ['%SpectraRangeCells: ' ncreadatt(ncRad,'/','RangeCells')];
    REMAP.OtherMetadata.Header{45,1} = ['%SpectraDopplerCells: ' ncreadatt(ncRad,'/','DopplerCells')];
    REMAP.OtherMetadata.Header{46,1} = ['%TableType: ' ncreadatt(ncRad,'/','TableType')];
    REMAP.OtherMetadata.Header{47,1} = ['%TableColumns: ' ncreadatt(ncRad,'/','TableColumns')];
    REMAP.OtherMetadata.Header{48,1} = ['%TableColumnTypes: ' ncreadatt(ncRad,'/','TableColumnsTypes')];
    REMAP.OtherMetadata.Header{49,1} = ['%TableRows: ' ncreadatt(ncRad,'/','TableRows')];
    
    % Fill REMAP structure
    remap_idx = 0;
    currentU = NC.U;
    arraySize = size(currentU(~isnan(currentU)),1);
    REMAP.LonLat = NaN*ones(arraySize,2);
    REMAP.RangeBearHead = NaN*ones(arraySize,3);
    REMAP.RadComp = NaN*ones(arraySize,1);
    REMAP.U = NaN*ones(arraySize,1);
    REMAP.V = NaN*ones(arraySize,1);
    REMAP.Flag = ones(arraySize,1);
    
    % Scan latitude and longitude grids, retrieve bearing and range and write
    % to the remapped structure
    for i=1:size(currentU,1)
        for j=1:size(currentU,2)
            if(~isnan(currentU(i,j)))
                remap_idx = remap_idx + 1;
                [distance,a12,a21] = m_idist(REMAP.SiteOrigin(1),REMAP.SiteOrigin(2),NC.Lon(i,j),NC.Lat(i,j),'wgs84');
                [val,idx_b]=min(abs(NC.bearing-a12));
                bearing = mod((90-NC.bearing(idx_b)),360);
                bearingRD(remap_idx,1) = NC.bearing(idx_b);
                [val,idx_r]=min(abs(NC.range-(distance/1000)));
                range = NC.range(idx_r);
                REMAP.LonLat(remap_idx,1) = NC.Lon(i,j);
                REMAP.LonLat(remap_idx,2) = NC.Lat(i,j);
                REMAP.RangeBearHead(remap_idx,1) = range;
                REMAP.RangeBearHead(remap_idx,2) = bearing;
                REMAP.RangeBearHead(remap_idx,3) = wrapTo360(bearing+180);
                REMAP.U(remap_idx,1) = NC.U(i,j);
                REMAP.V(remap_idx,1) = NC.V(i,j);
                REMAP.RadComp(remap_idx,1) = -NC.velo(i,j);
                REMAP.Error(remap_idx,1) = NC.etmp(i,j);
                espc(remap_idx,1) = NC.espc(i,j);
                etmp(remap_idx,1) = NC.etmp(i,j);
                maxv(remap_idx,1) = NC.maxv(i,j);
                minv(remap_idx,1) = NC.minv(i,j);
                ersc(remap_idx,1) = NC.ersc(i,j);
                ertc(remap_idx,1) = NC.ertc(i,j);
                xdst(remap_idx,1) = NC.xdst(i,j);
                ydst(remap_idx,1) = NC.ydst(i,j);
                rdva(remap_idx,1) = NC.velo(i,j);
                sprc(remap_idx,1) = NC.sprc(i,j);
            end
        end
    end
    
    % Write raw data
    REMAP.OtherMetadata.RawData = [REMAP.LonLat REMAP.U REMAP.V REMAP.Flag espc etmp maxv minv ersc ertc xdst ydst REMAP.RangeBearHead(:,1) bearingRD rdva wrapTo360(bearingRD+180) sprc];
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rMNn2r_err = 1;
    return
end

%%

if(rMNn2r_err==0)
    disp(['[' datestr(now) '] - - ' 'remapMetNoNCradial2ruv.m successfully executed.']);
end

return