%% gridTUV.m
% This function places total current data on a regular grid created
% according to the given geographical boundaries and the grid resolution.

% INPUT:
%         TUV: input TUV structure on an irregular grid.
%         networkData: cell array containing information about the network
%                      (metadata).
%         networkFields: field names of the cell array containing
%                       information about the network.

% OUTPUT:
%         gT_err: error flag (0 = correct, 1 = error).
%         TUV: TUV structure on a regular rectangular grid.


% Author: Lorenzo Corgnati
% Date: June 21, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [gT_err, TUV] = gridTUV(TUV,networkData,networkFields)

disp(['[' datestr(now) '] - - ' 'gridTUV.m started.']);

gT_err = 0;

warning('off', 'all');

%% Create a regular LonLat grid given the geographical boundaries and the grid resolution

try
    % Find the indices of the geospatial boundary fields
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    % Find the index of the grid resolution field
    grid_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'grid_resolution'))));
    
    [TUV.gridLon, TUV.gridLat] = LonLat_grid([networkData{geospatial_lon_minIndex},networkData{geospatial_lat_minIndex}], [networkData{geospatial_lon_maxIndex},networkData{geospatial_lat_maxIndex}], networkData{grid_resolutionIndex}, 'km');
    TUV.gridLat = flipud(TUV.gridLat);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    gT_err = 1;
end

%%

%% Fill the regular grid with the total current data

try
    % Prepare variables
    TUV.U_grid = NaN.*ones(size(TUV.gridLat));
    TUV.V_grid = NaN.*ones(size(TUV.gridLat));
    
    TUV.U_std = NaN.*ones(size(TUV.gridLat));
    TUV.V_std = NaN.*ones(size(TUV.gridLat));
    
    TUV.covariance = NaN.*ones(size(TUV.gridLat));
    
    TUV.GDOP = NaN.*ones(size(TUV.gridLat));
    
    TUV.DDENS = NaN.*ones(size(TUV.gridLat));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    gT_err = 1;
end

% Populate variables
try
    % Scan the regular lat/lon grid
    for latGrid_idx=1:size(TUV.gridLat,1)
        for lonGrid_idx=1:size(TUV.gridLat,2)
            % Set the cumulative variables
            Uvec = [];
            Vvec = [];
            stdUvec = [];
            stdVvec = [];
            COVvec = [];
            GDOPvec = [];
            DDENSvec = [];
            % Retrieve lower left and upper right corners of the current cell of the regular grid
            [URlon,URlat]=km2lonlat(TUV.gridLon(latGrid_idx,lonGrid_idx),TUV.gridLat(latGrid_idx,lonGrid_idx),networkData{grid_resolutionIndex}/2,networkData{grid_resolutionIndex}/2);
            [LLlon,LLlat]=km2lonlat(TUV.gridLon(latGrid_idx,lonGrid_idx),TUV.gridLat(latGrid_idx,lonGrid_idx),-networkData{grid_resolutionIndex}/2,-networkData{grid_resolutionIndex}/2);
            % Scan the LonLat matrix and find the positions included in each cell of the regular grid
            for vecLonLat_idx=1:size(TUV.LonLat,1)
                if((TUV.LonLat(vecLonLat_idx,1)>=LLlon) & (TUV.LonLat(vecLonLat_idx,1)<=URlon) & (TUV.LonLat(vecLonLat_idx,2)>=LLlat) & (TUV.LonLat(vecLonLat_idx,2)<=URlat))
                    Uvec = [Uvec TUV.U(vecLonLat_idx)];
                    Vvec = [Vvec TUV.V(vecLonLat_idx)];
                    stdUvec = [stdUvec TUV.ErrorEstimates(1,1).Uerr(vecLonLat_idx)];
                    stdVvec = [stdVvec TUV.ErrorEstimates(1,1).Verr(vecLonLat_idx)];
                    COVvec = [COVvec TUV.ErrorEstimates(1,1).UVCovariance(vecLonLat_idx)];
                    GDOPvec = [GDOPvec TUV.ErrorEstimates(1,1).TotalErrors(vecLonLat_idx)];
                    DDENSvec = [DDENSvec TUV.OtherMatrixVars.makeTotals_TotalsNumRads(vecLonLat_idx)];
                end
            end
            % Evaluate the average values to be inserted in the gridded variables
            % U and V components of current velocity
            TUV.U_grid(latGrid_idx,lonGrid_idx,1) = mean(Uvec)*0.01;
            TUV.V_grid(latGrid_idx,lonGrid_idx,1) = mean(Vvec)*0.01;
            % U and V standard errors
            TUV.U_std(latGrid_idx,lonGrid_idx,1) = mean(stdUvec)*0.01;
            TUV.V_std(latGrid_idx,lonGrid_idx,1) = mean(stdVvec)*0.01;
            % UV covariance
            TUV.covariance(latGrid_idx,lonGrid_idx,1) = mean(COVvec)*0.0001;
            % GDOP
            TUV.GDOP(latGrid_idx,lonGrid_idx,1) = mean(GDOPvec);
            % DDENS
            TUV.DDENS(latGrid_idx,lonGrid_idx,1) = mean(DDENSvec);
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    gT_err = 1;
end

if(gT_err==0)
    disp(['[' datestr(now) '] - - ' 'gridTUV.m successfully executed.']);
end

return