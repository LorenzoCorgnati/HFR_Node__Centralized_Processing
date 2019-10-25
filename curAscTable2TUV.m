%% curAscTable2TUV.m
% This function fills the TUV structure with data contained within the table
% loaded from input cur_asc total file.

% INPUT:
%         table: table loaded from input total file.
%         fields: field names of the table loaded from input total file.
%         ts: timestamp of the total file
%         gridLon: matrix containing the longitudes of the regular grid
%         gridLat: matrix containing the latitudes of the regular grid
%         siteLon: array containing the longitudes of the radial sites
%         siteLat: array containing the latitudes of the radial sites

% OUTPUT:
%         ca2T_err: error flag (0 = correct, 1 = error)
%         tuv: TUV structure containing total data


% Author: Lorenzo Corgnati
% Date: October 14, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [ca2T_err,tuv] = curAscTable2TUV(table,fields,ts,gridLon,gridLat,siteLon,siteLat)

disp(['[' datestr(now) '] - - ' 'curAscTable2TUV.m started.']);

ca2T_err = 0;

warning('off', 'all');


%% Create the TUV structure

try
    tuv = TUVstruct([size(table,1),1], 1);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    ca2T_err = 1;
end

%%

%% Fill the TUV structure with the total data

try
    % Create the arrays of longitude and latitude
    lonG = gridLon(1,:);
    latG = gridLat(:,1);
    
    % Find the index of the X indices field
    XindIndexC = strfind(fields, 'IX');
    XindIndex = find(not(cellfun('isempty', XindIndexC)));
    % Find the index of the Y indices field
    YindIndexC = strfind(fields, 'IY');
    YindIndex = find(not(cellfun('isempty', YindIndexC)));
    
    % Find the index of the U accuracy field
    uAccIndexC = strfind(fields, 'Acc_U[m/s]');
    uAccIndex = find(not(cellfun('isempty', uAccIndexC)));
    % Find the index of the V accuracy field
    vAccIndexC = strfind(fields, 'Acc_V[m/s]');
    vAccIndex = find(not(cellfun('isempty', vAccIndexC)));
    
    % Find the index of the U field
    uIndexC = strfind(fields, 'U[m/s]');
    uIndex = find(not(cellfun('isempty', uIndexC)));
    uIndex = uIndex(find(uIndex~=uAccIndex));
    % Find the index of the V field
    vIndexC = strfind(fields, 'V[m/s]');
    vIndex = find(not(cellfun('isempty', vIndexC)));
    vIndex = vIndex(find(vIndex~=vAccIndex));
    
    % Scan the data table and fill the TUV structure
    for tab_idx=1:size(table,1)
        % Map the X and Y indices to lat/lon coordinates
        tuv.LonLat(tab_idx,1) = lonG(table(tab_idx,XindIndex));
        tuv.LonLat(tab_idx,2) = latG(table(tab_idx,YindIndex));
    end
    
    % Insert the U and V data
    tuv.U = table(:,uIndex);
    tuv.V = table(:,vIndex);
    
    % Insert the U accuracy and V accuracy data
    tuv.ErrorEstimates(1,1).Uerr = table(:,uAccIndex);
    tuv.ErrorEstimates(1,1).Verr = table(:,vAccIndex);
    
    % Evaluate GDOP and add it
    % Evaluate radial angles for each grid cell
    for site_idx=1:length(siteLon)
        [s,radialAngles(:,site_idx),a21] = m_idist(siteLon(site_idx), siteLat(site_idx), tuv.LonLat(:,1), tuv.LonLat(:,2), 'wgs84');
    end
    % Evaluate GDOP for each grid cell (as the square root of the trace of the covariance matrix of radial angles)
    for cell_idx=1:size(table,1)
        tuv.ErrorEstimates(1,1).TotalErrors(cell_idx) = sqrt(trace(gdop_max_orthog(radialAngles(cell_idx,:))));
    end
    
    % Add the TimeStamp
    tsCell = strsplit(ts);
    tsVec = [];
    for tsCell_idx=1:length(tsCell)
        tsVec = [tsVec str2double(tsCell{tsCell_idx})];
    end
    tuv.TimeStamp = datenum(tsVec);
    
    % Add depth
    tuv.Depth = zeros(size(table,1),1);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    ca2T_err = 1;
end

%%

if(ca2T_err==0)
    disp(['[' datestr(now) '] - - ' 'curAscTable2TUV.m successfully executed.']);
end

return