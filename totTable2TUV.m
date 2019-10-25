%% totTable2TUV.m
% This function fills the TUV structure with data contained within the table
% loaded from input total file.

% INPUT:
%         table: table loaded from input total file.
%         fields: field names of the table loaded from input total file.
%         ts: timestamp of the total file
%         siteLon: array containing longitudes of the contributing radar
%                  sites
%         siteLat: array containing latitudes of the contributing radar
%                  sites

% OUTPUT:
%         t2T_err: error flag (0 = correct, 1 = error)
%         tuv: TUV structure containing total data


% Author: Lorenzo Corgnati
% Date: June 12, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [t2T_err,tuv] = totTable2TUV(table,fields,ts,siteLon,siteLat)

disp(['[' datestr(now) '] - - ' 'totTable2TUV.m started.']);

t2T_err = 0;

warning('off', 'all');

%% Create the TUV structure

try
    tuv = TUVstruct([size(table,1),1], 1);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    t2T_err = 1;
end

%%

%% Fill the TUV structure with the total data

try
    % Add the TimeStamp
    tsCell = strsplit(ts);
    tsVec = [];
    for tsCell_idx=1:length(tsCell)
        tsVec = [tsVec str2double(tsCell{tsCell_idx})];
    end
    tuv.TimeStamp = datenum(tsVec);
    
    % Add longitude and latitude
    IndexC = strfind(fields, 'LOND');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.LonLat(:,1) = table(:,Index);
    IndexC = strfind(fields, 'LATD');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.LonLat(:,2) = table(:,Index);
    
    % Add depth
    tuv.Depth = zeros(size(table,1),1);
    
    % Add U and V component of velocity
    IndexC = strfind(fields, 'VELU');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.U = table(:,Index);
    IndexC = strfind(fields, 'VELV');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.V = table(:,Index);
    
    % Add standard deviation of U and V component of velocity
    IndexC = strfind(fields, 'UQAL');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.ErrorEstimates(1,1).Uerr = table(:,Index);
    IndexC = strfind(fields, 'VQAL');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.ErrorEstimates(1,1).Verr = table(:,Index);
    
    % Add covariance
    IndexC = strfind(fields, 'CQAL');
    Index = find(not(cellfun('isempty', IndexC)));
    tuv.ErrorEstimates(1,1).UVCovariance = table(:,Index);
    
    % Retrieve the number of contributing radial sites
    SnCNindices = ~cellfun(@isempty,regexp(fields,'S[0-255]+CN'));
    numRadSite = sum(SnCNindices);
    
    % Add number of radial vectors which contributed to the total vector
    contSum = zeros(size(table,1),1);
    for cont_idx=1:numRadSite
        IndexC = strfind(fields, ['S' num2str(cont_idx) 'CN']);
        Index = find(not(cellfun('isempty', IndexC)));  
        contSum = contSum + table(:,Index);
    end
    tuv.OtherMatrixVars.makeTotals_TotalsNumRads = contSum;
    
    % Evaluate GDOP and add it
    % Evaluate radial angles for each grid cell
    for site_idx=1:length(siteLon)
        [s,radialAngles(:,site_idx),a21] = m_idist(siteLon(site_idx), siteLat(site_idx), tuv.LonLat(:,1), tuv.LonLat(:,2), 'wgs84');
    end
    % Evaluate GDOP for each grid cell (as the square root of the trace of the covariance matrix of radial angles)
    for cell_idx=1:size(table,1)
        tuv.ErrorEstimates(1,1).TotalErrors(cell_idx) = sqrt(trace(gdop_max_orthog(radialAngles(cell_idx,:))));
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    t2T_err = 1;
end

%%

if(t2T_err==0)
    disp(['[' datestr(now) '] - - ' 'totTable2TUV.m successfully executed.']);
end

return