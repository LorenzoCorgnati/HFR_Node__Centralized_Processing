%% cradAsciiTable2RAD.m
% This function fills the radial structure with data contained within the table
% loaded from input crad_ascii radial file.

% INPUT:
%         table: table loaded from input radial file.
%         fields: field names of the table loaded from input radial file.
%         ts: timestamp of the total file
%         gridLon: matrix containing the longitudes of the regular grid
%         gridLat: matrix containing the latitudes of the regular grid
%         siteLon: longitude of the radial site
%         siteLat: latitude of the radial site

% OUTPUT:
%         ca2T_err: error flag (0 = correct, 1 = error)
%         rad: radial structure containing radial data


% Author: Lorenzo Corgnati
% Date: May 14, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [ca2T_err,rad] = cradAsciiTable2RAD(table,fields,ts,gridLon,gridLat,siteLon,siteLat)

disp(['[' datestr(now) '] - - ' 'cradAsciiTable2RAD.m started.']);

ca2T_err = 0;

warning('off', 'all');

%%

%% Create the radial structure 

try
    % Add longitude and latitude
    rad.Lon = gridLon; %gridLon(1,:);
    rad.Lat = gridLat; %gridLat(:,1);
    
    % Prepare variables
    rad.rdva = NaN.*ones(length(gridLat),length(gridLon));
    rad.drva = NaN.*ones(length(gridLat),length(gridLon));
    rad.ewct = NaN.*ones(length(gridLat),length(gridLon));
    rad.nsct = NaN.*ones(length(gridLat),length(gridLon));
    rad.hcss = NaN.*ones(length(gridLat),length(gridLon));
    rad.eacc = NaN.*ones(length(gridLat),length(gridLon));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    ca2T_err = 1;
end

%%

%% Fill the radial structure with the radial data

try
    % Find the index of the LatC field
    latcIndexC = strfind(fields, 'LatC');
    latcIndex = find(not(cellfun('isempty', latcIndexC)));
    % Find the index of the LonC field
    loncIndexC = strfind(fields, 'LonC');
    loncIndex = find(not(cellfun('isempty', loncIndexC)));
    
    % Find the index of the KUR field
    kurIndexC = strfind(fields, 'KUR');
    kurIndex = find(not(cellfun('isempty', kurIndexC)));
    
    % Find the index of the SNV field
    snvIndexC = strfind(fields, 'SNV');
    snvIndex = find(not(cellfun('isempty', snvIndexC)));
    
    % Find the index of the ACC field
    snsIndexC = strfind(fields, 'SNS');
    snsIndex = find(not(cellfun('isempty', snsIndexC)));

    % Find the index of the SNR field
    snrIndexC = strfind(fields, 'SNR');
    snrIndex = find(not(cellfun('isempty', snrIndexC)));
    
    % Find the index of the PWR field
    pwrIndexC = strfind(fields, 'PWR');
    pwrIndex = find(not(cellfun('isempty', pwrIndexC)));
    
    % Scan the grid and fill the radial structure with table data
    for grd_y=1:length(gridLat)
        for grd_x=1:length(gridLon)
            % Retrieve table line number from Jan's formula: [x,y] -> (x-1)*100 + (y-1) + 1
%             lineNumber = (grd_x-1)*100 + (grd_y-1) + 1;
            lineNumber = (grd_x-1)*length(gridLat) + (grd_y-1) + 1;
            
%             if((table(lineNumber,latcIndex)~=-999) && (table(lineNumber,loncIndex)~=-999))
            if(table(lineNumber,kurIndex)~=0)
                
                % Insert velocity magnitude data
                rad.rdva(grd_y,grd_x) = table(lineNumber,snvIndex) / table(lineNumber,snrIndex);
                
                % Evaluate velocity bearing
                [s,rad.drva(grd_y,grd_x),a21] = m_idist(siteLon, siteLat, rad.Lon(grd_x), rad.Lat(grd_y), 'wgs84');
                
                % Evaluate U and V components of the velocity
                rad.ewct(grd_y,grd_x) = rad.rdva(grd_y,grd_x) * sin(rad.drva(grd_y,grd_x) * (pi/180));
                rad.nsct(grd_y,grd_x) = rad.rdva(grd_y,grd_x) * cos(rad.drva(grd_y,grd_x) * (pi/180));
                
                % Insert variance data
                rad.hcss(grd_y,grd_x) = table(lineNumber,snsIndex) / table(lineNumber,snrIndex);
                
                % Insert accuracy data
                rad.eacc(grd_y,grd_x) = rad.hcss(grd_y,grd_x) / sqrt(table(lineNumber,kurIndex));
            end
        end
    end
       
    % Add the TimeStamp
    tsCell = strsplit(ts);
    tsVec = [];
    for tsCell_idx=1:length(tsCell)
        tsVec = [tsVec str2double(tsCell{tsCell_idx})];
    end
    rad.TimeStamp = datenum(tsVec);
    
    % Add depth
    rad.depth = 0;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    ca2T_err = 1;
end

%%

if(ca2T_err==0)
    disp(['[' datestr(now) '] - - ' 'cradAsciiTable2RAD.m successfully executed.']);
end

return