%% curAscGridSpec.m
% This function retrieves top-left point coordinates of the first gridcell,
% the cell size and the number of lon and lat gridcells from the header of
% WERA cur_asc total files.

% INPUT:
%         header: header of the WERA cur_asc file.

% OUTPUT:
%         sCC_err: error flag (0 = correct, 1 = error)
%         tLlat: latitude of the top-left point of the first gridcell
%         tLlon: longitude of the top-left point of the first gridcell
%         cellSize: cell size in km
%         lonGrid: number of lon gridcells
%         latGrid: number of lat gridcells

% Author: Lorenzo Corgnati
% Date: October 5, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [sCC_err, tLlat, tLlon, cellSize, lonGrid,latGrid] = curAscGridSpec(header)

disp(['[' datestr(now) '] - - ' 'curAscGridSpec.m started.']);

sCC_err = 0;

warning('off', 'all');

try
    % Read the header
    for line_idx=1:length(header)
        splitLine = regexp(header{line_idx}, '[ \t]+', 'split');
        % Find the column names and indices
        if(length(splitLine)>1)
            if((any(strcmp(splitLine,'LAT(1,1)'))) && (any(strcmp(splitLine,'LON(1,1)'))) && (any(strcmp(splitLine,'DGT[km]'))) && (any(strcmp(splitLine,'NX'))) && (any(strcmp(splitLine,'NY'))))
                latLogical = strcmp(splitLine,'LAT(1,1)');
                lat_idx = find(latLogical, 1, 'first');
                lonLogical = strcmp(splitLine,'LON(1,1)');
                lon_idx = find(lonLogical, 1, 'first');
                dgtLogical = strcmp(splitLine,'DGT[km]');
                dgt_idx = find(dgtLogical, 1, 'first');
                nxLogical = strcmp(splitLine,'NX');
                nx_idx = find(nxLogical, 1, 'first');
                nyLogical = strcmp(splitLine,'NY');
                ny_idx = find(nyLogical, 1, 'first');
                % Find the data
                data_idx = line_idx + 1;
                expressionNum = '([0-9]+)';
                [startIndexNum,endIndexNum] = regexp(header{data_idx},expressionNum);
                while(isempty(startIndexNum))
                    data_idx = data_idx + 1;
                    [startIndexNum,endIndexNum] = regexp(header{data_idx},expressionNum);
                end
                splitLine = regexp(header{data_idx}, '[ \t]+', 'split');
                break;
            end
        end
    end
    % Fill the output variables
    tLlat = str2double(splitLine{lat_idx});
    tLlon = str2double(splitLine{lon_idx});
    cellSize = str2double(splitLine{dgt_idx});
    lonGrid = str2double(splitLine{nx_idx});
    latGrid = str2double(splitLine{ny_idx});
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sCC_err = 1;
end

if(sCC_err==0)
    disp(['[' datestr(now) '] - - ' 'curAscGridSpec.m successfully executed.']);
end

return

