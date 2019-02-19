%% tuvSiteCodeCoord.m
% This function retrieves the site codes and coordinates from the header of
% Codar tuv files.

% INPUT:
%         header: header of the Codar ruv or tuv file.
%         LLUVSpec: string containing the LLUVSpec version.

% OUTPUT:
%         sCC_err: error flag (0 = correct, 1 = error)
%         siteCodes: string array containing the site codes
%         siteLat: array containing the site codes
%         siteLon: array containing the site codes

% Author: Lorenzo Corgnati
% Date: July 11, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [sCC_err,siteCodes,siteLat,siteLon] = tuvSiteCodeCoord(header,LLUVSpec)

disp(['[' datestr(now) '] - - ' 'tuvSiteCodeCoord.m started.']);

sCC_err = 0;

warning('off', 'all');

try
    % Retrieve information from header
    if(str2double(LLUVSpec)<1.17)
        for header_idx=1:length(header)
            splitLine = regexp(header{header_idx}, ' ', 'split');
            
            % Retrieve site codes and coordinates
            if(strcmp(splitLine{1}, '%SiteSource:'))
                SiteSource = strrep(header{header_idx}(length('%SiteSource:')+2:length(header{header_idx})), '"', '');
                siteCellArray = strsplit(SiteSource);
                siteCodes(str2double(siteCellArray{2}),:) = siteCellArray{3};
                siteLat(str2double(siteCellArray{2})) = str2double(siteCellArray{4});
                siteLon(str2double(siteCellArray{2})) = str2double(siteCellArray{5});
            end
        end
    else
        for header_idx=1:length(header)
            splitLine = regexp(header{header_idx}, ' ', 'split');
            
            % Retrieve MRGS (merge sources) table structure
            if(strcmp(splitLine{1}, '%TableType:'))
                if(strcmp(splitLine{2}, 'MRGS'))
                    while(~strcmp(splitLine{1}, '%TableColumnTypes:'))
                        header_idx = header_idx + 1;
                        splitLine = regexp(header{header_idx}, ' ', 'split');
                    end
                    break
                end
            end
        end
        TableColumnTypes = strrep(header{header_idx}(length('%TableColumnTypes:')+2:length(header{header_idx})), '"', '');
        tableFields = strsplit(TableColumnTypes);
        % Retrieves fields indices
        %             IndexC = strfind(tableFields, 'SNDX');
        %             Index = find(not(cellfun('isempty', IndexC)));
        CodeC = strfind(tableFields, 'SITE');
        Code = find(not(cellfun('isempty', CodeC)));
        LatC = strfind(tableFields, 'OLAT');
        Lat = find(not(cellfun('isempty', LatC)));
        LonC = strfind(tableFields, 'OLON');
        Lon = find(not(cellfun('isempty', LonC)));
        % Find the MRGS table
        while(~strcmp(splitLine{1}, '%TableStart:'))
            header_idx = header_idx + 1;
            splitLine = regexp(header{header_idx}, ' ', 'split');
        end
        header_idx = header_idx + 1;
        splitLine = regexp(header{header_idx}, ' ', 'split');
        while(strcmp(splitLine{1}, '%%'))
            splitLine = regexp(header{header_idx}, ' ', 'split');
            header_idx = header_idx + 1;
        end
        site_idx = 0;
        while(~strcmp(splitLine{1}, '%TableEnd:'))
            site_idx = site_idx + 1;
            % Remove empty cells and cells containing '%'
            mrgsRow = splitLine(~cellfun('isempty',splitLine));
            mrgsRow(1)=[];
            % Build site codes and coordinates arrays
            tmpCode = mrgsRow{Code};
            siteCodes(site_idx,:) = strrep(tmpCode,'"','');
            siteLat(site_idx) = str2double(mrgsRow{Lat});
            siteLon(site_idx) = str2double(mrgsRow{Lon});
            splitLine = regexp(header{header_idx}, ' ', 'split');
            header_idx = header_idx + 1;
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sCC_err = 1;
end

if(sCC_err==0)
    disp(['[' datestr(now) '] - - ' 'tuvSiteCodeCoord.m successfully executed.']);
end

return

