%% curAscSiteCodeCoord.m
% This function retrieves the site codes and coordinates from the header of
% WERA cur_asc total files.

% INPUT:
%         header: header of the WERA cur_asc file.

% OUTPUT:
%         sCC_err: error flag (0 = correct, 1 = error)
%         siteCodes: string array containing the site codes
%         siteLat: array containing the site codes
%         siteLon: array containing the site codes

% Author: Lorenzo Corgnati
% Date: October 5, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [sCC_err,siteCodes,siteLat,siteLon] = curAscSiteCodeCoord(header)

disp(['[' datestr(now) '] - - ' 'curAscSiteCodeCoord.m started.']);

sCC_err = 0;

warning('off', 'all');

try
    % Read the header
    numSites = 0;
    for line_idx=1:length(header)
        splitLine = regexp(header{line_idx}, '[ \t]+', 'split');
        % Find the site lat, lon and codes
        if(length(splitLine)>1)
            if((any(strcmp(splitLine,'North'))) && (any(strcmp(splitLine,'East'))))
                numSites = numSites + 1;
                latLogical = strcmp(splitLine,'North');
                lat_idx = find(latLogical, 1, 'first') - 1;
                lonLogical = strcmp(splitLine,'East');
                lon_idx = find(lonLogical, 1, 'first') - 1;
                siteLogical = strcmp(splitLine,'UTC');
                site_idx = find(siteLogical, 1, 'first') + 1;
                % Fill the output variables
                siteLat(numSites) = str2double(splitLine{lat_idx});
                siteLon(numSites) = str2double(splitLine{lon_idx});
                siteCodes(numSites,:) = upper(splitLine{site_idx}(1:3));
            elseif ((any(strcmp(splitLine,'North'))) && (any(strcmp(splitLine,'West'))))
                numSites = numSites + 1;
                latLogical = strcmp(splitLine,'North');
                lat_idx = find(latLogical, 1, 'first') - 1;
                lonLogical = strcmp(splitLine,'West');
                lon_idx = find(lonLogical, 1, 'first') - 1;
                siteLogical = strcmp(splitLine,'UTC');
                site_idx = find(siteLogical, 1, 'first') + 1;
                % Fill the output variables
                siteLat(numSites) = str2double(splitLine{lat_idx});
                siteLon(numSites) = -str2double(splitLine{lon_idx});
                siteCodes(numSites,:) = upper(splitLine{site_idx}(1:3));
                
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sCC_err = 1;
end

if(sCC_err==0)
    disp(['[' datestr(now) '] - - ' 'curAscSiteCodeCoord.m successfully executed.']);
end

return

