%% dataTableColumnNames.m
% This function retrieves the column names of the data table from the header of
% Codar tuv files.

% INPUT:
%         header: header of the Codar ruv or tuv file.
%         LLUVSpec: string containing the LLUVSpec version.

% OUTPUT:
%         dTCN_err: error flag (0 = correct, 1 = error)
%         tableFields: string array containing the site codes

% Author: Lorenzo Corgnati
% Date: July 11, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [dTCN_err,tableFields] = dataTableColumnNames(header,LLUVSpec)

disp(['[' datestr(now) '] - - ' 'dataTableColumnNames.m started.']);

dTCN_err = 0;

warning('off', 'all');

try
    % Retrieve information from header
    for header_idx=1:length(header)
        splitLine = regexp(header{header_idx}, ' ', 'split');
        
        % Retrieve site codes and coordinates
        if(strcmp(splitLine{1}, '%TableType:'))
            if(strcmp(splitLine{2}, 'LLUV'))
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
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    dTCN_err = 1;
end

if(dTCN_err==0)
    disp(['[' datestr(now) '] - - ' 'dataTableColumnNames.m successfully executed.']);
end

return

