%% curAscHeaderDataTable.m
% This function retrieves the file header and the data table from WERA
% cur_asc total files.

% INPUT:
%         ascFile: cell array of strings containing the WERA cur_asc total file.

% OUTPUT:
%         cAHDT_err: error flag (0 = correct, 1 = error)
%         header: cell array of strings containing the file header
%         columnNames: cell array of strings containing the data table cloumn names
%         dataTable: matrix containing the data table

% Author: Lorenzo Corgnati
% Date: October 5, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [cAHDT_err, header, columnNames, dataTable] = curAscHeaderDataTable(ascFile)

disp(['[' datestr(now) '] - - ' 'curAscHeaderDataTable.m started.']);

cAHDT_err = 0;

warning('off', 'all');

try
    % Read the file and look for the data table column names
    % (they are the beginning of the data tabel, thus separating header and data table)
    for line_idx=1:length(ascFile)
        splitLine = regexp(ascFile{line_idx}, '[ \t]+', 'split');
        if(length(splitLine)>1)
            if((any(strcmp(splitLine,'IX'))) && (any(strcmp(splitLine,'IY'))) && (any(strcmp(splitLine,'U[m/s]'))) && (any(strcmp(splitLine,'V[m/s]'))) && (any(strcmp(splitLine,'Acc_U[m/s]'))) && (any(strcmp(splitLine,'Acc_V[m/s]'))))
                endHeader = line_idx - 1;
                % Delete empty lines before data
                data_idx = line_idx + 1;
                expressionNum = '([0-9]+)';
                [startIndexNum,endIndexNum] = regexp(ascFile{data_idx},expressionNum);
                while(isempty(startIndexNum))
                    data_idx = data_idx + 1;
                    [startIndexNum,endIndexNum] = regexp(ascFile{data_idx},expressionNum);
                end
                startDataTable = data_idx;
                break;
            end
        end
    end
    % Retrieve header
    header = ascFile(1:endHeader);
    % Retrieve data table column names
    splitColumnNames = regexp(ascFile(line_idx), '[ \t]+', 'split');
    columnNames = splitColumnNames{1,1};
    % Retrieve data table
    dataTableCellArrayStr = ascFile(startDataTable:length(ascFile));
    dataTable = NaN.*ones(length(dataTableCellArrayStr),length(columnNames));
    for dT_idx=1:length(dataTableCellArrayStr)
        [startIndexNum,endIndexNum] = regexp(dataTableCellArrayStr{dT_idx},expressionNum);
        if(~isempty(startIndexNum))
            splitRow = regexp(dataTableCellArrayStr{dT_idx}, '[ \t]+', 'split');
            for sR_idx=1:length(splitRow)
                dataTable(dT_idx,sR_idx) = str2double(splitRow{sR_idx});
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    cAHDT_err = 1;
end

if(cAHDT_err==0)
    disp(['[' datestr(now) '] - - ' 'curAscHeaderDataTable.m successfully executed.']);
end

return

