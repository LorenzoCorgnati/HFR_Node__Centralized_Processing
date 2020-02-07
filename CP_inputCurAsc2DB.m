%% CP_inputCurAsc2DB.m
% This application lists the input cur_asc (WERA totals) files pushed by the HFR data providers
% and insert into the HFR database the information needed for the conversion of
% the total data files into the European standard data model.

% Author: Lorenzo Corgnati
% Date: October 3, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

iCurDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputCurAsc2DB started.']);

startDateNum = datenum(startDate);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

%%

%% Query the database for retrieving network data

% Set and exectute the query
try
    network_selectquery = 'SELECT * FROM network_tb WHERE EU_HFR_processing_flag=1';
    network_curs = exec(conn,network_selectquery);
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table for retrieving network data successfully executed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
    disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

% Close cursor
try
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

%%

%% Scan the networks, list the related total files and insert information into the database

try
    % Find the index of the network_id field
    network_idIndexC = strfind(network_columnNames, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the index of the input file path field
    inputPathIndexC = strfind(network_columnNames, 'total_input_folder_path');
    inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        iCurDB_err = 0;
        if(~isempty(network_data{network_idx,inputPathIndex}))
            % Trim heading and trailing whitespaces from folder path
            network_data{network_idx,inputPathIndex} = strtrim(network_data{network_idx,inputPathIndex});
            % List the input cur_asc files
            try
                ascFiles = rdir([network_data{network_idx,inputPathIndex} filesep '**' filesep '*.cur_asc'],'datenum>floor(now-3)');
                disp(['[' datestr(now) '] - - ' 'Total files from ' network_data{network_idx,network_idIndex} ' network successfully listed.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                iCurDB_err = 1;
            end
            % Insert information about the cur_asc file into the database (if not yet present)
            for asc_idx=1:length(ascFiles)
                iCurDB_err = 0;
                % Retrieve the filename
                [pathstr,name,ext]=fileparts(ascFiles(asc_idx).name);
                noFullPathName=[name ext];
                % Check if the current cur_asc file is already present on the database
                try
                    dbTotals_selectquery = ['SELECT * FROM total_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND filename = ' '''' noFullPathName ''' ORDER BY timestamp'];
                    dbTotals_curs = exec(conn,dbTotals_selectquery);
                    disp(['[' datestr(now) '] - - ' 'Query to total_input_tb table for checking if ' noFullPathName ' total file is already present in the database successfully executed.']);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iCurDB_err = 1;
                end
                
                % Fetch data
                try
                    dbTotals_curs = fetch(dbTotals_curs);
                    disp(['[' datestr(now) '] - - ' 'Data about the presence of ' noFullPathName ' total file in the database successfully fetched from total_input_tb table.']);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iCurDB_err = 1;
                end
                                
                if(rows(dbTotals_curs) == 0)
                    % Retrieve information about the cur_asc file
                    try
%                         % Load the total file as text
%                         ascFile = textread(ascFiles(asc_idx).name,  '%s', 'whitespace', '\n');
%                         % Read the file header and look for timestamp
%                         for line_idx=1:length(ascFile)
%                             splitLine = regexp(ascFile{line_idx}, '[ \t]+', 'split');
%                             if(length(splitLine)>1)
%                                 expressionDate = '([0-9]{2}-[A-Z]{3}-[0-9]{4})';
%                                 expressionTime = '([0-9]{2}:[0-9]{2})';
%                                 [startIndexDate,endIndexDate] = regexp(splitLine{1},expressionDate);
%                                 [startIndexTime,endIndexTime] = regexp(splitLine{2},expressionTime);
%                                 if((~isempty(startIndexDate)) && (~isempty(startIndexTime)) && (sum(splitLine{3}=='UTC')==3))
%                                     date = splitLine{1}(startIndexDate:endIndexDate);
%                                     time = splitLine{2}(startIndexTime:endIndexTime);
%                                     TimeStampVec = datevec([date ' ' time]);
%                                     TimeStamp = [num2str(TimeStampVec(1)) ' ' num2str(TimeStampVec(2),'%02d') ' ' num2str(TimeStampVec(3),'%02d') ' ' num2str(TimeStampVec(4),'%02d') ' ' num2str(TimeStampVec(5),'%02d') ' ' num2str(TimeStampVec(6),'%02d')];
%                                     break;
%                                 end
%                             end
%                         end
                        [date,time] = textread(ascFiles(asc_idx).name, '%11c %*0c %5c',1, 'headerlines',1);
                        TimeStampVec = datevec([date ' ' time]);
                        TimeStamp = [num2str(TimeStampVec(1)) ' ' num2str(TimeStampVec(2),'%02d') ' ' num2str(TimeStampVec(3),'%02d') ' ' num2str(TimeStampVec(4),'%02d') ' ' num2str(TimeStampVec(5),'%02d') ' ' num2str(TimeStampVec(6),'%02d')];
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCurDB_err = 1;
                    end
                    
                    % Evaluate datetime from, Time Stamp
                    try
                        [t2d_err,DateTime] = timestamp2datetime(TimeStamp);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCurDB_err = 1;
                    end
                    
                    % Retrieve information about the cur_asc file
                    try
                        ascFileInfo = dir(ascFiles(asc_idx).name);
                        ascFilesize = ascFileInfo.bytes/1024;
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCurDB_err = 1;
                    end
                    
                    % Write cur_asc info in total_input_tb table
                    try
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'filepath' 'network_id' 'timestamp' 'datetime' 'reception_date' 'filesize' 'extension' 'NRT_processed_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {noFullPathName,pathstr,network_data{network_idx,network_idIndex},TimeStamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),ascFilesize,ext,0};
                        
                        % Append the product data into the total_input_tb table on the database.
                        tablename = 'total_input_tb';
                        datainsert(conn,tablename,addColnames,addData);
                        disp(['[' datestr(now) '] - - ' noFullPathName ' total file information successfully inserted into total_input_tb table.']);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCurDB_err = 1;
                    end                                       
                end
                
                % Close cursor to total_input_tb table
                try
                    close(dbTotals_curs);
                    disp(['[' datestr(now) '] - - ' 'Cursor to total_input_tb table successfully closed.']);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iCurDB_err = 1;
                end                
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

%%

%% Close connection

try
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCurDB_err = 1;
end

%%

if(iCurDB_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_inputCurAsc2DB successfully executed.']);
end