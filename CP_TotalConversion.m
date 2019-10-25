%% CP_TotalConversion_main.m
% This application reads the HFR database for collecting information about
% the total data files to be converted into the European standard data
% model and calls the conversion functions.

% Author: Lorenzo Corgnati
% Date: October 16, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

TC_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_TotalConversion started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
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
    TC_err = 1;
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
    disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

% Close cursor
try
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

%%

%% Scan the networks and convert the related total files

try
    % Find the index of the network_id field
    network_idIndexC = strfind(network_columnNames, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the index of the input file path field
    inputPathIndexC = strfind(network_columnNames, 'total_input_folder_path');
    inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        TC_err = 0;
        % Retrieve information on the stations belonging to the current network
        try
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
            station_curs = exec(conn,station_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of the ' network_data{network_idx,network_idIndex} ' network successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Fetch data
        try
            station_curs = fetch(station_curs);
            station_data = station_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the stations of the ' network_data{network_idx,network_idIndex} ' network successfully fetched from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve column names
        try
            station_columnNames = columnnames(station_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve the number of stations belonging to the current network
        try
            numStations = rows(station_curs);
            disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the ' network_data{network_idx,network_idIndex} ' network successfully retrieved from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Close cursor to station_tb table
        try
            close(station_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve the total files to be converted
        try
            toBeConvertedTotals_selectquery = ['SELECT * FROM total_input_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND NRT_processed_flag = 0'];
            toBeConvertedTotals_curs = exec(conn,toBeConvertedTotals_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to total_input_tb table for retrieving the total files from ' network_data{network_idx,network_idIndex} ' network to be converted successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Fetch data
        try
            toBeConvertedTotals_curs = fetch(toBeConvertedTotals_curs);
            toBeConvertedTotals_data = toBeConvertedTotals_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the total files from ' network_data{network_idx,network_idIndex} ' network to be converted successfully fetched from total_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve column names
        try
            toBeConvertedTotals_columnNames = columnnames(toBeConvertedTotals_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from total_input_tb table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve the number of totals to be converted for the current network
        try
            numToBeConvertedTotals = rows(toBeConvertedTotals_curs);
            disp(['[' datestr(now) '] - - ' 'Number of total files from ' network_data{network_idx,network_idIndex} ' network to be converted successfully retrieved from total_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Close cursor to total_input_tb table
        try
            close(toBeConvertedTotals_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to total_input_tb table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        try
            % Find the index of the filename field
            filenameIndexC = strfind(toBeConvertedTotals_columnNames, 'filename');
            filenameIndex = find(not(cellfun('isempty', filenameIndexC)));
            
            % Find the index of the filepath field
            filepathIndexC = strfind(toBeConvertedTotals_columnNames, 'filepath');
            filepathIndex = find(not(cellfun('isempty', filepathIndexC)));
            
            % Find the index of the extension field
            extensionIndexC = strfind(toBeConvertedTotals_columnNames, 'extension');
            extensionIndex = find(not(cellfun('isempty', extensionIndexC)));
            
            % Find the index of the timestamp field
            timestampIndexC = strfind(toBeConvertedTotals_columnNames, 'timestamp');
            timestampIndex = find(not(cellfun('isempty', timestampIndexC)));
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Scan the tuv files to be converted
        for toBeConverted_idx=1:numToBeConvertedTotals
            TC_err = 0;
            try
                if (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.tuv')) % Codar data
                    % v2.1.1
                    [TC_err, network_data(network_idx,:), outputFilename,outputFilesize] = tuv2netCDF_v32([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},network_data(network_idx,:),network_columnNames,station_data,station_columnNames);
                elseif (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.cur_asc')) % WERA data
                    % v2.1.1
                    [TC_err, network_data(network_idx,:), outputFilename,outputFilesize] = curAsc2netCDF_v32([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},network_data(network_idx,:),network_columnNames,station_data,station_columnNames);
                elseif (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.asc')) % WERA data
                    % v2.1.1
                    [TC_err, network_data(network_idx,:), outputFilename, outputFilesize] = ascTot2netCDF_v32([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},network_data(network_idx,:),network_columnNames,station_data,station_columnNames);
                end
                disp(['[' datestr(now) '] - - ' outputFilename ' total netCDF v2.1.1 file successfully created and stored.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Update NRT_processed_flag in total_input_tb table
            try
                if(TC_err==0)
                    % Find the index of the NRT_processed_flag field
                    NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeConvertedTotals_columnNames, 'NRT_processed_flag'))));
                    % Define a cell array containing the column names to be updated
                    updateColnames = {'NRT_processed_flag'};
                    
                    % Define a cell array that contains the data for insertion
                    updateData = {1};
                    
                    % Update the total_input_tb table on the database
                    tablename = 'total_input_tb';
                    whereclause = ['WHERE filename = ' '''' toBeConvertedTotals_data{toBeConverted_idx,filenameIndex} ''''];
                    update(conn,tablename,updateColnames,updateData,whereclause);
                    disp(['[' datestr(now) '] - - ' 'total_input_tb table successfully updated with NRT processed flag for timestamp ' toBeConvertedTotals_data{toBeConverted_idx,timestampIndex} '.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Insert converted total info in total_HFRnetCDF_tb table
            try
                if((TC_err==0) && (exist('outputFilename','var') ~= 0))
                    % Evaluate datetime from, Time Stamp
                    [t2d_err,DateTime] = timestamp2datetime(toBeConvertedTotals_data{toBeConverted_idx,timestampIndex});
                    
                    % Define a cell array containing the column names to be added
                    addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'input_filename' 'check_flag'};
                    
                    % Define a cell array that contains the data for insertion
                    addData = {outputFilename,network_data{network_idx,network_idIndex},toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),outputFilesize,toBeConvertedTotals_data{toBeConverted_idx,filenameIndex},0};
                    
                    % Append the product data into the total_HFRnetCDF_tb table on the database.
                    tablename = 'total_HFRnetCDF_tb';
                    datainsert(conn,tablename,addColnames,addData);
                    disp(['[' datestr(now) '] - - ' outputFilename 'total file information successfully inserted into total_HFRnetCDF_tb table.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            clear outputFilename outputFilesize;
            
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

%%

%% Close connection

try
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

%%

if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_TotalConversion successfully executed.']);
end

pause(20);