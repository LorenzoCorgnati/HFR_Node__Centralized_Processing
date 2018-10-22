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

disp(['[' datestr(now) '] - - ' 'CP_TotalConversion_main started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
end

%%

%% Query the database for retrieving network data

% Set and exectute the query
try
    network_selectquery = 'SELECT * FROM network_tb WHERE EU_HFR_processing_flag=1';
    network_curs = exec(conn,network_selectquery);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table successfully executed.']);
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Data from network_tb table successfully fetched.']);
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Number of networks from network_tb table successfully retrieved.']);
end

% Close cursor
try
    close(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
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
        try
            toBeConvertedTotals_selectquery = ['SELECT * FROM total_input_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND NRT_processed_flag = 0'];
            toBeConvertedTotals_curs = exec(conn,toBeConvertedTotals_selectquery);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'Query to total_input_tb table successfully executed.']);
        end
        
        % Fetch data
        try
            toBeConvertedTotals_curs = fetch(toBeConvertedTotals_curs);
            toBeConvertedTotals_data = toBeConvertedTotals_curs.Data;
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'Data from total_input_tb table successfully fetched.']);
        end
        
        % Retrieve column names
        try
            toBeConvertedTotals_columnNames = columnnames(toBeConvertedTotals_curs,true);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'Column names from total_input_tb table successfully retrieved.']);
        end
        
        % Retrieve the number of totals to be converted for the current network
        try
            numToBeConvertedTotals = rows(toBeConvertedTotals_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'Number of totals to be converted per the current network from total_input_tb table successfully retrieved.']);
        end
        
        % Close cursor to total_input_tb table
        try
            close(toBeConvertedTotals_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'Cursor to total_input_tb table successfully closed.']);
        end
        
        try
            % Find the index of the filename field
            filenameIndexC = strfind(toBeConvertedTotals_columnNames, 'filename');
            filenameIndex = find(not(cellfun('isempty', filenameIndexC)));
            
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
            try
                [yMDF_err,yearFolder,monthFolder,dayFolder] = yearMonthDayFolder(toBeConvertedTotals_data{toBeConverted_idx,timestampIndex});
                if (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, 'tuv')) % Codar data
                    [TC_err, network_data(network_idx,:), outputFilename,outputFilesize] = tuv2netCDF_v31([network_data{network_idx,inputPathIndex} filesep dayFolder filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},network_data(network_idx,:),network_columnNames);
                elseif (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, 'cur_asc')) % WERA data
                    [TC_err, network_data(network_idx,:), outputFilename,outputFilesize] = curAsc2netCDF_v31([network_data{network_idx,inputPathIndex} filesep dayFolder filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},network_data(network_idx,:),network_columnNames);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            if (TC_err == 0)
                disp(['[' datestr(now) '] - - ' 'File ' toBeConvertedTotals_data{toBeConverted_idx,filenameIndex} ' successfully converted.']);
            end
            
            % Update NRT_processed_flag in total_input_tb table
            try
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
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            if(TC_err==0)
                disp(['[' datestr(now) '] - - ' 'total_input_tb table successfully updated with NRT processed flag.']);
            end
            
            % Insert converted total info in total_HFRnetCDF_tb table
            try
                % Evaluate datetime from, Time Stamp
                [t2d_err,DateTime] = timestamp2datetime(toBeConvertedTotals_data{toBeConverted_idx,timestampIndex});
                
                % Define a cell array containing the column names to be added
                addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'filesize' 'input_filename' 'check_flag'};
                
                % Define a cell array that contains the data for insertion
                addData = {outputFilename,network_data{network_idx,network_idIndex},toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},DateTime,outputFilesize,toBeConvertedTotals_data{toBeConverted_idx,filenameIndex},0};
                
                % Append the product data into the total_HFRnetCDF_tb table on the database.
                tablename = 'total_HFRnetCDF_tb';
                datainsert(conn,tablename,addColnames,addData);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            if(TC_err==0)
                disp(['[' datestr(now) '] - - ' 'Total converted file information successfully inserted into total_HFRnetCDF_tb table.']);
            end
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
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end
if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
end

%%

pause(2700);