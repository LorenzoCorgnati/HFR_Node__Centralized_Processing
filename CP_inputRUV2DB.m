%% CP_inputRUV2DB.m
% This application lists the input ruv files pushed by the HFR data providers
% and insert into the HFR database the information needed for the
% combination of radial files into totals and for the generation of the
% radial and total data files into the European standard data model.

% Author: Lorenzo Corgnati
% Date: October 16, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

iRDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputRUV2DB started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
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
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table successfully executed.']);
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Data from network_tb table successfully fetched.']);
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Number of networks from network_tb table successfully retrieved.']);
end

% Close cursor
try
    close(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
end

%%

%% Scan the networks, find the stations, list the related radial files and insert information into the database

try
    % Find the index of the network_id field
    network_idIndexC = strfind(network_columnNames, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        try
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
            station_curs = exec(conn,station_selectquery);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        if(iRDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Query to station_tb table successfully executed.']);
        end
        
        % Fetch data
        try
            station_curs = fetch(station_curs);
            station_data = station_curs.Data;
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        if(iRDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Data from station_tb table successfully fetched.']);
        end
        
        % Retrieve column names
        try
            station_columnNames = columnnames(station_curs,true);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        if(iRDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
        end
        
        % Retrieve the number of stations belonging to the current network
        try
            numStations = rows(station_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        if(iRDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the current network from station_tb table successfully retrieved.']);
        end
        
        % Close cursor to station_tb table
        try
            close(station_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        if(iRDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
        end
        
        try
            % Find the index of the input file path field
            inputPathIndexC = strfind(station_columnNames, 'radial_input_folder_path');
            inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
            
            % Find the index of the station_id field
            station_idIndexC = strfind(station_columnNames, 'station_id');
            station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Scan the stations
        for station_idx=1:numStations
            % List the input ruv files for the current station
            try
                ruvFiles = rdir([station_data{station_idx,inputPathIndex} '/*/*/*/*.ruv']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                iRDB_err = 1;
            end
            
            % Insert information about the ruv file into the database (if not yet present)
            for ruv_idx=1:length(ruvFiles)
                try
                    % Retrieve the filename
                    noFullPathName = ruvFiles(ruv_idx).name(length(ruvFiles(ruv_idx).folder)+2:length(ruvFiles(ruv_idx).name));
                    % Check if the current ruv file is already present on the database
                    dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND filename = ' '''' noFullPathName ''''];
                    dbRadials_curs = exec(conn,dbRadials_selectquery);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iRDB_err = 1;
                end
                if(iRDB_err==0)
                    disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table successfully executed.']);
                end
                % Fetch data
                try
                    dbRadials_curs = fetch(dbRadials_curs);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iRDB_err = 1;
                end
                if(iRDB_err==0)
                    disp(['[' datestr(now) '] - - ' 'Data from radial_input_tb table successfully fetched.']);
                end
                
                if(rows(dbRadials_curs) == 0)
                    % Retrieve information about the ruv file
                    try
                        % Load the radial file as structure
                        radStruct = loadRDLFile(ruvFiles(ruv_idx).name,'false','warning');
                        % Read the file header
                        radHeader = radStruct.OtherMetadata.Header;
                        % Retrieve information from header
                        for header_idx=1:length(radHeader)
                            splitLine = regexp(radHeader{header_idx}, ' ', 'split');
                            % Retrieve TimeStamp
                            if(strcmp(splitLine{1}, '%TimeStamp:'))
                                TimeStamp = strrep(radHeader{header_idx}(length('%TimeStamp:')+2:length(radHeader{header_idx})), '"', '');
                                break;
                            end
                        end
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iRDB_err = 1;
                    end
                    
                    try
                        % Evaluate datetime from, Time Stamp
                        [t2d_err,DateTime] = timestamp2datetime(TimeStamp);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iRDB_err = 1;
                    end
                    
                    % Retrieve information about the ruv file
                    try
                        ruvFileInfo = dir(ruvFiles(ruv_idx).name);
                        ruvFilesize = ruvFileInfo.bytes*0.001;
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iRDB_err = 1;
                    end
                    
                    % Write ruv info in radial_input_tb table
                    try
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'network_id' 'station_id' 'timestamp' 'datetime' 'filesize' 'extension' 'NRT_processed_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {noFullPathName,network_data{network_idx,network_idIndex},station_data{station_idx,station_idIndex},TimeStamp,DateTime,ruvFilesize,'ruv',0};
                        
                        % Append the product data into the radial_input_tb table on the database.
                        tablename = 'radial_input_tb';
                        datainsert(conn,tablename,addColnames,addData);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iRDB_err = 1;
                    end
                    if(iRDB_err==0)
                        disp(['[' datestr(now) '] - - ' 'Radial input file information successfully inserted into radial_input_tb table.']);
                    end
                end
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

%%

%% Close connection

try
    close(conn);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
end

%%