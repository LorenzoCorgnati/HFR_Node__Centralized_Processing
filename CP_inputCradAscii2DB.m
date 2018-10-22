%% CP_inputCradAscii2DB.m
% This application lists the input crad_ascii (WERA radials) files pushed by the HFR data providers
% and insert into the HFR database the information needed for the
% combination of radial files into totals and for the generation of the
% radial and total data files into the European standard data model.

% Author: Lorenzo Corgnati
% Date: October 16, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

iCradDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputCradAscii2DB started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
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
    iCradDB_err = 1;
end
if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table successfully executed.']);
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Data from network_tb table successfully fetched.']);
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Number of networks from network_tb table successfully retrieved.']);
end

% Close cursor
try
    close(network_curs);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
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
    iCradDB_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        try
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
            station_curs = exec(conn,station_selectquery);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        if(iCradDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Query to station_tb table successfully executed.']);
        end
        
        % Fetch data
        try
            station_curs = fetch(station_curs);
            station_data = station_curs.Data;
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        if(iCradDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Data from station_tb table successfully fetched.']);
        end
        
        % Retrieve column names
        try
            station_columnNames = columnnames(station_curs,true);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        if(iCradDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
        end
        
        % Retrieve the number of stations belonging to the current network
        try
            numStations = rows(station_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        if(iCradDB_err==0)
            disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the current network from station_tb table successfully retrieved.']);
        end
        
        % Close cursor to station_tb table
        try
            close(station_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        if(iCradDB_err==0)
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
            iCradDB_err = 1;
        end
        
        % Scan the stations
        for station_idx=1:numStations
            % List the input crad_ascii files for the current station
            try
                cradFiles = rdir([station_data{station_idx,inputPathIndex} '/*/*/*/*.crad_ascii']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                iCradDB_err = 1;
            end
            
            % Insert information about the crad_ascii file into the database (if not yet present)
            for crad_idx=1:length(cradFiles)
                try
                    % Retrieve the filename
                    noFullPathName = cradFiles(crad_idx).name(length(cradFiles(crad_idx).folder)+2:length(cradFiles(crad_idx).name));
                    % Check if the current crad_ascii file is already present on the database
                    dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND filename = ' '''' noFullPathName ''''];
                    dbRadials_curs = exec(conn,dbRadials_selectquery);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iCradDB_err = 1;
                end
                if(iCradDB_err==0)
                    disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table successfully executed.']);
                end
                % Fetch data
                try
                    dbRadials_curs = fetch(dbRadials_curs);
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iCradDB_err = 1;
                end
                if(iCradDB_err==0)
                    disp(['[' datestr(now) '] - - ' 'Data from radial_input_tb table successfully fetched.']);
                end
                
                if(rows(dbRadials_curs) == 0)
                    % Retrieve information about the crad_ascii file
                    try
                        % Load the total file as text
                        cradFile = textread(cradFiles(crad_idx).name,  '%s', 'whitespace', '\n');
                        % Read the file header and look for timestamp
                        for line_idx=1:length(cradFile)
                            splitLine = regexp(cradFile{line_idx}, '[ \t]+', 'split');
                            if(length(splitLine)>1)
                                expressionDate = '([0-9]{2}-[A-Z]{3}-[0-9]{2})';
                                expressionTime = '([0-9]{2}:[0-9]{2})';
                                expressionUTC = 'UTC';
                                [startIndexDate,endIndexDate] = regexp(cradFile{line_idx},expressionDate);
                                [startIndexTime,endIndexTime] = regexp(cradFile{line_idx},expressionTime);
                                [startIndexUTC,endIndexUTC] = regexp(cradFile{line_idx},expressionUTC);
                                if((~isempty(startIndexDate)) && (~isempty(startIndexTime)) && (~isempty(startIndexUTC)))
                                    date = cradFile{line_idx}(startIndexDate:endIndexDate);
                                    time = cradFile{line_idx}(startIndexTime:endIndexTime);
                                    TimeStampVec = datevec([date ' ' time]);
                                    TimeStamp = [num2str(TimeStampVec(1)) ' ' num2str(TimeStampVec(2),'%02d') ' ' num2str(TimeStampVec(3),'%02d') ' ' num2str(TimeStampVec(4),'%02d') ' ' num2str(TimeStampVec(5),'%02d') ' ' num2str(TimeStampVec(6),'%02d')];
                                    break;
                                end
                            end
                        end
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCradDB_err = 1;
                    end
                    
                    % Evaluate datetime from, Time Stamp
                    try
                        [t2d_err,DateTime] = timestamp2datetime(TimeStamp);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCradDB_err = 1;
                    end
                    
                    % Retrieve information about the crad_ascii file
                    try
                        cradFileInfo = dir(cradFiles(crad_idx).name);
                        cradFilesize = cradFileInfo.bytes*0.001;
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCradDB_err = 1;
                    end
                    
                    % Write crad_ascii info in radial_input_tb table
                    try
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'network_id' 'station_id' 'timestamp' 'datetime' 'filesize' 'extension' 'NRT_processed_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {noFullPathName,network_data{network_idx,network_idIndex},station_data{station_idx,station_idIndex},TimeStamp,DateTime,cradFilesize,'crad_ascii',0};
                        
                        % Append the product data into the radial_input_tb table on the database.
                        tablename = 'radial_input_tb';
                        datainsert(conn,tablename,addColnames,addData);
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iCradDB_err = 1;
                    end
                    if(iCradDB_err==0)
                        disp(['[' datestr(now) '] - - ' 'Radial input file information successfully inserted into radial_input_tb table.']);
                    end
                end
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end

%%

%% Close connection

try
    close(conn);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iCradDB_err = 1;
end
if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
end

%%