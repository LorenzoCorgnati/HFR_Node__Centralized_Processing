%% CP_inputMetNoRadials2DB.m
% This application lists the input ruv files pushed by the HFR data providers
% and insert into the HFR database the information needed for the
% combination of radial files into totals and for the generation of the
% radial and total data files into the European standard data model.

% Author: Lorenzo Corgnati
% Date: November 4, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

iRDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputMetNoRadials2DB started.']);

startDateNum = datenum(startDate);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

%%

%% Query the database for retrieving network data

% Set and exectute the query
try
    network_selectquery = 'SELECT * FROM network_tb WHERE EU_HFR_processing_flag=1 AND (network_id=''HFR-Skagerrak'' OR network_id=''HFR-Finnmark'')';
%     network_selectquery = 'SELECT * FROM network_tb WHERE (network_id=''HFR-Skagerrak'' OR network_id=''HFR-Finnmark'')';
    network_curs = exec(conn,network_selectquery);
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table for retrieving network data successfully executed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
    disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Close cursor
try
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

%%

%% Scan the networks, find the stations, list the related radial files and insert information into the database

try
    % Find the index of the network_id field
    network_idIndexC = strfind(network_columnNames, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the index of the TDS_root_url field
    TDS_root_urlIndexC = strfind(network_columnNames, 'TDS_root_url');
    TDS_root_urlIndex = find(not(cellfun('isempty', TDS_root_urlIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        iRDB_err = 0;
        try
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
            station_curs = exec(conn,station_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of the ' network_data{network_idx,network_idIndex} ' network successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Fetch data
        try
            station_curs = fetch(station_curs);
            station_data = station_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the stations of the ' network_data{network_idx,network_idIndex} ' network successfully fetched from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Retrieve column names
        try
            station_columnNames = columnnames(station_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Retrieve the number of stations belonging to the current network
        try
            numStations = rows(station_curs);
            disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the ' network_data{network_idx,network_idIndex} ' network successfully retrieved from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Close cursor to station_tb table
        try
            close(station_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
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
            % Build the list of radials for the time period to be processed
            try
                [iRDB_err,matRad] = MetNoTDSradListBuilder(network_data{network_idx,TDS_root_urlIndex},station_data{station_idx,station_idIndex},startDate,station_data{station_idx,inputPathIndex});
                if(iRDB_err == 0)
                    disp(['[' datestr(now) '] - - ' 'Radials files from ' station_data{station_idx,station_idIndex} ' station successfully listed.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                iRDB_err = 1;
            end
            
            % Scan the radials of the time period to be processed
            if(iRDB_err == 0)
                for rad_idx=1:length(matRad.mat)
                    iRDB_err = 0;
                    % Check if the current radial was already remapped
                    try
                        if (exist(matRad.mat{rad_idx}, 'file') ~= 2)
                            % Remap the radial nc file to RUV struct and save the mat file
                            try
                                [iRDB_err,ruvRad] = remapMetNoNCradial2ruv(matRad.nc{rad_idx},matRad.timeIndex{rad_idx});
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                iRDB_err = 1;
                            end
                            
                            % Save radial mat file
                            if(iRDB_err == 0)
                                try
                                    % Retrieve the filename
                                    [pathstr,name,ext]=fileparts(matRad.mat{rad_idx});
                                    noFullPathName=[name ext];
                                    save(matRad.mat{rad_idx}, 'ruvRad');
                                    disp(['[' datestr(now) '] - - ' noFullPathName ' file successfully saved.']);
                                catch err
                                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                    iRDB_err = 1;
                                end
                                
                                % Insert information about the ruv file into the database (if not yet present)
                                try
                                    % Check if the current ruv file is already present on the database
                                    dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND filename = ' '''' noFullPathName ''' ORDER BY timestamp'];
                                    dbRadials_curs = exec(conn,dbRadials_selectquery);
                                    disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table for checking if ' noFullPathName ' radial file is already present in the database successfully executed.']);
                                catch err
                                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                    iRDB_err = 1;
                                end
                                
                                % Fetch data
                                try
                                    dbRadials_curs = fetch(dbRadials_curs);
                                    disp(['[' datestr(now) '] - - ' 'Data about the presence of ' noFullPathName ' radial file in the database successfully fetched from radial_input_tb table.']);
                                catch err
                                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                    iRDB_err = 1;
                                end
                                
                                if(rows(dbRadials_curs) == 0)
                                    % Retrieve information about the ruv file
                                    try
                                        % Read the file header
                                        radHeader = ruvRad.OtherMetadata.Header;
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
                                        ruvFileInfo = dir(matRad.mat{rad_idx});
                                        ruvFilesize = ruvFileInfo.bytes/1024;
                                    catch err
                                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                        iRDB_err = 1;
                                    end
                                    
                                    % Write ruv info in radial_input_tb table
                                    try
                                        % Define a cell array containing the column names to be added
                                        addColnames = {'filename' 'filepath' 'network_id' 'station_id' 'timestamp' 'datetime' 'reception_date' 'filesize' 'extension' 'NRT_processed_flag'};
                                        
                                        % Define a cell array that contains the data for insertion
                                        addData = {noFullPathName,pathstr,network_data{network_idx,network_idIndex},station_data{station_idx,station_idIndex},TimeStamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),ruvFilesize,ext,0};
                                        
                                        % Append the product data into the radial_input_tb table on the database.
                                        tablename = 'radial_input_tb';
                                        datainsert(conn,tablename,addColnames,addData);
                                        disp(['[' datestr(now) '] - - ' noFullPathName ' radial file information successfully inserted into radial_input_tb table.']);
                                    catch err
                                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                        iRDB_err = 1;
                                    end
                                end
                            end
                        end
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        iRDB_err = 1;
                    end
                    
                    clear ruvRad
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
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

%%

if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_inputMetNoRadials2DB successfully executed.']);
end