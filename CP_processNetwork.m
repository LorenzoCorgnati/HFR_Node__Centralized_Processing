%% CP_processNetwork.m
% This application processes the radial and total files pushed by the HFR
% data providers for generating radial and total files according to the
% European standard data model.

% The first processing step consists in the listing of the input files
% (both radial and total) pushed by the HFR data providers for inserting
% into the HFR database the information needed for the % combination of
% radial files into totals and for the generation of the radial and total
% data files into the European standard data model.

% The second processing step consists in reading the HFR database for
% collecting information about the radial data files to be combined into
% totals and in the combination and generation of radial and total data
% files according to the European standard data model.

% The third processing step consists in reading the HFR database for
% collecting information about the total data files to be converted into
% the European standard data model and in the generation of total data
% files according to the European standard data model.

% INPUT:
%         networkID: network ID of the network to be processed

% OUTPUT:
%         pN_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: February 21, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [pN_err] = CP_processNetwork(networkID)

%% Setup

warning('off', 'all');

% Setup netCDF toolbox
setup_nctoolbox;

% Setup JBDC driver for MySQL
% javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');
javaaddpath('/home/radarcombine/Libraries/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

% Setup map colormap
set(0,'DefaultFigureColormap',feval('jet'));

pN_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_processNetwork started.']);

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    pN_err = 1;
end

%%

%% Query the database for retrieving network data

if(isopen(conn))
    % Set and exectute the query
    try
        network_selectquery = ['SELECT * FROM network_tb WHERE network_id=''' networkID ''''];
        network_curs = exec(conn,network_selectquery);
        disp(['[' datestr(now) '] - - Query to network_tb table for retrieving network data successfully executed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Retrieve column names
    try
        network_curs = fetch(network_curs);
        network_columnNames = columnnames(network_curs,true);
        disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Find the index of the network_id field
    try
        network_idIndexC = strfind(network_columnNames, 'network_id');
        network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Fetch data
    try
        network_data = network_curs.Data;
        disp(['[' datestr(now) '] - - ' networkID ' network data successfully fetched from network_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Close cursor
    try
        close(network_curs);
        disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
end

%%

%% Query the database for retrieving station data

if(isopen(conn))
    % Set and exectute the query
    try
        % Manage the case of the ISMAR-LaMMA integrated network (HFR-WesternItaly)
        if(strcmp(network_data{1,network_idIndex},'HFR-WesternItaly'))
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''HFR-TirLig'' OR network_id = ' '''HFR-LaMMA'''];
        else
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{1,network_idIndex} ''''];
        end
        station_curs = exec(conn,station_selectquery);
        disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of the ' network_data{1,network_idIndex} ' network successfully executed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Fetch data
    try
        station_curs = fetch(station_curs);
        station_data = station_curs.Data;
        disp(['[' datestr(now) '] - - ' 'Data of the stations of the ' network_data{1,network_idIndex} ' network successfully fetched from station_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Retrieve column names
    try
        station_columnNames = columnnames(station_curs,true);
        disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Close cursor to station_tb table
    try
        close(station_curs);
        disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    %%
    
    %% Process data
    
    % Set datetime of the starting date of the processing period
    try
        startDate = startCombinationDate(now);
%         startDate = '2012-01-29';
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    % Check the connection with database
    while (~isopen(conn))
        try
            conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
            disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            pN_err = 1;
        end
    end
    
    % Process
    try
        % Input radials
        if(strcmp(network_data{1,network_idIndex},'HFR-Skagerrak') || strcmp(network_data{1,network_idIndex},'HFR-Finnmark'))
            iRDB_err = CP_inputMetNoRadials2DB(conn,startDate,network_data,network_columnNames,station_data,station_columnNames);
        else
            iRDB_err = CP_inputRUV2DB(conn,startDate,network_data{1,network_idIndex},station_data,station_columnNames);
            iAscRadDB_err = CP_inputAscRad2DB(conn,startDate,network_data{1,network_idIndex},station_data,station_columnNames);
            %     CP_inputCradAscii2DB;
        end
        
        % Combine radials into totals
        HFRC_err = CP_HFRCombiner(conn,startDate,network_data,network_columnNames,station_data,station_columnNames);
        
        % Input totals
        iTDB_err = CP_inputTUV2DB(conn,startDate,network_data,network_columnNames);
        iAscTotDB_err = CP_inputAscTot2DB(conn,startDate,network_data,network_columnNames);
        iCurDB_err = CP_inputCurAsc2DB(conn,startDate,network_data,network_columnNames);
        
        % Convert totals
        TC_err = CP_TotalConversion(conn,startDate,network_data,network_columnNames,station_data,station_columnNames);
        
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    %%
    
end

if(pN_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_processNetwork.m successfully executed.']);
end

return