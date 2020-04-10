%% CP_EU_HFR_Node_Launcher.m
% This wrapper launches the scripts for inserting into the HFR database
% the information about radial and totala files (both Codar and WERA)
% pushed by the data providers and for combining radials into totals and
% converting radials and totals to netCDF files according to the European
% standard data model.

% This version implements parallel computing by launching a separate
% process per each network to be processed.

% Author: Lorenzo Corgnati
% Date: February 7, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

%% Setup

warning('off', 'all');

clear all
close all
clc

% Setup JBDC driver for MySQL
javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

EHNP_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Launcher started.']);

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

%% Connect to database

if(EHNP_err==0)
    try
        conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
        disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
end

%%

%% Query the database for retrieving network data

if(isopen(conn))
    % Set and exectute the query
    try
        network_selectquery = 'SELECT * FROM network_tb WHERE EU_HFR_processing_flag=1';
        network_selectquery = 'SELECT * FROM `network_tb` WHERE `network_id` LIKE ''%_TEST''';
        network_curs = exec(conn,network_selectquery);
        disp(['[' datestr(now) '] - - ' 'Query to network_tb table for retrieving network data successfully executed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
    % Fetch data
    try
        network_curs = fetch(network_curs);
        network_data = network_curs.Data;
        disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
    % Retrieve column names
    try
        network_columnNames = columnnames(network_curs,true);
        disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
    % Retrieve the number of networks
    try
        numNetworks = rows(network_curs);
        disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
    % Close cursor
    try
        close(network_curs);
        disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
end

%%

%% Close connection

try
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    EHNP_err = 1;
end

%%

%% Launch the functions for processing the networks in simultaneous Matlab instances

if(EHNP_err==0)
    % Find the index of the network_id field
    try
        network_idIndexC = strfind(network_columnNames, 'network_id');
        network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        pN_err = 1;
    end
    
    for network_idx=1:numNetworks
        logFilename = ['/var/log/EU_HFR_NODE_Processor_CP/EU_HFR_Node_CP-' network_data{network_idx,network_idIndex} '.log'];
        system(['/usr/local/MATLAB/R2018a/bin/matlab -nodesktop -r ''/mnt/data/CNR/RADAR/Script/MATLAB/HFR_Node__Centralized_Processing/CP_processNetwork(''' network_data{network_idx,network_idIndex} '''); quit'' -logfile ' logFilename ' &']);
        sleep(20);
    end
    
    % DEVO POI CAMBIARE startDate IN CP_processNetwork
    
end

%%

if(pN_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Launcher successfully executed.']);
end

quit