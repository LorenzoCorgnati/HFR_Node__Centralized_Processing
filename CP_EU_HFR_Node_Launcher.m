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
% javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');
javaaddpath('/home/radarcombine/Libraries/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

% Setup m_map folder path
% mMapFolder = '/home/lorenz/Toolboxes/Matlab_HFR_AddOn/m_map';
mMapFolder = '/home/radarcombine/Libraries/Matlab_HFR_AddOn/m_map';

EHNP_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Launcher started.']);

%%

%% Set database parameters -- INSERT DATABASE user AND password

sqlConfig.user = 'XXXXXXX';
sqlConfig.password = 'XXXXXXXX';
sqlConfig.host = '150.145.136.104';
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
%         network_selectquery = 'SELECT * FROM `network_tb` WHERE `network_id` LIKE ''%_TEST''';
%         network_selectquery = 'SELECT * FROM `network_tb` WHERE `network_id` LIKE ''HFR-US_%''';
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

%% Launch the functions for processing the networks in simultaneous Matlab processes

if(EHNP_err==0)
    % Find the index of the network_id field
    try
        network_idIndexC = strfind(network_columnNames, 'network_id');
        network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
        
        % Create a cluster object
        parallel.defaultClusterProfile('local');
        c = parcluster();
        % Create and submit one job per network
        for network_idx=1:numNetworks
            jobs{network_idx} = createJob(c); % cell array containing the jobs
            jobs{network_idx}.AttachedFiles = {mMapFolder};
            % Create tasks (one per network)
            createTask(jobs{network_idx}, @CP_processNetwork, 1, {network_data{network_idx,network_idIndex}}, 'CaptureDiary',true);
            % Submit the job to the cluster
            submit(jobs{network_idx});
            disp(['[' datestr(now) '] - - ' 'Job for processing ' network_data{network_idx,network_idIndex} ' network submitted.']);
            %             pause(60);
        end
        
        % Continuously check if a job finished or failed and relaunch it
        kk = 5;
        while(kk>0)
            % Find finished or failed jobs and tasks
            jIndex = [];
            for job_idx=1:length(jobs)
                % Jobs
                if((contains(jobs{job_idx}.State,'finished')) || (contains(jobs{job_idx}.State,'failed')))
                    if((contains(jobs{job_idx}.Tasks.State,'finished')) || (contains(jobs{job_idx}.Tasks.State,'failed')))
                        jIndex = [jIndex;job_idx];
                    end
                end
            end
            
            % Retrieve input arguments for jobs to be relaunched
            jIndex = unique(jIndex);
            for job_idx=1:length(jIndex)
                finishedNetworks{job_idx} = cell2mat([jobs{jIndex(job_idx)}.Tasks.InputArguments]);
            end            
            
            % Print logs and clean the job cell array for finished jobs
            if(~isempty(jIndex))
                for job_idx=length(jIndex):-1:1
                    if (exist(['/var/log/EU_HFR_NODE_Processor_CP/' cell2mat([jobs{jIndex(job_idx)}.Tasks.InputArguments])], 'dir') ~= 7)
                        mkdir(['/var/log/EU_HFR_NODE_Processor_CP/' cell2mat([jobs{jIndex(job_idx)}.Tasks.InputArguments])]);
                    end                    
                    logFilename = ['/var/log/EU_HFR_NODE_Processor_CP/' cell2mat([jobs{jIndex(job_idx)}.Tasks.InputArguments]) filesep 'EU_HFR_Node_CP-' cell2mat([jobs{jIndex(job_idx)}.Tasks.InputArguments]) '.log'];
                    fid = fopen(logFilename,'a');
                    fprintf(fid,'%s\n',jobs{jIndex(job_idx)}.Tasks.Diary);
                    fclose(fid);
                    delete(jobs{jIndex(job_idx)});
                    jobs(jIndex(job_idx)) = [];
                end
            end
            
            % Relaunch finished jobs
            if(~isempty(jIndex))
                for job_idx=1:length(jIndex)
                    % Create e job
                    newJob_idx = length(jobs)+1;
                    jobs{newJob_idx} = createJob(c);
                    jobs{newJob_idx}.AttachedFiles = {mMapFolder};
                    % Create task
                    createTask(jobs{newJob_idx}, @CP_processNetwork, 1, {finishedNetworks{job_idx}}, 'CaptureDiary',true);
                    % Submit new job
                    submit(jobs{newJob_idx});
                    disp(['[' datestr(now) '] - - ' 'Job for processing ' finishedNetworks{job_idx} ' network resubmitted.']);
                    %                 pause(60);
                end
            end
            
            clear finishedNetworks;
            pause (15);
            
        end
        
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
end

%%

if(EHNP_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Launcher successfully executed.']);
end
