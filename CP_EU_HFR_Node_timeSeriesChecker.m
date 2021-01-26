%% CP_EU_HFR_Node_timeSeriesChecker.m
% This application checks the completeness of the HFR total and radial data
% time series and plots the time series in time for the selected HFR
% network.

% Author: Lorenzo Corgnati
% Date: January 19, 2021

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');
clear all
close all
clc
EHNtsc_err = 0;
disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_timeSeriesChecker started.']);

%% Setup

% Setup netCDF toolbox
setup_nctoolbox;
% Setup JBDC driver for MySQL
% javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');
javaaddpath('/home/radarcombine/Libraries/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

%%

%% Select time series to check

% Set the network name
networkID = 'HFR-TirLig';
% Set the format version
format_version = 'v2.2';

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

try
    %% Create the folder and filenames for the time series check results
    
    % Set the date ofthe check
    checkDate = datestr(now,'dd-mmm-yyyy');
    % Set the folder name
    checkFolder = ['TimeSeries_Check' filesep networkID filesep checkDate];
    % Create the folder
    if (exist(checkFolder, 'dir') ~= 7)
        mkdir(checkFolder);
    end
    % Set the filenames for total and radial time series check results
    TresFile = [checkFolder filesep 'totalTimeSeries.txt'];
    RresFile = [checkFolder filesep 'radialTimeSeries.txt'];
    plotResFile = [checkFolder filesep 'timeSeriesPlot.png'];
    figResFile = [checkFolder filesep 'timeSeriesPlot.fig'];
    
    %%
    
    %% Connect to database
    
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
    
    %%
    
    %% Query the database for retrieving data about network and stations
    
    % Set and exectute the query to network_tb
    network_selectquery = ['SELECT * FROM network_tb WHERE network_id = ''' networkID ''''];
    network_curs = exec(conn,network_selectquery);
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table for retrieving data of ' networkID ' network successfully executed.']);
    
    % Fetch data
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Data of ' networkID ' network successfully fetched from network_tb table.']);
    
    % Retrieve column names
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
    
    % Retrieve the number of networks
    numNetworks = rows(network_curs);
    
    % Close cursor
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
    
    % Set and exectute the query to station_tb
    station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' networkID ''''];
    station_curs = exec(conn,station_selectquery);
    disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of ' networkID ' network successfully executed.']);
    
    % Fetch data
    station_curs = fetch(station_curs);
    station_data = station_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Data of the stations of ' networkID ' network successfully fetched from station_tb table.']);
    
    % Retrieve column names
    station_columnNames = columnnames(station_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
    
    % Retrieve the number of stations belonging to the current network
    numStations = rows(station_curs);
    disp(['[' datestr(now) '] - - ' 'Number of stations belonging to ' networkID ' network successfully retrieved from station_tb table.']);
    
    % Close cursor to station_tb table
    close(station_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
    
    %%
    
    %% Check time series for total data
    
    % Find the index of the output file path field
    ToutputPathIndexC = strfind(network_columnNames, 'total_HFRnetCDF_folder_path');
    ToutputPathIndex = find(not(cellfun('isempty', ToutputPathIndexC)));
    
    %% DEBUG
    
%     network_data{ToutputPathIndex} = ['../' networkID '/Totals_nc'];
    
    %%
    
    % List nc files to be checked
    TncFiles = rdir([network_data{ToutputPathIndex} filesep format_version filesep '**' filesep '*.nc']);
    
    if(~isempty(TncFiles))
        % Read time stamps of present data
        TdataTS = NaN .* ones(size(TncFiles));
        for tot_idx=1:length(TncFiles)         
            curTimeStr = TncFiles(tot_idx).name(end-17:end-3);
            yyyy = str2double(curTimeStr(1:4));
            mm = str2double(curTimeStr(6:7));
            dd = str2double(curTimeStr(9:10));
            hh = str2double(curTimeStr(12:13));
            mi = str2double(curTimeStr(14:15));
            curTS = datenum(yyyy, mm, dd, hh, mi, 0);
            TdataTS(tot_idx) = curTS;
        end
        
        % Find the index of the data temporal resolution field
        TtemporalResolutionIndexC = strfind(network_columnNames, 'temporal_resolution');
        TtemporalResolutionIndex = find(not(cellfun('isempty', TtemporalResolutionIndexC)));
        
        % Create the complete time frame for the time series
        TstartTime = datetime(datevec(min(TdataTS)));
        TendTime = datetime(datevec(max(TdataTS)));
        TtempResMinutes = minutes(network_data{TtemporalResolutionIndex});
%         TtimeFrame = datenum(TstartTime:TtempResMinutes:TendTime);
        TtimeFrame = TstartTime:TtempResMinutes:TendTime;
        minTS = TtimeFrame(1);
        maxTS = TtimeFrame(end);        
        
        % Find missing dates in the time series
        TdataDates = round(datevec(TdataTS));
        TfullDates = round(datevec(TtimeFrame));
        [~, Tindex] = ismember(TdataDates, TfullDates, 'rows');
        Tfull = zeros(size(TtimeFrame));
        Tfull(Tindex) = 1;
        
        % Save the list of missing dates in the time series
        TmissingDates = TtimeFrame(Tfull==0);
        Tfid = fopen(TresFile,'w');
        fprintf(Tfid,'%s\n',networkID);
        if(isempty(TmissingDates))            
            fprintf(Tfid,'%s',['Total time series COMPLETE in the period ' datestr(TstartTime) ' to ' datestr(TendTime)]);
        else
            fprintf(Tfid,'%s\n\n',['Total time series NOT COMPLETE in the period ' datestr(TstartTime) ' to ' datestr(TendTime)]);
            fprintf(Tfid,'%s\n','The missing dates are:');
            for ii=1:length(TmissingDates)
                fprintf(Tfid,'%s\n',datestr(TmissingDates(ii)));
            end
        end
        fclose(Tfid);
        
        disp(['[' datestr(now) '] - - ' 'Total data time series for network ' networkID ' successfully checked.']);
        
    else
        disp(['[' datestr(now) '] - - ' 'No total data found for ' networkID ' network.']);
    end
    
    %% Check time series for radial data
    
    % Find the index of the station id field
    station_idIndexC = strfind(station_columnNames, 'station_id');
    station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
        
    % Find the index of the output file path field
    RoutputPathIndexC = strfind(station_columnNames, 'radial_HFRnetCDF_folder_path');
    RoutputPathIndex = find(not(cellfun('isempty', RoutputPathIndexC)));
    
    % Scan the radial stations and check time series for each one
    for st_idx=1:numStations
        %% DEBUG
        
%         station_data{st_idx,RoutputPathIndex} = ['../' networkID '/Radials_nc'];
        
        %%
        
        % Retrieve the station id
        stationID = station_data{st_idx,station_idIndex};
        
        % List nc files to be checked
        RncFiles = rdir([station_data{st_idx,RoutputPathIndex} filesep format_version filesep stationID filesep '**' filesep '*.nc']);
        
        if(~isempty(RncFiles))
            % Read time stamps of present data
            RdataTS = NaN .* ones(size(RncFiles));
            for rad_idx=1:length(RncFiles)
                curTimeStr = RncFiles(rad_idx).name(end-17:end-3);
                yyyy = str2double(curTimeStr(1:4));
                mm = str2double(curTimeStr(6:7));
                dd = str2double(curTimeStr(9:10));
                hh = str2double(curTimeStr(12:13));
                mi = str2double(curTimeStr(14:15));
                curTS = datenum(yyyy, mm, dd, hh, mi, 0);
                RdataTS(rad_idx) = curTS;
            end
            
            % Find the index of the data temporal resolution field
            RtemporalResolutionIndexC = strfind(station_columnNames, 'temporal_resolution');
            RtemporalResolutionIndex = find(not(cellfun('isempty', RtemporalResolutionIndexC)));
            
            % Create the complete time frame for the time series
            RstartTime = datetime(datevec(min(RdataTS)));
            RendTime = datetime(datevec(max(RdataTS)));
            RtempResMinutes = minutes(station_data{st_idx,RtemporalResolutionIndex});
%             RtimeFrame = datenum(RstartTime:RtempResMinutes:RendTime);
            RtimeFrame = RstartTime:RtempResMinutes:RendTime;
            
            % Check minimum and maximum time stamps for plots
            if(RtimeFrame(1)<minTS)
                minTS = RtimeFrame(1);
            end
            if(RtimeFrame(end)>maxTS)
                maxTS = RtimeFrame(end);
            end
            
            % Find missing dates in the time series
            RdataDates = round(datevec(RdataTS));
            RfullDates = round(datevec(RtimeFrame));
            [~, Rindex] = ismember(RdataDates, RfullDates, 'rows');
            Rfull = zeros(size(RtimeFrame));
            Rfull(Rindex) = 1;
            
            % Save the list of missing dates in the time series
            RmissingDates = RtimeFrame(Rfull==0);
            Rfid = fopen(RresFile,'a');
            fprintf(Rfid,'%s\n',[networkID '-' stationID]);
            if(isempty(RmissingDates))
                fprintf(Rfid,'%s\n\n\n',['Radial time series COMPLETE in the period ' datestr(RstartTime) ' to ' datestr(RendTime)]);
            else
                fprintf(Rfid,'%s\n\n',['Radial time series NOT COMPLETE in the period ' datestr(RstartTime) ' to ' datestr(RendTime)]);
                fprintf(Rfid,'%s\n','The missing dates are:');
                for jj=1:length(RmissingDates)
                    fprintf(Rfid,'%s\n',datestr(RmissingDates(jj)));
                end
                fprintf(Rfid,'%s\n\n','');
            end
            fclose(Rfid);
            
            disp(['[' datestr(now) '] - - ' 'Radial data time series for station ' stationID ' successfully checked.']);
            
            % Record data for plots
            radialTimeSeries(st_idx).RtimeFrame = RtimeFrame;
            radialTimeSeries(st_idx).Rfull = Rfull;
            
            clear RdataTS RstartTime RendTime RtempResMinutes RtimeFrame RdataDates RfullDates Rfull Rindex RmissingDates;
            
        else
            disp(['[' datestr(now) '] - - ' 'No radial data found for station ' stationID '.']);
        end
    end
    
    %%
    
    %% Create and save plot
    
    % Set number of subfigure
    numSubPlot = numNetworks + numStations;
    % Create plot
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(numSubPlot,1,1);
        
    % Plot the timeseries temporal evolution for total data
    plot(TtimeFrame,Tfull,'-b','linewidth',2);
    hold on;
    title(['NRT total data time series for ' networkID ' network'], 'FontSize', 16);
    xlabel('Dates');
    ylabel('Time series');
    xlim([minTS maxTS]);
%     datetick('x','mmm-yyyy','keepticks');
    
    % Plot the timeseries temporal evolution for radial data
    for st_idx=1:length(radialTimeSeries)
        subplot(numSubPlot,1,st_idx+1);
        plot(radialTimeSeries(st_idx).RtimeFrame,radialTimeSeries(st_idx).Rfull,'-r','linewidth',2);
        hold on;
        title(['NRT radial data time series for ' station_data{st_idx,station_idIndex} ' station'], 'FontSize', 16);
        xlabel('Dates');
        ylabel('Time series');
        xlim([minTS maxTS]);
%         datetick('x','mmm-yyyy','keepticks');
    end
    
    % Save the time series plot
    saveas(gcf,plotResFile);
    savefig(gcf,figResFile);
    close;
    
    %%   
    
    %% Close connection to database
    
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
    
    %%
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    EHNtsc_err = 1;
end

%%
if(EHNtsc_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_timeSeriesChecker successfully executed.']);
else
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_timeSeriesChecker exited with an error.']);
end
%%