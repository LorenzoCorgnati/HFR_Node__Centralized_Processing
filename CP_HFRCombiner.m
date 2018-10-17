%% HFRCombiner.m
% This application reads the HFR database for collecting information about
% the radial data files to be combined into totals and performs the
% combination and the generation of radial and total data files into the
% European standard data model.

% Author: Lorenzo Corgnati
% Date: July 17, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

% Setup netCDF toolbox
setup_nctoolbox;

% Setup JBDC driver for MySQL
javaaddpath('/Users/reverendo/Toolboxes/mysql-connector-java-5.1.17.jar');

% Setup map colormap
set(0,'DefaultFigureColormap',feval('jet'));

HFRC_err = 0;

disp(['[' datestr(now) '] - - ' 'HFRCombiner started.']);

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

% Set the infinite loop for continuous operation
kk = 5;
while(kk>0)
    
    %% Set the flag for updating the network_tb table of the database
    network_tbUpdateFlag = 0;
    
    %% Set datetime of the starting date of the combination period
    
    startDate = startCombinationDate(now);
    
    %%
    
    %% Connect to database
    
    try
        conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    if(HFRC_err==0)
        disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
    end
    
    %%
    
    %% Query the database for retrieving network data
    
    if(HFRC_err==0)
        % Set and exectute the query
        try
            network_selectquery = 'SELECT * FROM network_tb WHERE EU_HFR_processing_flag=1';
            network_curs = exec(conn,network_selectquery);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        if(HFRC_err==0)
            disp(['[' datestr(now) '] - - ' 'Query to network_tb table successfully executed.']);
        end
        
        % Fetch data
        try
            network_curs = fetch(network_curs);
            network_data = network_curs.Data;
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        if(HFRC_err==0)
            disp(['[' datestr(now) '] - - ' 'Data from network_tb table successfully fetched.']);
        end
        
        % Retrieve column names
        try
            network_columnNames = columnnames(network_curs,true);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        if(HFRC_err==0)
            disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
        end
        
        % Retrieve the number of networks
        try
            numNetworks = rows(network_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        if(HFRC_err==0)
            disp(['[' datestr(now) '] - - ' 'Number of networks from network_tb table successfully retrieved.']);
        end
        
        % Close cursor
        try
            close(network_curs);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        if(HFRC_err==0)
            disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
        end
    end
    
    %%
    
    %% Scan the networks and combine the related radial files
    
    if(HFRC_err==0)
        % Find the index of the network_id field
        network_idIndexC = strfind(network_columnNames, 'network_id');
        network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
        
        % Find the indices of the geospatial boundaries fields
        geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'geospatial_lat_min'))));
        geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'geospatial_lat_max'))));
        geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'geospatial_lon_min'))));
        geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'geospatial_lon_max'))));
        
        % Find the index of the grid resolution field
        grid_resolutionIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'grid_resolution'))));
        
        % Find the index of the of the combination search radius field
        matPathIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'total_mat_folder_path'))));
        
        % Find the index of the of the output mat files folder path field
        spatthreshIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'combination_search_radius'))));
    end
    
    if(HFRC_err==0)
        % Scan the networks
        for network_idx=1:numNetworks
            % Build the regular LonLat grid given the geographical boundaries and the grid resolution for the radial combination into total
            [gridLon, gridLat] = LonLat_grid([network_data{network_idx,geospatial_lon_minIndex},network_data{network_idx,geospatial_lat_minIndex}], [network_data{network_idx,geospatial_lon_maxIndex},network_data{network_idx,geospatial_lat_maxIndex}], network_data{network_idx,grid_resolutionIndex}, 'km');
            gridLat = flipud(gridLat);
            lon = gridLon(1,:);
            lat = gridLat(:,1);
            length_lon=length(lon);
            length_lat=length(lat);
            for i=1:length_lon
                lonG(1+(i-1)*length_lat:(i-1)*length_lat+length_lat) = lon(i)*ones(1,length_lat);
                latG(1+(i-1)*length_lat:(i-1)*length_lat+length_lat) = lat;
            end
            Grid(:,1) = lonG';
            Grid(:,2) = latG';
            
            disp(['[' datestr(now) '] - - ' 'Grid for radial combination successfully generated.']);
            
            % Build the mask for the total masking
            try
                if (exist(['Total_Masks' filesep network_data{network_idx,network_idIndex}], 'dir') ~= 7)
                    mkdir(['Total_Masks' filesep network_data{network_idx,network_idIndex}]);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            try
                fname = makeCoast([network_data{network_idx,geospatial_lon_minIndex},network_data{network_idx,geospatial_lon_maxIndex}],[network_data{network_idx,geospatial_lat_minIndex},network_data{network_idx,geospatial_lat_maxIndex}],'transverse mercator',['Total_Masks' filesep network_data{network_idx,network_idIndex} filesep network_data{network_idx,network_idIndex} '_MaskMap.mat'],5);
                load(['Total_Masks' filesep network_data{network_idx,network_idIndex} filesep network_data{network_idx,network_idIndex} '_MaskMap.mat']);
                disp(['[' datestr(now) '] - - ' 'Mask area successfully generated.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            
            % Retrieve the radial files to be combined and converted
            try
                station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
                station_curs = exec(conn,station_selectquery);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Query to station_tb table successfully executed.']);
            end
            
            % Fetch data
            try
                station_curs = fetch(station_curs);
                station_data = station_curs.Data;
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Data from station_tb table successfully fetched.']);
            end
            
            % Retrieve column names
            try
                station_columnNames = columnnames(station_curs,true);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
            end
            
            % Retrieve the number of stations belonging to the current network
            try
                numStations = rows(station_curs);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the current network from station_tb table successfully retrieved.']);
            end
            
            % Close cursor to station_tb table
            try
                close(station_curs);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
            end
            
            % Find the index of the input file path field
            inputPathIndexC = strfind(station_columnNames, 'radial_input_folder_path');
            inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
            
            % Find the index of the station_id field in the station_tb table
            STstation_idIndexC = strfind(station_columnNames, 'station_id');
            STstation_idIndex = find(not(cellfun('isempty', STstation_idIndexC)));
            
            % Retrieve the radial files related to the current station to be combined
            try
                toBeCombinedRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND NRT_processed_flag = 0 AND datetime>' '''' startDate ''' ORDER BY timestamp'];
                toBeCombinedRadials_curs = exec(conn,toBeCombinedRadials_selectquery);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Query to radial_input_folder_path table successfully executed.']);
            end
            
            % Fetch data
            try
                toBeCombinedRadials_curs = fetch(toBeCombinedRadials_curs);
                toBeCombinedRadials_data = toBeCombinedRadials_curs.Data;
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Data from radial_input_folder_path table successfully fetched.']);
            end
            
            % Retrieve column names
            try
                toBeCombinedRadials_columnNames = columnnames(toBeCombinedRadials_curs,true);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Column names from radial_input_folder_path table successfully retrieved.']);
            end
            
            % Retrieve the number of radials to be combined
            try
                numToBeCombinedRadials = rows(toBeCombinedRadials_curs);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Number of radials to be combined belonging to the current network from radial_input_folder_path table successfully retrieved.']);
            end
            
            % Close cursor to radial_input_folder_path table
            try
                close(toBeCombinedRadials_curs);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                HFRC_err = 1;
            end
            if(HFRC_err==0)
                disp(['[' datestr(now) '] - - ' 'Cursor to radial_input_folder_path table successfully closed.']);
            end
            
            % Find the index of the time stamp field
            timeStampIndexC = strfind(toBeCombinedRadials_columnNames, 'timestamp');
            timeStampIndex = find(not(cellfun('isempty', timeStampIndexC)));
            
            % Find the index of the filename field
            filenameIndexC = strfind(toBeCombinedRadials_columnNames, 'filename');
            filenameIndex = find(not(cellfun('isempty', filenameIndexC)));
            
            % Find the index of the station_id field in the radial_input_tb table
            RIstation_idIndexC = strfind(toBeCombinedRadials_columnNames, 'station_id');
            RIstation_idIndex = find(not(cellfun('isempty', RIstation_idIndexC)));
            
            % Find the index of the extension field
            extensionIndexC = strfind(toBeCombinedRadials_columnNames, 'extension');
            extensionIndex = find(not(cellfun('isempty', extensionIndexC)));
            
            % Find the index of the NRT_processed_flag field
            NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeCombinedRadials_columnNames, 'NRT_processed_flag'))));
            
            if(HFRC_err==0)
                % Scan the radials to be combined, group them by timestamp and combine them
                for radial_idx=1:numToBeCombinedRadials
                    if(toBeCombinedRadials_data{radial_idx,NRT_processed_flagIndex} == 0)
                        % Find the indices of the radial files of the current timestamp to be combined
                        toBeCombinedRadialIndicesC = strfind(toBeCombinedRadials_data(:,timeStampIndex), toBeCombinedRadials_data{radial_idx,timeStampIndex});
                        toBeCombinedRadialIndices = find(not(cellfun('isempty', toBeCombinedRadialIndicesC)));
                        % Build the path string according to the folder structure
                        [yMDF_err,yearFolder,monthFolder,dayFolder] = yearMonthDayFolder(toBeCombinedRadials_data{radial_idx,timeStampIndex});
                        try
                            % Build the radial file paths
                            for indices_idx=1:length(toBeCombinedRadialIndices)
                                toBeCombinedStationIndexC = strfind(station_data(:,STstation_idIndex), toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),RIstation_idIndex});
                                toBeCombinedStationIndex = find(not(cellfun('isempty', toBeCombinedStationIndexC)));
                                radFiles(indices_idx) = {[station_data{toBeCombinedStationIndex,inputPathIndex} filesep dayFolder filesep toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),filenameIndex}]};
                            end
                        catch err
                            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        try
                            % Load the radial files to be combined
                            if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                                disp(['[' datestr(now) '] - - ' 'loadRDLfile loading ...']);
                                RADIAL = loadRDLFile(radFiles, 'false', 'warning');
                            elseif(strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'crad_ascii')) % WERA data
                                % TO BE DONE
                            end
                        catch err
                            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Convert the radial files into netCDF according to the European standard data model
                        for ruv_idx=1:length(toBeCombinedRadialIndices)
                            toBeCombinedStationIndexC = strfind(station_data(:,STstation_idIndex), toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),RIstation_idIndex});
                            toBeCombinedStationIndex = find(not(cellfun('isempty', toBeCombinedStationIndexC)));
                            try
                                if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                                    [R2C_err,network_data(network_idx,:),station_data(toBeCombinedStationIndex,:),radOutputFilename,radOutputFilesize,sitePatternDate(ruv_idx,:)] = ruv2netCDF_v31(RADIAL(ruv_idx),network_data(network_idx,:),network_columnNames,station_data(toBeCombinedStationIndex,:),station_columnNames,toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),timeStampIndex});
                                    % LINES BELOW TO BE COMMENTED WHEN THE WERA FILE CONVERTER IS RUNNING
                                    disp(['[' datestr(now) '] - - ' radOutputFilename ' radial netCDF v2.1 file successfully created and stored.']);
                                    contrSitesIndices(ruv_idx) = toBeCombinedStationIndex;
                                elseif (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'crad_ascii')) % WERA data
                                    % TO BE DONE
                                end
                                %                                 % LINES BELOW TO BE UNCOMMENTED WHEN THE WERA FILE CONVERTER IS RUNNING
                                %                                 disp(['[' datestr(now) '] - - ' radOutputFilename ' radial netCDF v2.1 file successfully created and stored.']);
                                %                                 contrSitesIndices(ruv_idx) = toBeCombinedStationIndex;
                            catch err
                                display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                            % Insert radial info in radial_HFRnetCDF_tb table
                            if(HFRC_err==0)
                                try
                                    if(exist('radOutputFilename','var') ~= 0)
                                        % Delete the eventually present entry with the same name from the database
                                        radial_deletequery = ['DELETE FROM radial_HFRnetCDF_tb WHERE filename = ' '''' radOutputFilename ''''];
                                        radial_deletecurs = exec(conn,radial_deletequery);
                                        
                                        % Evaluate datetime from, Time Stamp
                                        [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),timeStampIndex});
                                        
                                        % Define a cell array containing the column names to be added
                                        addColnames = {'filename' 'network_id' 'station_id' 'timestamp' 'datetime' 'filesize' 'input_filename' 'check_flag'};
                                        
                                        % Define a cell array that contains the data for insertion
                                        addData = {radOutputFilename,network_data{network_idx,network_idIndex},station_data{toBeCombinedStationIndex,STstation_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,radOutputFilesize,toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),filenameIndex},0};
                                        
                                        % Append the product data into the radial_HFRnetCDF_tb table on the database.
                                        tablename = 'radial_HFRnetCDF_tb';
                                        datainsert(conn,tablename,addColnames,addData);
                                    end
                                catch err
                                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                    HFRC_err = 1;
                                end
                            end
                            if(HFRC_err==0)
                                disp(['[' datestr(now) '] - - ' 'Radial file information successfully inserted into radial_HFRnetCDF_tb table.']);
                            end
                        end
                        
                        % Evaluate last pattern date
                        try
                            if (HFRC_err == 0)
                                if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                                    for rd_idx=1:length(RADIAL)
                                        patternTS(rd_idx) = datenum([str2double(sitePatternDate(rd_idx, 1:5)) str2double(sitePatternDate(rd_idx, 6:8)) str2double(sitePatternDate(rd_idx, 9:11)) str2double(sitePatternDate(rd_idx, 13:15)) str2double(sitePatternDate(rd_idx, 16:18)) str2double(sitePatternDate(rd_idx, 19:20))]);
                                    end
                                    lastPatternDate = max(patternTS);
                                end
                            end
                        catch err
                            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Check if the pattern date is correct on the database and update it if needed
                        if(HFRC_err==0)
                            try
                                if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                                    last_calibration_dateIndex = find(not(cellfun('isempty', strfind(network_columnNames, 'last_calibration_date'))));
                                    if(datenum(network_data{network_idx,last_calibration_dateIndex})<lastPatternDate)
                                        network_tbUpdateFlag = 1;
                                        network_data{network_idx,last_calibration_dateIndex} = datestr(lastPatternDate, 'yyyy-mm-dd');
                                    end
                                end
                            catch err
                                display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                        end
                        
                        % Combine the Codar radial files into total
                        if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                            try
                                disp(['[' datestr(now) '] - - ' 'makeTotals combining ...']);
                                [TUV,R] = makeTotals(RADIAL, 'Grid', Grid, 'TimeStamp', RADIAL(1,1).TimeStamp, 'spatthresh', network_data{network_idx,spatthreshIndex}, 'tempthresh', 1/24);
                                % Totals setting on a regular grid
                                [TUVgrid,DIM,I] = gridTotals( TUV, 'true', 'true');
                                % Totals masking
                                [TUVmask,I] = maskTotals(TUVgrid,ncst,0);
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                            
                            % Save the total mat file
                            if(HFRC_err==0)
                                try
                                    ts = datevec(TUVmask.TimeStamp);
                                    time_str = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
                                    [tFB_err, matFilePath] = totalFolderBuilder_v21(network_data{network_idx,matPathIndex}, toBeCombinedRadials_data{radial_idx,timeStampIndex});
                                    save([matFilePath filesep network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat'], 'TUVmask');
                                    disp(['[' datestr(now) '] - - ' network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat' ' file successfully saved.']);
                                catch err
                                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                    HFRC_err = 1;
                                end
                            end
                        end
                        
                        % Create the total netCDF file according to the European standard data model
                        if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                            [T2C_err,network_data(network_idx,:),station_data(contrSitesIndices,:),totOutputFilename,totOutputFilesize] = tot2netCDF_v31(TUVmask,network_data(network_idx,:),network_columnNames,station_data(contrSitesIndices,:),station_columnNames,toBeCombinedRadials_data{radial_idx,timeStampIndex});
                            % LINE BELOW TO BE COMMENTED WHEN THE WERA FILE CONCERTER IS RUNNING
                            disp(['[' datestr(now) '] - - ' totOutputFilename ' total netCDF v2.1 file successfully created and stored.']);
                        elseif (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'crad_ascii')) % WERA data
                            % TO BE DONE
                        end
                        %                         % LINE BELOW TO BE UNCOMMENTED WHEN THE WERA FILE CONCERTER IS RUNNING
                        %                         disp(['[' datestr(now) '] - - ' totOutputFilename ' total netCDF v2.1 file successfully created and stored.']);
                        
                        % Update NRT_processed_flag in the local radial table
                        if(HFRC_err==0)
                            try
                                toBeCombinedRadials_data(toBeCombinedRadialIndices,NRT_processed_flagIndex)={1};
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                        end
                        
                        % Update NRT_processed_flag in radial_input_tb table
                        if(HFRC_err==0)
                            try
                                if(length(toBeCombinedRadialIndices) == numStations)
                                    % Define a cell array containing the column names to be updated
                                    updateColnames = {'NRT_processed_flag'};
                                    
                                    % Define a cell array that contains the data for insertion
                                    updateData = {1};
                                    
                                    % Update the radial_input_tb table on the database
                                    tablename = 'radial_input_tb';
                                    whereclause = ['WHERE timestamp = ' '''' toBeCombinedRadials_data{radial_idx,timeStampIndex} ''''];
                                    update(conn,tablename,updateColnames,updateData,whereclause);
                                end
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                            if(HFRC_err==0)
                                disp(['[' datestr(now) '] - - ' 'radial_input_tb table successfully updated with NRT processed flag.']);
                            end
                        end
                        
                        % Insert converted total info in total_HFRnetCDF_tb table
                        if(HFRC_err==0)
                            try
                                if(exist('totOutputFilename','var') ~= 0)
                                    % Delete the eventually present entry with the same name from the database
                                    total_deletequery = ['DELETE FROM total_HFRnetCDF_tb WHERE filename = ' '''' totOutputFilename ''''];
                                    total_deletecurs = exec(conn,total_deletequery);
                                    
                                    % Evaluate datetime from, Time Stamp
                                    [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{radial_idx,timeStampIndex});
                                    
                                    % Define a cell array containing the column names to be added
                                    addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'filesize' 'mat_filename' 'check_flag'};
                                    
                                    % Define a cell array that contains the data for insertion
                                    addData = {totOutputFilename,network_data{network_idx,network_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,totOutputFilesize,[network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat'],0};
                                    
                                    % Append the product data into the total_HFRnetCDF_tb table on the database.
                                    tablename = 'total_HFRnetCDF_tb';
                                    datainsert(conn,tablename,addColnames,addData);
                                end
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                        end
                        if(HFRC_err==0)
                            disp(['[' datestr(now) '] - - ' 'Total combined file information successfully inserted into total_HFRnetCDF_tb table.']);
                        end
                        
                    end
                    
                end
                
            end
            
            % Update the last calibration date on the network_tb table of the database (if needed)
            if(HFRC_err==0)
                try
                    if(network_tbUpdateFlag == 1)
                        % Define a cell array containing the column names to be updated
                        updateColnames = {'last_calibration_date'};
                        
                        % Define a cell array that contains the data for insertion
                        updateData = network_data(network_idx,last_calibration_dateIndex);
                        
                        % Update the network_tb table on the database
                        tablename = 'network_tb';
                        whereclause = ['WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
                        update(conn,tablename,updateColnames,updateData,whereclause);
                        disp(['[' datestr(now) '] - - ' 'network_tb table successfully updated with last calibration date.']);
                    end
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    HFRC_err = 1;
                end
            end
            
            clear Grid gridLon gridLat lonG latG lon lat ruvFiles;
            
        end
        
    end
    
    %%
    
    %% Close connection
    
    try
        close(conn);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    if(HFRC_err==0)
        disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
    end
    
    %%
    
    pause(2700);
    
end