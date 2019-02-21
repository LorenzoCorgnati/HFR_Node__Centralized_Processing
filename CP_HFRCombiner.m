%% CP_HFRCombiner.m
% This application reads the HFR database for collecting information about
% the radial data files to be combined into totals and performs the
% combination and the generation of radial and total data files into the
% European standard data model.

% Author: Lorenzo Corgnati
% Date: July 17, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

HFRC_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_HFRCombiner started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
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
    HFRC_err = 1;
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
    disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

% Close cursor
try
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

%%

%% Scan the networks and combine the related radial files

try
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
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

% Scan the networks
try
    for network_idx=1:numNetworks
        try
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
            disp(['[' datestr(now) '] - - ' 'Grid for radial combination for ' network_data{network_idx, network_idIndex} ' successfully generated.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
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
            disp(['[' datestr(now) '] - - ' 'Mask area for ' network_data{network_idx, network_idIndex} ' successfully generated.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        
        % Retrieve the radial files to be combined and converted
        try
            station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''''];
            station_curs = exec(conn,station_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of the ' network_data{network_idx,network_idIndex} ' network successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Fetch data
        try
            station_curs = fetch(station_curs);
            station_data = station_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the stations of the ' network_data{network_idx,network_idIndex} ' network successfully fetched from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Retrieve column names
        try
            station_columnNames = columnnames(station_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Retrieve the number of stations belonging to the current network
        try
            numStations = rows(station_curs);
            disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the ' network_data{network_idx,network_idIndex} ' network successfully retrieved from station_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Close cursor to station_tb table
        try
            close(station_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        try
            % Find the index of the input file path field
            inputPathIndexC = strfind(station_columnNames, 'radial_input_folder_path');
            inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
            
            % Find the index of the station_id field in the station_tb table
            STstation_idIndexC = strfind(station_columnNames, 'station_id');
            STstation_idIndex = find(not(cellfun('isempty', STstation_idIndexC)));
            
            % Find the index of the last calibration date field
            last_calibration_dateIndexC = strfind(station_columnNames, 'last_calibration_date');
            last_calibration_dateIndex = find(not(cellfun('isempty', last_calibration_dateIndexC)));
            
            % Find the index of the end of operation date field
            operational_toIndexC = strfind(station_columnNames, 'operational_to');
            operational_toIndex = find(not(cellfun('isempty', operational_toIndexC)));
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        
        % Retrieve the number of operational stations
        try
            numActiveStations = numStations;
            for station_idx=1:numStations
                if(size(station_data{station_idx,operational_toIndex},2)~=4)
                    numActiveStations = numActiveStations - 1;
                end
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        
        % Retrieve the radial files related to the current station to be combined
        try
            toBeCombinedRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND NRT_processed_flag = 0 ORDER BY timestamp'];
            toBeCombinedRadials_curs = exec(conn,toBeCombinedRadials_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table for retrieving the radial files from ' network_data{network_idx,network_idIndex} ' network to be combined successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Fetch data
        try
            toBeCombinedRadials_curs = fetch(toBeCombinedRadials_curs);
            toBeCombinedRadials_data = toBeCombinedRadials_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the radial files from ' network_data{network_idx,network_idIndex} ' network to be combined successfully fetched from radial_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Retrieve column names
        try
            toBeCombinedRadials_columnNames = columnnames(toBeCombinedRadials_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from radial_input_folder_path table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Retrieve the number of radials to be combined
        try
            numToBeCombinedRadials = rows(toBeCombinedRadials_curs);
            disp(['[' datestr(now) '] - - ' 'Number of the radial files from ' network_data{network_idx,network_idIndex} ' network to be combined successfully retrieved from radial_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        % Close cursor to radial_input_folder_path table
        try
            close(toBeCombinedRadials_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to radial_input_folder_path table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
                
        try
            % Find the index of the time stamp field
            timeStampIndexC = strfind(toBeCombinedRadials_columnNames, 'timestamp');
            timeStampIndex = find(not(cellfun('isempty', timeStampIndexC)));
            
            % Find the index of the filename field
            filenameIndexC = strfind(toBeCombinedRadials_columnNames, 'filename');
            filenameIndex = find(not(cellfun('isempty', filenameIndexC)));
            
            % Find the index of the filepath field
            filepathIndexC = strfind(toBeCombinedRadials_columnNames, 'filepath');
            filepathIndex = find(not(cellfun('isempty', filepathIndexC)));
            
            % Find the index of the station_id field in the radial_input_tb table
            RIstation_idIndexC = strfind(toBeCombinedRadials_columnNames, 'station_id');
            RIstation_idIndex = find(not(cellfun('isempty', RIstation_idIndexC)));
            
            % Find the index of the extension field
            extensionIndexC = strfind(toBeCombinedRadials_columnNames, 'extension');
            extensionIndex = find(not(cellfun('isempty', extensionIndexC)));
            
            % Find the index of the NRT_processed_flag field
            NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeCombinedRadials_columnNames, 'NRT_processed_flag'))));                        
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        
        try
            % Scan the radials to be combined, group them by timestamp and combine them
            for radial_idx=1:numToBeCombinedRadials
                if(toBeCombinedRadials_data{radial_idx,NRT_processed_flagIndex} == 0)
                    % Find the indices of the radial files of the current timestamp to be combined
                    toBeCombinedRadialIndicesC = strfind(toBeCombinedRadials_data(:,timeStampIndex), toBeCombinedRadials_data{radial_idx,timeStampIndex});
                    toBeCombinedRadialIndices = find(not(cellfun('isempty', toBeCombinedRadialIndicesC)));                                       
                    try
                        % Build the radial file paths
                        for indices_idx=1:length(toBeCombinedRadialIndices)
                            toBeCombinedStationIndexC = strfind(station_data(:,STstation_idIndex), toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),RIstation_idIndex});
                            toBeCombinedStationIndex = find(not(cellfun('isempty', toBeCombinedStationIndexC)));
                            station_data{toBeCombinedStationIndex,inputPathIndex} = strtrim(station_data{toBeCombinedStationIndex,inputPathIndex});
                            radFiles(indices_idx) = {[toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),filepathIndex} filesep toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),filenameIndex}]};                                                        
                        end
                    catch err
                        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        HFRC_err = 1;
                    end
                    
                    try
                        % Load the radial files to be combined
                        if (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, 'ruv')) % Codar data
                            disp(['[' datestr(now) '] - - ' 'loadRDLfile loading ...']);
                            RADIAL = loadRDLFile(radFiles, 'false', 'warning');
                        elseif(strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, 'crad_ascii')) % WERA data
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
                            if (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, 'ruv')) % Codar data
                                [R2C_err,network_data(network_idx,:),station_data(toBeCombinedStationIndex,:),radOutputFilename,radOutputFilesize,station_tbUpdateFlag] = ruv2netCDF_v31(RADIAL(ruv_idx),network_data(network_idx,:),network_columnNames,station_data(toBeCombinedStationIndex,:),station_columnNames,toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),timeStampIndex});
                                % LINES BELOW TO BE COMMENTED WHEN THE WERA FILE CONVERTER IS RUNNING
                                disp(['[' datestr(now) '] - - ' radOutputFilename ' radial netCDF v2.1 file successfully created and stored.']);
                                contrSitesIndices(ruv_idx) = toBeCombinedStationIndex;
                            elseif (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, 'crad_ascii')) % WERA data
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
                        try
                            if(exist('radOutputFilename','var') ~= 0)
                                % Delete the eventually present entry with the same name from the database
                                radial_deletequery = ['DELETE FROM radial_HFRnetCDF_tb WHERE filename = ' '''' radOutputFilename ''''];
                                radial_deletecurs = exec(conn,radial_deletequery);
                                
                                % Evaluate datetime from, Time Stamp
                                [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),timeStampIndex});
                                
                                % Define a cell array containing the column names to be added
                                addColnames = {'filename' 'network_id' 'station_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'input_filename' 'check_flag'};
                                
                                % Define a cell array that contains the data for insertion
                                addData = {radOutputFilename,network_data{network_idx,network_idIndex},station_data{toBeCombinedStationIndex,STstation_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),radOutputFilesize,toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),filenameIndex},0};
                                
                                % Append the product data into the radial_HFRnetCDF_tb table on the database.
                                tablename = 'radial_HFRnetCDF_tb';
                                datainsert(conn,tablename,addColnames,addData);
                                disp(['[' datestr(now) '] - - ' radOutputFilename ' radial file information successfully inserted into radial_HFRnetCDF_tb table.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                                                
                        % Update the last calibration date on the station_tb table of the database (if needed)
                        try
                            if(station_tbUpdateFlag == 1)
                                % Define a cell array containing the column names to be updated
                                updateColnames = {'last_calibration_date'};
                                
                                % Define a cell array that contains the data for insertion
                                updateData = station_data(toBeCombinedStationIndex,last_calibration_dateIndex);
                                
                                % Update the network_tb table on the database
                                tablename = 'station_tb';
                                whereclause = ['WHERE station_id = ' '''' station_data{toBeCombinedStationIndex,STstation_idIndex} ''''];
                                update(conn,tablename,updateColnames,updateData,whereclause);
                                station_tbUpdateFlag = 0;
                                disp(['[' datestr(now) '] - - ' 'station_tb table successfully updated with last calibration date.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                    end
                    
                    % Combine the Codar radial files into total
                    if(size(radFiles,2)>1)
                        if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                            try
                                disp(['[' datestr(now) '] - - ' 'makeTotals combining radials...']);
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
                            try
                                ts = datevec(TUVmask.TimeStamp);
                                time_str = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
                                network_data{network_idx,matPathIndex} = strtrim(network_data{network_idx,matPathIndex});
                                [tFB_err, matFilePath] = totalFolderBuilder_v21(network_data{network_idx,matPathIndex}, toBeCombinedRadials_data{radial_idx,timeStampIndex});
                                save([matFilePath filesep network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat'], 'TUVmask');
                                disp(['[' datestr(now) '] - - ' network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat' ' file successfully saved.']);
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                HFRC_err = 1;
                            end
                        end
                        
                        % Create the total netCDF file according to the European standard data model
                        try
                            if (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'ruv')) % Codar data
                                [T2C_err,network_data(network_idx,:),station_data(contrSitesIndices,:),totOutputFilename,totOutputFilesize] = tot2netCDF_v31(TUVmask,network_data(network_idx,:),network_columnNames,station_data(contrSitesIndices,:),station_columnNames,toBeCombinedRadials_data{radial_idx,timeStampIndex});
                                % LINE BELOW TO BE COMMENTED WHEN THE WERA FILE CONCERTER IS RUNNING
                                disp(['[' datestr(now) '] - - ' totOutputFilename ' total netCDF v2.1 file successfully created and stored.']);
                            elseif (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'crad_ascii')) % WERA data
                                % TO BE DONE
                            end
                            %                         % LINE BELOW TO BE UNCOMMENTED WHEN THE WERA FILE CONCERTER IS RUNNING
                            %                         disp(['[' datestr(now) '] - - ' totOutputFilename ' total netCDF v2.1 file successfully created and stored.']);
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Update NRT_processed_flag in the local radial table
                        try
                            toBeCombinedRadials_data(toBeCombinedRadialIndices,NRT_processed_flagIndex)={1};
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Update NRT_processed_flag in radial_input_tb table
                        try
                            if(length(toBeCombinedRadialIndices) == numActiveStations)
                                % Define a cell array containing the column names to be updated
                                updateColnames = {'NRT_processed_flag'};
                                
                                % Define a cell array that contains the data for insertion
                                updateData = {1};
                                
                                % Update the radial_input_tb table on the database
                                tablename = 'radial_input_tb';
                                whereclause = ['WHERE timestamp = ' '''' toBeCombinedRadials_data{radial_idx,timeStampIndex} ''''];
                                update(conn,tablename,updateColnames,updateData,whereclause);
                                disp(['[' datestr(now) '] - - ' 'radial_input_tb table successfully updated with NRT processed flag for timestamp ' toBeCombinedRadials_data{radial_idx,timeStampIndex} '.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                                                
                        % Insert converted total info in total_HFRnetCDF_tb table
                        try
                            if(exist('totOutputFilename','var') ~= 0)
                                % Delete the eventually present entry with the same name from the database
                                total_deletequery = ['DELETE FROM total_HFRnetCDF_tb WHERE filename = ' '''' totOutputFilename ''''];
                                total_deletecurs = exec(conn,total_deletequery);
                                
                                % Evaluate datetime from, Time Stamp
                                [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{radial_idx,timeStampIndex});
                                
                                % Define a cell array containing the column names to be added
                                addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'mat_filename' 'check_flag'};
                                
                                % Define a cell array that contains the data for insertion
                                addData = {totOutputFilename,network_data{network_idx,network_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),totOutputFilesize,[network_data{network_idx,network_idIndex} '_TOTL_' time_str '.mat'],0};
                                
                                % Append the product data into the total_HFRnetCDF_tb table on the database.
                                tablename = 'total_HFRnetCDF_tb';
                                datainsert(conn,tablename,addColnames,addData);
                                disp(['[' datestr(now) '] - - ' totOutputFilename 'total combined file information successfully inserted into total_HFRnetCDF_tb table.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end                        
                    end
                    
                    clear radFiles RADIAL contrSitesIndices TUV TUVgrid TUVmask radOutputFilename totOutputFilename;
                    
                end
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            HFRC_err = 1;
        end
        
        clear Grid gridLon gridLat lonG latG lon lat;
        
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

%%

%% Close connection

try
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

%%

if(HFRC_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_HFRCombiner successfully executed.']);
end

pause(1200);