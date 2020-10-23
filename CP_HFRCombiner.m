%% CP_HFRCombiner.m
% This function reads the HFR database for collecting information about
% the radial data files to be combined into totals and performs the
% combination and the generation of radial and total data files into the
% European standard data model.

% INPUT:
%         conn: connection to database
%         startDate: processing start date
%         networkData: cell array containing information about the processed
%                      network (metadata)
%         networkFields: field names of the cell array containing
%                       information about the processed network.
%         stationData: cell array containing information about the stations
%                      of the processed network (metadata)
%         stationFields: field names of the cell array containing
%                       information about the stations of the processed network.

% OUTPUT:
%         HFRC_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: February 10, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [HFRC_err] = CP_HFRCombiner(conn,startDate,networkData,networkFields,stationData,stationFields)

warning('off', 'all');

HFRC_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_HFRCombiner started.']);

%%

%% Scan the input radial files, convert them and combine them into totals

try
    % Find the index of the network_id field
    network_idIndexC = strfind(networkFields, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the indices of the geospatial boundaries fields
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    
    % Find the index of the grid resolution field
    grid_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'grid_resolution'))));
    
    % Find the index of the of the combination search radius field
    matPathIndex = find(not(cellfun('isempty', strfind(networkFields, 'total_mat_folder_path'))));
    
    % Find the index of the of the output mat files folder path field
    spatthreshIndex = find(not(cellfun('isempty', strfind(networkFields, 'combination_search_radius'))));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

try
    HFRC_err = 0;
    try
        % Build the regular LonLat grid given the geographical boundaries and the grid resolution for the radial combination into total
        [gridLon, gridLat] = LonLat_grid([networkData{1,geospatial_lon_minIndex},networkData{1,geospatial_lat_minIndex}], [networkData{1,geospatial_lon_maxIndex},networkData{1,geospatial_lat_maxIndex}], networkData{1,grid_resolutionIndex}, 'km');
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
        disp(['[' datestr(now) '] - - ' 'Grid for radial combination for ' networkData{1, network_idIndex} ' successfully generated.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    % Build the mask for the total masking
    try
        if (exist(['Total_Masks' filesep networkData{1,network_idIndex}], 'dir') ~= 7)
            mkdir(['Total_Masks' filesep networkData{1,network_idIndex}]);
        end
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    try
        fname = makeCoast([networkData{1,geospatial_lon_minIndex},networkData{1,geospatial_lon_maxIndex}],[networkData{1,geospatial_lat_minIndex},networkData{1,geospatial_lat_maxIndex}],'transverse mercator',['Total_Masks' filesep networkData{1,network_idIndex} filesep networkData{1,network_idIndex} '_MaskMap.mat'],5);
        load(['Total_Masks' filesep networkData{1,network_idIndex} filesep networkData{1,network_idIndex} '_MaskMap.mat']);
        disp(['[' datestr(now) '] - - ' 'Mask area for ' networkData{1, network_idIndex} ' successfully generated.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    % Retrieve the radial files to be combined and converted
    try
        % Find the index of the input file path field
        inputPathIndexC = strfind(stationFields, 'radial_input_folder_path');
        inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
        
        % Find the index of the station_id field in the station_tb table
        STstation_idIndexC = strfind(stationFields, 'station_id');
        STstation_idIndex = find(not(cellfun('isempty', STstation_idIndexC)));
        
        % Find the index of the last calibration date field
        last_calibration_dateIndexC = strfind(stationFields, 'last_calibration_date');
        last_calibration_dateIndex = find(not(cellfun('isempty', last_calibration_dateIndexC)));
        
        % Find the index of the end of operation date field
        operational_toIndexC = strfind(stationFields, 'operational_to');
        operational_toIndex = find(not(cellfun('isempty', operational_toIndexC)));
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    % Retrieve the number of operational stations
    try
        numStations = size(stationData,1);
        numActiveStations = numStations;
        for station_idx=1:numStations
            if(size(stationData{station_idx,operational_toIndex},2)~=4)
                numActiveStations = numActiveStations - 1;
            end
        end
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    % Retrieve the radial files related to the current station to be combined
    try
        % Manage the case of the ISMAR-LaMMA integrated network (HFR-WesternItaly)
        if(strcmp(networkData{1,network_idIndex},'HFR-WesternItaly'))
            toBeCombinedRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND (network_id = ' '''HFR-TirLig'' OR network_id = ' '''HFR-LaMMA'') AND NRT_processed_flag_integrated_network = 0 ORDER BY timestamp'];
        else
            toBeCombinedRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' networkData{1,network_idIndex} ''' AND NRT_processed_flag = 0 ORDER BY timestamp'];
        end
        toBeCombinedRadials_curs = exec(conn,toBeCombinedRadials_selectquery);
        disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table for retrieving the radial files from ' networkData{1,network_idIndex} ' network to be combined successfully executed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    % Fetch data
    try
        toBeCombinedRadials_curs = fetch(toBeCombinedRadials_curs);
        toBeCombinedRadials_data = toBeCombinedRadials_curs.Data;
        disp(['[' datestr(now) '] - - ' 'Data of the radial files from ' networkData{1,network_idIndex} ' network to be combined successfully fetched from radial_input_tb table.']);
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
        disp(['[' datestr(now) '] - - ' 'Number of the radial files from ' networkData{1,network_idIndex} ' network to be combined successfully retrieved from radial_input_tb table.']);
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
        if(strcmp(networkData{1,network_idIndex},'HFR-WesternItaly'))
            NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeCombinedRadials_columnNames, 'NRT_processed_flag_integrated_network'))));
        else
            NRT_processed_flagIndex = strmatch('NRT_processed_flag',toBeCombinedRadials_columnNames,'exact');
            %                 NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeCombinedRadials_columnNames, 'NRT_processed_flag'))));
        end
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    try
        % Scan the radials to be combined, group them by timestamp and combine them
        for radial_idx=1:numToBeCombinedRadials
            HFRC_err = 0;
            if(toBeCombinedRadials_data{radial_idx,NRT_processed_flagIndex} == 0)
                % Find the indices of the radial files of the current timestamp to be combined
                toBeCombinedRadialIndicesC = strfind(toBeCombinedRadials_data(:,timeStampIndex), toBeCombinedRadials_data{radial_idx,timeStampIndex});
                toBeCombinedRadialIndices = find(not(cellfun('isempty', toBeCombinedRadialIndicesC)));
                try
                    % Build the radial file paths
                    for indices_idx=1:length(toBeCombinedRadialIndices)
                        toBeCombinedStationIndexC = strfind(stationData(:,STstation_idIndex), toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),RIstation_idIndex});
                        toBeCombinedStationIndex = find(not(cellfun('isempty', toBeCombinedStationIndexC)));
                        stationData{toBeCombinedStationIndex,inputPathIndex} = strtrim(stationData{toBeCombinedStationIndex,inputPathIndex});
                        radFiles(indices_idx) = {[toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),filepathIndex} filesep toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),filenameIndex}]};
                    end
                catch err
                    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    HFRC_err = 1;
                end
                
                try
                    % Load the radial files to be combined
                    if (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.ruv')) % Codar data
                        disp(['[' datestr(now) '] - - ' 'loadRDLfile loading ...']);
                        RADIAL = loadRDLFile(radFiles, 'false', 'warning');
%                    elseif(strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.mat')) % MetNo data
%                        for mat_idx=1:length(radFiles)
%                            load(radFiles{mat_idx});
%                            RADIAL(mat_idx) = ruvRad;
%                            clear ruvRad
%                        end
                    elseif(strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.crad_ascii')) % WERA data
                        % NOTHING TO DO
                    end
                catch err
                    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    HFRC_err = 1;
                end
                
                % Convert the radial files into netCDF according to the European standard data model
                for ruv_idx=1:length(toBeCombinedRadialIndices)
                    station_tbUpdateFlag = 0;
                    toBeCombinedStationIndexC = strfind(stationData(:,STstation_idIndex), toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),RIstation_idIndex});
                    toBeCombinedStationIndex = find(not(cellfun('isempty', toBeCombinedStationIndexC)));
                    try
                        if (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.ruv')) % Codar data
                            if(~strcmp(networkData{1,network_idIndex},'HFR-WesternItaly'))
                                % v2.2
                                [R2C_err,networkData(1,:),stationData(toBeCombinedStationIndex,:),radOutputFilename,radOutputFilesize,station_tbUpdateFlag] = ruv2netCDF_v22(RADIAL(ruv_idx),networkData(1,:),networkFields,stationData(toBeCombinedStationIndex,:),stationFields,toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),timeStampIndex});
                                disp(['[' datestr(now) '] - - ' radOutputFilename ' radial netCDF v2.2 file successfully created and stored.']);
                            else
                                station_tbUpdateFlag = 0;
                            end
%                        elseif (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.mat')) % MetNo data
%                            % v2.2
%                            [R2C_err,networkData(1,:),stationData(toBeCombinedStationIndex,:),radOutputFilename,radOutputFilesize,station_tbUpdateFlag] = mat2netCDF_v33(RADIAL(ruv_idx),networkData(1,:),networkFields,stationData(toBeCombinedStationIndex,:),stationFields,toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),timeStampIndex});
                        elseif (strcmp(toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),extensionIndex}, '.crad_ascii')) % WERA data
                            [R2C_err,networkData(1,:),radOutputFilename,radOutputFilesize] = cradAscii2netCDF_v22(radFiles{ruv_idx},networkData(1,:),networkFields,stationData(toBeCombinedStationIndex,:),stationFields,toBeCombinedRadials_data{toBeCombinedRadialIndices(indices_idx),timeStampIndex});
                            numActiveStations = length(toBeCombinedRadialIndices); % WERA radials are not combined
                            disp(['[' datestr(now) '] - - ' radOutputFilename ' radial netCDF v2.2 file successfully created and stored.']);
                        end
                        contrSitesIndices(ruv_idx) = toBeCombinedStationIndex;
                    catch err
                        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        HFRC_err = 2;
                    end
                    
                    % Insert radial info in radial_HFRnetCDF_tb table
                    try
                        if((HFRC_err==0) && (exist('radOutputFilename','var') ~= 0))
                            % Delete the eventually present entry with the same name from the database
                            radial_deletequery = ['DELETE FROM radial_HFRnetCDF_tb WHERE filename = ' '''' radOutputFilename ''''];
                            radial_deletecurs = exec(conn,radial_deletequery);
                            
                            % Evaluate datetime from, Time Stamp
                            [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),timeStampIndex});
                            
                            % Define a cell array containing the column names to be added
                            addColnames = {'filename' 'network_id' 'station_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'input_filename' 'check_flag'};
                            
                            % Define a cell array that contains the data for insertion
                            addData = {radOutputFilename,networkData{1,network_idIndex},stationData{toBeCombinedStationIndex,STstation_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),radOutputFilesize,toBeCombinedRadials_data{toBeCombinedRadialIndices(ruv_idx),filenameIndex},0};
                            
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
                        if((HFRC_err==0) && (station_tbUpdateFlag == 1))
                            % Define a cell array containing the column names to be updated
                            updateColnames = {'last_calibration_date'};
                            
                            % Define a cell array that contains the data for insertion
                            updateData = stationData(toBeCombinedStationIndex,last_calibration_dateIndex);
                            
                            % Update the network_tb table on the database
                            tablename = 'station_tb';
                            whereclause = ['WHERE station_id = ' '''' stationData{toBeCombinedStationIndex,STstation_idIndex} ''''];
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
                if((size(radFiles,2)>1) && (HFRC_err~=2))
                    % Check if all the radials have the same extension
                    for ext_idx=1:length(toBeCombinedRadialIndices)
                        extensions{ext_idx} = toBeCombinedRadials_data{toBeCombinedRadialIndices(ext_idx),extensionIndex};
                    end
                    extensions = uniqueStrCell(extensions);
                    % Check the extension
                    if ((length(extensions)==1) && ((strcmp(extensions, '.ruv')) || (strcmp(extensions, '.mat')))) % Codar data
                        try
                            disp(['[' datestr(now) '] - - ' 'makeTotals combining radials...']);
                            [TUV,R] = makeTotals(RADIAL, 'Grid', Grid, 'TimeStamp', RADIAL(1,1).TimeStamp, 'spatthresh', networkData{1,spatthreshIndex}, 'tempthresh', 1/24);
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
                            networkData{1,matPathIndex} = strtrim(networkData{1,matPathIndex});
                            % v2.2
                            [tFB_err, matFilePath] = totalFolderBuilder_v22(networkData{1,matPathIndex}, toBeCombinedRadials_data{radial_idx,timeStampIndex});
                            save([matFilePath filesep networkData{1,network_idIndex} '-Total_' time_str '.mat'], 'TUVmask');
                            disp(['[' datestr(now) '] - - ' networkData{1,network_idIndex} '-Total_' time_str '.mat' ' file successfully saved.']);
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Plot the current map for HFR-TirLig network
                        try
                            if(strcmp(networkData{1,network_idIndex},'HFR-TirLig'))
                                % Totals cleaning for GDOP
                                gdop_sP = sqrt(6.25);
                                maxspd_sP = 500;
                                [TUVclean,I] = cleanTotals(TUVmask,maxspd_sP,{'GDOPMaxOrthog','TotalErrors',gdop_sP});
                                % Plot
                                shadePlot_TirLig;
                                % Save the map file
                                saveas(gcf,['/home/radarcombine/EU_HFR_NODE/HFR_TirLig/Totals_map/' time_str '.jpg']);
                                close
                                disp(['[' datestr(now) '] - - ' time_str ' map for HFR-TirLig network successfully saved.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                        
                        % Plot the current map for HFR-WesternItaly network
                        try
                            if(strcmp(networkData{1,network_idIndex},'HFR-WesternItaly'))
                                % Totals cleaning for GDOP
                                gdop_sP = sqrt(6.25);
                                maxspd_sP = 500;
                                [TUVclean,I] = cleanTotals(TUVmask,maxspd_sP,{'GDOPMaxOrthog','TotalErrors',gdop_sP});
                                % Plot
                                shadePlot_WesternItaly;
                                % Save the map file
                                saveas(gcf,['/home/radarcombine/EU_HFR_NODE/HFR_WesternItaly/Totals_map/' time_str '.jpg']);
                                close
                                disp(['[' datestr(now) '] - - ' time_str ' map for HFR-TirLig network successfully saved.']);
                            end
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            HFRC_err = 1;
                        end
                    end
                    
                    % Create the total netCDF file according to the European standard data model
                    try
                        if (strcmp(extensions, '.ruv') || strcmp(extensions, '.mat')) % Codar data
                            % v2.2
                            [T2C_err,networkData(1,:),stationData(contrSitesIndices,:),totOutputFilename,totOutputFilesize] = tot2netCDF_v22(TUVmask,networkData(1,:),networkFields,stationData(contrSitesIndices,:),stationFields,toBeCombinedRadials_data{radial_idx,timeStampIndex},stationData);
                            disp(['[' datestr(now) '] - - ' totOutputFilename ' total netCDF v2.2 file successfully created and stored.']);
                        elseif (strcmp(toBeCombinedRadials_data{toBeCombinedStationIndex,extensionIndex}, 'crad_ascii')) % WERA data
                            % NOTHING TO DO
                        end
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        HFRC_err = 1;
                    end
                    
                    % Insert converted total info in total_HFRnetCDF_tb table
                    try
                        if((HFRC_err==0) && (exist('totOutputFilename','var') ~= 0))
                            % Delete the eventually present entry with the same name from the database
                            total_deletequery = ['DELETE FROM total_HFRnetCDF_tb WHERE filename = ' '''' totOutputFilename ''''];
                            total_deletecurs = exec(conn,total_deletequery);
                            
                            % Evaluate datetime from, Time Stamp
                            [t2d_err,DateTime] = timestamp2datetime(toBeCombinedRadials_data{radial_idx,timeStampIndex});
                            
                            % Define a cell array containing the column names to be added
                            addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'mat_filename' 'check_flag'};
                            
                            % Define a cell array that contains the data for insertion
                            addData = {totOutputFilename,networkData{1,network_idIndex},toBeCombinedRadials_data{radial_idx,timeStampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),totOutputFilesize,[networkData{1,network_idIndex} '-Total_' time_str '.mat'],0};
                            
                            % Append the product data into the total_HFRnetCDF_tb table on the database.
                            tablename = 'total_HFRnetCDF_tb';
                            datainsert(conn,tablename,addColnames,addData);
                            disp(['[' datestr(now) '] - - ' totOutputFilename ' total combined file information successfully inserted into total_HFRnetCDF_tb table.']);
                        end
                    catch err
                        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                        HFRC_err = 1;
                    end
                end
                
                % Update NRT_processed_flag in the local radial table
                try
                    if(HFRC_err==0)
                        toBeCombinedRadials_data(toBeCombinedRadialIndices,NRT_processed_flagIndex)={1};
                    end
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    HFRC_err = 1;
                end
                
                % Update NRT_processed_flag in radial_input_tb table
                try
                    if((HFRC_err==0) && (length(toBeCombinedRadialIndices) == numActiveStations))
                        % Define a cell array containing the column names to be updated
                        if(strcmp(networkData{1,network_idIndex},'HFR-WesternItaly'))
                            updateColnames = {'NRT_processed_flag_integrated_network'};
                            whereclause = ['WHERE timestamp = ' '''' toBeCombinedRadials_data{radial_idx,timeStampIndex} ''' AND (network_id = ' '''HFR-TirLig'' OR network_id = ' '''HFR-LaMMA'')'];
                        else
                            updateColnames = {'NRT_processed_flag'};
                            whereclause = ['WHERE timestamp = ' '''' toBeCombinedRadials_data{radial_idx,timeStampIndex} ''' AND network_id = ' '''' networkData{1,network_idIndex} ''''];
                        end
                        
                        % Define a cell array that contains the data for insertion
                        updateData = {1};
                        
                        % Update the radial_input_tb table on the database
                        tablename = 'radial_input_tb';
                        update(conn,tablename,updateColnames,updateData,whereclause);
                        disp(['[' datestr(now) '] - - ' 'radial_input_tb table successfully updated with NRT processed flag for timestamp ' toBeCombinedRadials_data{radial_idx,timeStampIndex} '.']);
                    end
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    HFRC_err = 1;
                end
                
                clear radFiles RADIAL contrSitesIndices TUV TUVgrid TUVmask radOutputFilename totOutputFilename;
                
            end
        end
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        HFRC_err = 1;
    end
    
    %     clear Grid gridLon gridLat lonG latG lon lat;
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    HFRC_err = 1;
end

%%

if(HFRC_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_HFRCombiner successfully executed.']);
end

% pause(20);

return
