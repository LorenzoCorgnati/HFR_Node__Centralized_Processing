%% CP_TotalConversion_main.m
% This function reads the HFR database for collecting information about
% the total data files to be converted into the European standard data
% model and calls the conversion functions.

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
%         TC_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: February 11, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [TC_err] = CP_TotalConversion(conn,startDate,networkData,networkFields,stationData,stationFields)

warning('off', 'all');

TC_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_TotalConversion started.']);

%%

%% Scan the input total files and convert them

try
    % Find the index of the network_id field
    network_idIndexC = strfind(networkFields, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the index of the input file path field
    inputPathIndexC = strfind(networkFields, 'total_input_folder_path');
    inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

% Scan the networks
try
    TC_err = 0;
    % Process HFR-MATROOS network
    if(strcmp(networkData{1,network_idIndex},'HFR-MATROOS'))
        % Read data via OpenDAP
        [TC_err, OpenDAPncData] = MATROOSreadOpenDAP(startDate,networkData(1,:),networkFields);
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' 'HFR-MATROOS data successfully read via OpenDAP.']);
            % Retrieve the timestamps of the files already converted
            try
                MATROOSconverted_selectquery = ['SELECT datetime FROM total_HFRnetCDF_tb WHERE datetime>' '''' startDate ''' AND network_id = ''HFR-MATROOS'''];
                MATROOSconverted_curs = exec(conn,MATROOSconverted_selectquery);
                disp(['[' datestr(now) '] - - ' 'Query to total_HFRnetCDF_tb table for retrieving the timestamps from HFR-MATROOS network already converted successfully executed.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Fetch data
            try
                MATROOSconverted_curs = fetch(MATROOSconverted_curs);
                MATROOSconverted_data = MATROOSconverted_curs.Data;
                disp(['[' datestr(now) '] - - ' 'Data of the timestamps from HFR-MATROOS network already converted successfully fetched from total_HFRnetCDF_tb table.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Retrieve the number of already converted timestamps for HFR-MATROOS network
            try
                numMATROOSconverted = rows(MATROOSconverted_curs);
                disp(['[' datestr(now) '] - - ' 'Number of already converted timestamps for HFR-MATROOS network successfully retrieved from total_HFRnetCDF_tb table.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Close cursor to station_tb table
            try
                close(MATROOSconverted_curs);
                disp(['[' datestr(now) '] - - ' 'Cursor to total_HFRnetCDF_tb table successfully closed.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Retrieve the indices (related to the netCDF data array) of the timestamps to be converted for HFR-MATROOS network
            if(numMATROOSconverted~=0)
                [TC_err, MATROOStoBeConvertedIndices] = MATROOScompareTimeStamps(cell2mat(MATROOSconverted_data),OpenDAPncData.time);
            else
                MATROOStoBeConvertedIndices = 1:length(OpenDAPncData.time);
            end
            
            % Convert files
            for MtBC_idx=1:length(MATROOStoBeConvertedIndices)
                % v2.2
                [TC_err, networkData(1,:), outputFilename,outputFilesize] = MATROOS2netCDF_v22(OpenDAPncData,MATROOStoBeConvertedIndices(MtBC_idx),networkData(1,:),networkFields,stationData,stationFields);
                if(TC_err==0)
                    disp(['[' datestr(now) '] - - ' outputFilename ' total netCDF v2.2 file successfully created and stored.']);
                end
                
                % Insert converted total info in total_HFRnetCDF_tb table
                try
                    if((TC_err==0) && (exist('outputFilename','var') ~= 0))
                        % Define timestamp
                        ts = datevec(OpenDAPncData.time(MATROOStoBeConvertedIndices(MtBC_idx)));
                        DBtimestamp = sprintf('%.4d %.2d %.2d %.2d %.2d 00',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
                        
                        % Evaluate datetime from, Time Stamp
                        [t2d_err,DateTime] = timestamp2datetime(DBtimestamp);
                        
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'check_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {outputFilename,networkData{1,network_idIndex},DBtimestamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),outputFilesize,0};
                        
                        % Append the product data into the total_HFRnetCDF_tb table on the database.
                        tablename = 'total_HFRnetCDF_tb';
                        datainsert(conn,tablename,addColnames,addData);
                        disp(['[' datestr(now) '] - - ' outputFilename ' total file information successfully inserted into total_HFRnetCDF_tb table.']);
                    end
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    TC_err = 1;
                end
                
                clear outputFilename outputFilesize;
                
            end
            
        else
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> Something went wrong in reading HFR-MATROOS data via OpenDAP.']);
        end
    % Process HFR-US networks    
    elseif(contains(networkData{1,network_idIndex},'HFR-US'))
        % Read data via OpenDAP
        [TC_err, OpenDAPncData] = USreadOpenDAP(startDate,networkData(1,:),networkFields);
        if(TC_err==0)
            disp(['[' datestr(now) '] - - ' networkData{1,network_idIndex} ' data successfully read via OpenDAP.']);
            % Retrieve the timestamps of the files already converted
            try
                USconverted_selectquery = ['SELECT datetime FROM total_HFRnetCDF_tb WHERE datetime>' '''' startDate ''' AND network_id =' '''' networkData{1,network_idIndex} ''''];
                USconverted_curs = exec(conn,USconverted_selectquery);
                disp(['[' datestr(now) '] - - ' 'Query to total_HFRnetCDF_tb table for retrieving the timestamps from ' networkData{1,network_idIndex} ' network already converted successfully executed.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Fetch data
            try
                USconverted_curs = fetch(USconverted_curs);
                USconverted_data = USconverted_curs.Data;
                disp(['[' datestr(now) '] - - ' 'Data of the timestamps from ' networkData{1,network_idIndex} ' network already converted successfully fetched from total_HFRnetCDF_tb table.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Retrieve the number of already converted timestamps for HFR-MATROOS network
            try
                numUSconverted = rows(USconverted_curs);
                disp(['[' datestr(now) '] - - ' 'Number of already converted timestamps for ' networkData{1,network_idIndex} ' network successfully retrieved from total_HFRnetCDF_tb table.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Close cursor to station_tb table
            try
                close(USconverted_curs);
                disp(['[' datestr(now) '] - - ' 'Cursor to total_HFRnetCDF_tb table successfully closed.']);
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Retrieve the indices (related to the netCDF data array) of the timestamps to be converted for HFR-US network
            if(numUSconverted~=0)
                [TC_err, UStoBeConvertedIndices] = UScompareTimeStamps(cell2mat(USconverted_data),OpenDAPncData.time);
            else
                UStoBeConvertedIndices = 1:length(OpenDAPncData.time);
            end
            
            % Convert files
            for UtBC_idx=1:length(UStoBeConvertedIndices)
                % v2.2
                [TC_err, networkData(1,:), outputFilename,outputFilesize] = US2netCDF_v22(OpenDAPncData,UStoBeConvertedIndices(UtBC_idx),networkData(1,:),networkFields,stationData,stationFields);
                if(TC_err==0)
                    disp(['[' datestr(now) '] - - ' outputFilename ' total netCDF v2.2 file successfully created and stored.']);
                end
                
                % Insert converted total info in total_HFRnetCDF_tb table
                try
                    if((TC_err==0) && (exist('outputFilename','var') ~= 0))
                        % Define timestamp
                        ts = datevec(OpenDAPncData.time(UStoBeConvertedIndices(UtBC_idx)));
                        DBtimestamp = sprintf('%.4d %.2d %.2d %.2d %.2d 00',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
                        
                        % Evaluate datetime from, Time Stamp
                        [t2d_err,DateTime] = timestamp2datetime(DBtimestamp);
                        
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'check_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {outputFilename,networkData{1,network_idIndex},DBtimestamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),outputFilesize,0};
                        
                        % Append the product data into the total_HFRnetCDF_tb table on the database.
                        tablename = 'total_HFRnetCDF_tb';
                        datainsert(conn,tablename,addColnames,addData);
                        disp(['[' datestr(now) '] - - ' outputFilename ' total file information successfully inserted into total_HFRnetCDF_tb table.']);
                    end
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    TC_err = 1;
                end
                
                clear outputFilename outputFilesize;
                
            end
            
        else
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> Something went wrong in reading ' networkData{1,network_idIndex} ' data via OpenDAP.']);
        end
    % Process networks having standard input files
    else        
        % Retrieve the total files to be converted
        try
            toBeConvertedTotals_selectquery = ['SELECT * FROM total_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' networkData{1,network_idIndex} ''' AND NRT_processed_flag = 0'];
            toBeConvertedTotals_curs = exec(conn,toBeConvertedTotals_selectquery);
            disp(['[' datestr(now) '] - - ' 'Query to total_input_tb table for retrieving the total files from ' networkData{1,network_idIndex} ' network to be converted successfully executed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Fetch data
        try
            toBeConvertedTotals_curs = fetch(toBeConvertedTotals_curs);
            toBeConvertedTotals_data = toBeConvertedTotals_curs.Data;
            disp(['[' datestr(now) '] - - ' 'Data of the total files from ' networkData{1,network_idIndex} ' network to be converted successfully fetched from total_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve column names
        try
            toBeConvertedTotals_columnNames = columnnames(toBeConvertedTotals_curs,true);
            disp(['[' datestr(now) '] - - ' 'Column names from total_input_tb table successfully retrieved.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Retrieve the number of totals to be converted for the current network
        try
            numToBeConvertedTotals = rows(toBeConvertedTotals_curs);
            disp(['[' datestr(now) '] - - ' 'Number of total files from ' networkData{1,network_idIndex} ' network to be converted successfully retrieved from total_input_tb table.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Close cursor to total_input_tb table
        try
            close(toBeConvertedTotals_curs);
            disp(['[' datestr(now) '] - - ' 'Cursor to total_input_tb table successfully closed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        try
            % Find the index of the filename field
            filenameIndexC = strfind(toBeConvertedTotals_columnNames, 'filename');
            filenameIndex = find(not(cellfun('isempty', filenameIndexC)));
            
            % Find the index of the filepath field
            filepathIndexC = strfind(toBeConvertedTotals_columnNames, 'filepath');
            filepathIndex = find(not(cellfun('isempty', filepathIndexC)));
            
            % Find the index of the extension field
            extensionIndexC = strfind(toBeConvertedTotals_columnNames, 'extension');
            extensionIndex = find(not(cellfun('isempty', extensionIndexC)));
            
            % Find the index of the timestamp field
            timestampIndexC = strfind(toBeConvertedTotals_columnNames, 'timestamp');
            timestampIndex = find(not(cellfun('isempty', timestampIndexC)));
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            TC_err = 1;
        end
        
        % Scan the tuv files to be converted
        for toBeConverted_idx=1:numToBeConvertedTotals
            TC_err = 0;
            try
                if (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.tuv')) % Codar data
                    % v2.2
                    [TC_err, networkData(1,:), outputFilename,outputFilesize] = tuv2netCDF_v22([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},networkData(1,:),networkFields,stationData,stationFields);
                elseif (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.cur_asc')) % WERA data
                    % v2.2
                    [TC_err, networkData(1,:), outputFilename,outputFilesize] = curAsc2netCDF_v22([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},networkData(1,:),networkFields,stationData,stationFields);
                elseif (strcmp(toBeConvertedTotals_data{toBeConverted_idx,extensionIndex}, '.asc')) % WERA data
                    % v2.2
                    [TC_err, networkData(1,:), outputFilename, outputFilesize] = ascTot2netCDF_v22([toBeConvertedTotals_data{toBeConverted_idx,filepathIndex} filesep toBeConvertedTotals_data{toBeConverted_idx,filenameIndex}],toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},networkData(1,:),networkFields,stationData,stationFields);
                end
                if(TC_err==0)
                    disp(['[' datestr(now) '] - - ' outputFilename ' total netCDF v2.2 file successfully created and stored.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Update NRT_processed_flag in total_input_tb table
            try
                if(TC_err==0)
                    % Find the index of the NRT_processed_flag field
                    NRT_processed_flagIndex = find(not(cellfun('isempty', strfind(toBeConvertedTotals_columnNames, 'NRT_processed_flag'))));
                    % Define a cell array containing the column names to be updated
                    updateColnames = {'NRT_processed_flag'};
                    
                    % Define a cell array that contains the data for insertion
                    updateData = {1};
                    
                    % Update the total_input_tb table on the database
                    tablename = 'total_input_tb';
                    whereclause = ['WHERE filename = ' '''' toBeConvertedTotals_data{toBeConverted_idx,filenameIndex} ''''];
                    update(conn,tablename,updateColnames,updateData,whereclause);
                    disp(['[' datestr(now) '] - - ' 'total_input_tb table successfully updated with NRT processed flag for timestamp ' toBeConvertedTotals_data{toBeConverted_idx,timestampIndex} '.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            % Insert converted total info in total_HFRnetCDF_tb table
            try
                if((TC_err==0) && (exist('outputFilename','var') ~= 0))
                    % Evaluate datetime from, Time Stamp
                    [t2d_err,DateTime] = timestamp2datetime(toBeConvertedTotals_data{toBeConverted_idx,timestampIndex});
                    
                    % Define a cell array containing the column names to be added
                    addColnames = {'filename' 'network_id' 'timestamp' 'datetime' 'creation_date' 'filesize' 'input_filename' 'check_flag'};
                    
                    % Define a cell array that contains the data for insertion
                    addData = {outputFilename,networkData{1,network_idIndex},toBeConvertedTotals_data{toBeConverted_idx,timestampIndex},DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),outputFilesize,toBeConvertedTotals_data{toBeConverted_idx,filenameIndex},0};
                    
                    % Append the product data into the total_HFRnetCDF_tb table on the database.
                    tablename = 'total_HFRnetCDF_tb';
                    datainsert(conn,tablename,addColnames,addData);
                    disp(['[' datestr(now) '] - - ' outputFilename ' total file information successfully inserted into total_HFRnetCDF_tb table.']);
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                TC_err = 1;
            end
            
            clear outputFilename outputFilesize;
            
        end
        
    end
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TC_err = 1;
end

%%

if(TC_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_TotalConversion successfully executed.']);
end

% pause(20);

return