%% CP_inputRUV2DB.m
% This function lists the input ruv files pushed by the HFR data providers
% and insert into the HFR database the information needed for the
% combination of radial files into totals and for the generation of the
% radial and total data files into the European standard data model.

% INPUT:
%         conn: connection to database
%         startDate: processing start date
%         networkID: network_id of the processed network
%         stationData: cell array containing information about the stations
%                      of the processed network (metadata)
%         stationFields: field names of the cell array containing
%                       information about the stations of the processed network.

% OUTPUT:
%         iRDB_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: February 7, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [iRDB_err] = CP_inputRUV2DB(conn,startDate,networkID,stationData,stationFields)

warning('off', 'all');

iRDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputRUV2DB started.']);

%% Scan the stations

try
    % Find the index of the input file path field
    inputPathIndexC = strfind(stationFields, 'radial_input_folder_path');
    inputPathIndex = find(not(cellfun('isempty', inputPathIndexC)));
    
    % Find the index of the station_id field
    station_idIndexC = strfind(stationFields, 'station_id');
    station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end
for station_idx=1:size(stationData,1)
    % List the input ruv files for the current station
    if(~isempty(stationData{station_idx,inputPathIndex}))
        % Trim heading and trailing whitespaces from folder path
        stationData{station_idx,inputPathIndex} = strtrim(stationData{station_idx,inputPathIndex});
        % List the input ruv files for the current station
        try
            ruvFiles = rdir([stationData{station_idx,inputPathIndex} filesep '**' filesep '*.ruv'],'datenum>floor(now-3)');
            disp(['[' datestr(now) '] - - ' 'Radials files from ' stationData{station_idx,station_idIndex} ' station successfully listed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Insert information about the ruv file into the database (if not yet present)
        for ruv_idx=1:length(ruvFiles)
            iRDB_err = 0;
            try
                % Retrieve the filename
                [pathstr,name,ext]=fileparts(ruvFiles(ruv_idx).name);
                noFullPathName=[name ext];
                % Check if the current ruv file is already present on the database
                dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' networkID ''' AND filename = ' '''' noFullPathName ''' ORDER BY timestamp'];
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
                    % Load the radial file as structure
                    radStruct = loadRDLFile(ruvFiles(ruv_idx).name,'false','warning');
                    % Read the file header
                    radHeader = radStruct.OtherMetadata.Header;
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
                    ruvFileInfo = dir(ruvFiles(ruv_idx).name);
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
                    addData = {noFullPathName,pathstr,networkID,stationData{station_idx,station_idIndex},TimeStamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),ruvFilesize,ext,0};
                    
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

end

%%

if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_inputRUV2DB successfully executed.']);
end

return