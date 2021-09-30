%% CP_inputCradAscii2DB.m
% This function lists the input crad_ascii (WERA radials) files pushed by
% the HFR data providers and insert into the HFR database the information
% needed for the combination of radial files into totals and for the
% generation of the radial and total data files into the European standard
% data model.

% INPUT:
%         conn: connection to database
%         startDate: processing start date
%         networkID: network_id of the processed network
%         stationData: cell array containing information about the stations
%                      of the processed network (metadata)
%         stationFields: field names of the cell array containing
%                       information about the stations of the processed network.

% OUTPUT:
%         iCradDB_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: May 20, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [iCradDB_err] = CP_inputCradAscii2DB(conn,startDate,networkID,stationData,stationFields)

warning('off', 'all');

iCradDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputCradAscii2DB started.']);

startDateNum = datenum(startDate);

%%

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
    % List the input crad_ascii files for the current station
    if(~isempty(stationData{station_idx,inputPathIndex}))
        % Trim heading and trailing whitespaces from folder path
        stationData{station_idx,inputPathIndex} = strtrim(stationData{station_idx,inputPathIndex});
        % List the input crad_ascii files for the current station
        try
            cradFiles = rdir([stationData{station_idx,inputPathIndex} filesep '**' filesep '*.crad_ascii'],'datenum>floor(now-3)');
            disp(['[' datestr(now) '] - - ' 'Radial files from ' stationData{station_idx,station_idIndex} ' station successfully listed.']);
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iCradDB_err = 1;
        end
        
        % Insert information about the crad_ascii file into the database (if not yet present)
        for crad_idx=1:length(cradFiles)
            try
                iCradDB_err = 0;
                % Retrieve the filename
                [pathstr,name,ext]=fileparts(cradFiles(crad_idx).name);
                noFullPathName=[name ext];
                % Check if the file is related to current data (and not to wave data)
                if(~contains(noFullPathName, 'WAV'))
                    % Check if the current crad_ascii file is already present on the database
                    dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' networkID ''' AND filename = ' '''' noFullPathName ''' ORDER BY timestamp'];
                    dbRadials_curs = exec(conn,dbRadials_selectquery);
                    disp(['[' datestr(now) '] - - ' 'Query to radial_input_tb table for checking if ' noFullPathName ' radial file is already present in the database successfully executed.']);
                    
                    % Fetch data
                    dbRadials_curs = fetch(dbRadials_curs);
                    disp(['[' datestr(now) '] - - ' 'Data about the presence of ' noFullPathName ' radial file in the database successfully fetched from radial_input_tb table.']);
                    
                    if(rows(dbRadials_curs) == 0)
                        % Retrieve information about the crad_ascii file
                        % Read the timestamp from the header
                        [date,time] = textread(cradFiles(crad_idx).name, '%*15c %9c %*0c %5c',1);
                        TimeStampVec = datevec([date ' ' time]);
                        TimeStamp = [num2str(TimeStampVec(1)) ' ' num2str(TimeStampVec(2),'%02d') ' ' num2str(TimeStampVec(3),'%02d') ' ' num2str(TimeStampVec(4),'%02d') ' ' num2str(TimeStampVec(5),'%02d') ' ' num2str(TimeStampVec(6),'%02d')];
                        
                        % Evaluate datetime from, Time Stamp
                        [t2d_err,DateTime] = timestamp2datetime(TimeStamp);
                        
                        % Retrieve information about the crad_ascii file
                        cradFileInfo = dir(cradFiles(crad_idx).name);
                        cradFilesize = cradFileInfo.bytes/1024;
                        
                        % Write crad_ascii info in radial_input_tb table
                        % Define a cell array containing the column names to be added
                        addColnames = {'filename' 'filepath' 'network_id' 'station_id' 'timestamp' 'datetime' 'reception_date' 'filesize' 'extension' 'NRT_processed_flag'};
                        
                        % Define a cell array that contains the data for insertion
                        addData = {noFullPathName,pathstr,networkID,stationData{station_idx,station_idIndex},TimeStamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),cradFilesize,ext,0};
                        
                        % Append the product data into the radial_input_tb table on the database.
                        tablename = 'radial_input_tb';
                        datainsert(conn,tablename,addColnames,addData);
                        disp(['[' datestr(now) '] - - ' noFullPathName ' radial file information successfully inserted into radial_input_tb table.']);
                    end
                end
            catch err
                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                iCradDB_err = 1;
            end
        end
    end
end

%%

if(iCradDB_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_inputCradAscii2DB successfully executed.']);
end

return
