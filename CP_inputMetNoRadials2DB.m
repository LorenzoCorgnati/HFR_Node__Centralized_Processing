%% CP_inputMetNoRadials2DB.m
% This function lists the MetNo hourly radials read from OpenDAP urls
% of the aggregated netCDF files from MetNo TDS and inserts into the HFR
% database the information needed for the combination of radial files into
% totals and for the generation of the radial and total data files
% according to the European standard data model.

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
%         iRDB_err: error flag (0 = correct, 1 = error)

% Author: Lorenzo Corgnati
% Date: February 7, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [iRDB_err] = CP_inputMetNoRadials2DB(conn,startDate,networkData,networkFields,stationData,stationFields)

warning('off', 'all');

iRDB_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_inputMetNoRadials2DB started.']);

%% Scan the stations

try
    % Find the index of the network_id field
    network_idIndexC = strfind(networkFields, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
    
    % Find the index of the TDS_root_url field
    TDS_root_urlIndexC = strfind(networkFields, 'TDS_root_url');
    TDS_root_urlIndex = find(not(cellfun('isempty', TDS_root_urlIndexC)));
    
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

try
    for station_idx=1:size(stationData,1)
        % Build the list of radials for the time period to be processed
        try
            [iRDB_err,matRad] = MetNoTDSradListBuilder(networkData{1,TDS_root_urlIndex},stationData{station_idx,station_idIndex},startDate,stationData{station_idx,inputPathIndex});
            if(iRDB_err == 0)
                disp(['[' datestr(now) '] - - ' 'Radials files from ' stationData{station_idx,station_idIndex} ' station successfully listed.']);
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            iRDB_err = 1;
        end
        
        % Scan the radials of the time period to be processed
        if(iRDB_err == 0)
            for rad_idx=1:length(matRad.mat)
                iRDB_err = 0;
                % Check if the current radial was already remapped
                try
                    if (exist(matRad.mat{rad_idx}, 'file') ~= 2)
                        % Remap the radial nc file to RUV struct and save the mat file
                        try
                            [iRDB_err,ruvRad] = remapMetNoNCradial2ruv(matRad.nc{rad_idx},matRad.timeIndex{rad_idx});
                        catch err
                            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                            iRDB_err = 1;
                        end
                        
                        % Save radial mat file
                        if(iRDB_err == 0)
                            try
                                % Retrieve the filename
                                [pathstr,name,ext]=fileparts(matRad.mat{rad_idx});
                                noFullPathName=[name ext];
                                save(matRad.mat{rad_idx}, 'ruvRad');
                                disp(['[' datestr(now) '] - - ' noFullPathName ' file successfully saved.']);
                            catch err
                                disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                                iRDB_err = 1;
                            end
                            
                            % Insert information about the ruv file into the database (if not yet present)
                            try
                                % Check if the current ruv file is already present on the database
                                dbRadials_selectquery = ['SELECT * FROM radial_input_tb WHERE datetime>' '''' startDate ''' AND network_id = ' '''' networkData{1,network_idIndex} ''' AND filename = ' '''' noFullPathName ''' ORDER BY timestamp'];
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
                                    % Read the file header
                                    radHeader = ruvRad.OtherMetadata.Header;
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
                                    ruvFileInfo = dir(matRad.mat{rad_idx});
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
                                    addData = {noFullPathName,pathstr,networkData{1,network_idIndex},stationData{station_idx,station_idIndex},TimeStamp,DateTime,(datestr(now,'yyyy-mm-dd HH:MM:SS')),ruvFilesize,ext,0};
                                    
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
                catch err
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
                    iRDB_err = 1;
                end
                
                clear ruvRad
            end
        end
        
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    iRDB_err = 1;
end

%%

if(iRDB_err==0)
    disp(['[' datestr(now) '] - - ' 'CP_inputMetNoRadials2DB successfully executed.']);
end

return