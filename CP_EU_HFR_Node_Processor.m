%% CP_EU_HFR_Node_Processor.m
% This wrapper launches the scripts for inserting into the HFR database
% the information about radial and totala files (both Codar and WERA)
% pushed by the data providers and for combining radials into totals and
% converting radials and totals to netCDF files according to the European
% standard data model.

% Author: Lorenzo Corgnati
% Date: October 20, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

% Setup netCDF toolbox
setup_nctoolbox;

% Setup JBDC driver for MySQL
javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

% Setup map colormap
set(0,'DefaultFigureColormap',feval('jet'));

EHNP_err = 0;

disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Processor started.']);

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

%% Processing
% Set the infinite loop for continuous operation
kk = 5;
while(kk>0)
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Processor loop started.']);
    
    %% Set datetime of the starting date of the processing period
    
    try
        startDate = startCombinationDate(now);
        startDate = '2012-01-29';
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        EHNP_err = 1;
    end
    
    %%
    
    % RADIALS COMBINATION & RADIALS AND TOTALS CONVERSION
%     CP_inputRUV2DB;
    %     CP_inputCradAscii2DB;
    CP_HFRCombiner;
    
    % TOTALS CONVERSION
    CP_inputTUV2DB;
    CP_inputCurAsc2DB;
    CP_TotalConversion;
    
    disp(['[' datestr(now) '] - - ' 'CP_EU_HFR_Node_Processor loop ended.']);
end

%%