%% twoPastHours.m
% This function builds the strings in the form YYYY_MM_DD_HHHH related to
% the time stamps of the two previous hours with respect to the input time stamp.
% The function also builds the time-related part of the folder path for the
% files of two previous hours, in the form /YYYY/YYYY_MM/YYYY_MM_DD/.

% INPUT:
%         present: time stamp of the current hour
%         timeStep: temporal resolution of the current dataset

% OUTPUT:
%         past2: structure containing the time stamp and the time-related part
%                of the folder path of the file related to 2 time steps before the current hour
%         past1: structure containing the time stamp and the time-related part
%                of the folder path of the file related to 1 time step before the current hour


% Author: Lorenzo Corgnati
% Date: June 15, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%


function [past2, past1] = twoPastHours(present,timeStep)

disp(['[' datestr(now) '] - - ' 'twoPastHours.m started.']);

tPH_err = 0;

warning('off', 'all');

try
    % Evaluate the timestamps of the two past hours
    past2_num = present - 2*((timeStep)/(24*60));
    past1_num = present - timeStep/(24*60);
    
    past2_vec = datevec(past2_num);
    past1_vec = datevec(past1_num);
    
    % Retrieve the time elements
    year1 = num2str(past1_vec(1));
    month1 = sprintf('%02d',past1_vec(2));
    day1 = sprintf('%02d',past1_vec(3));
    hour1 = sprintf('%02d',past1_vec(4));
    minutes1 = sprintf('%02d',past1_vec(5));
    
    year2 = num2str(past2_vec(1));
    month2 = sprintf('%02d',past2_vec(2));
    day2 = sprintf('%02d',past2_vec(3));
    hour2 = sprintf('%02d',past2_vec(4));
    minutes2 = sprintf('%02d',past2_vec(5));
    
    % Build the time stamps of the two past hours
    past1.TS = [year1 '_' month1 '_' day1 '_' hour1 minutes1];
    past2.TS = [year2 '_' month2 '_' day2 '_' hour2 minutes2];
    
    % Build the time-related part of the folder paths
    past1.fP = [year1 filesep year1 '_' month1 filesep year1 '_' month1 '_' day1];
    past2.fP = [year2 filesep year2 '_' month2 filesep year2 '_' month2 '_' day2];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tPH_err = 1;
end

if(tPH_err==0)
    disp(['[' datestr(now) '] - - ' 'twoPastHours.m successfully executed.']);
end

return