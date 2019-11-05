%% MetNoTDSradListBuilder.m
% This function builds the list of MetNo hourly radials from OpenDAP urls
% of the aggregated netCDF files from MetNo TDS.

% INPUT:
%         TDSroot: MetNo TDS OpenDAP url root path.
%         siteCode: code of the radial site.
%         starDate: start combination date.
%         matFolderPath: path of the radial mat file folder.

% OUTPUT:
%         dMNuB_err: error flag (0 = correct, 1 = error)
%         radialList: structure containing the list of radial mat files to
%                       be remapped, the OpenDAP url of the corresponding netCDF file and
%                       the time index within the corresponding aggregated netCDF file.


% Author: Lorenzo Corgnati
% Date: November 4, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [dMNuB_err,radialList] = MetNoTDSradListBuilder(TDSroot,siteCode,startDate,matFolderPath)

disp(['[' datestr(now) '] - - ' 'MetNoTDSradListBuilder.m started.']);

dMNuB_err = 0;

warning('off', 'all');

% Evaluate the number of days of the time period to be processed
try
    dayNum =ceil(now - datenum(startDate));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    dMNuB_err = 1;
    return
end

% Create the daily TDS urls
for day_idx=1:dayNum
    try
        curDay = datestr(now-(dayNum-day_idx),'yyyy-mm-dd');
        yearFolder = curDay(1:4);
        monthFolder = curDay(6:7);
        dayFolder = curDay(9:10);
        TDSurls{day_idx} = [TDSroot yearFolder '/' monthFolder '/' dayFolder '/' siteCode '/' 'RDLm_' siteCode '_' yearFolder '_' monthFolder '_' dayFolder '.nc'];
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        dMNuB_err = 1;
        return
    end
end

% Create the radial list
try
    rad_idx = 0;
    for day_idx=1:dayNum
        [pathstr,name,ext]=fileparts(TDSurls{day_idx});
        NC.time = (double(ncread(TDSurls{day_idx},'time'))/86400) + datenum('01-01-1970');
        for time_idx=1:length(NC.time)
            rad_idx = rad_idx + 1;
            curHour = strrep(datestr(NC.time(time_idx),'HH-MM'),'-','');
            radialList.mat{rad_idx} = [matFolderPath name '_' curHour '.mat'];
            radialList.nc{rad_idx} = TDSurls{day_idx};
            radialList.timeIndex{rad_idx} = time_idx;
        end
    end
catch err
    if(contains(err.message, 'file not found'))
        radialList = {};
        disp(['[' datestr(now) '] - - ' TDSurls{day_idx} ' file not present']);
    else
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    end
    dMNuB_err = 1;
    return
end

if(dMNuB_err==0)
    disp(['[' datestr(now) '] - - ' 'MetNoTDSradListBuilder.m successfully executed.']);
end

return