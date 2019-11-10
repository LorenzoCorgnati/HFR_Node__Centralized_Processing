%% curAscTotalQCtests_v11.m
% This function performs the QC tests on total velocity data. The tests are
% the ones defined for the European common data and metadata model. This
% function is suited for the conversion of native .cur_asc WERA files in to
% netCDF files compliant to the European common data and metadata model.
% In particular, the following tests are performed:
%       - Velocity threshold
%       - GDOP threshold
%       - Variance threshold
%       - Temporal Derivative
%       - Data Density threshold
%       - Balance of contributing radials from different sites

% INPUT:
%         mat_tot: structure containing total file in Codar format
%         Total_QC_params: structure containing parameters for total QC tests

% OUTPUT:
%         overall: overall quality flag (good data value is assigned
%                         if and only if all QC tests are passed
%         varThr: Variance threshold quality flags
%         GDOPThr: GDOP threshold quality flags
%         dataDens: Data Density quality flags
%         velThr: Velocity threshold quality flags

% Author: Lorenzo Corgnati
% Date: June 6, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, varThr, GDOPThr, dataDens, velThr] = curAscTotalQCtests_v11(mat_tot, Total_QC_params)

disp(['[' datestr(now) '] - - ' 'curAscTotalQCtests_v11.m started.']);

TQC_err = 0;

%% Prepare QC flag variables

try
    overall = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(mat_tot.U_grid)));
    varThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(mat_tot.U_grid)));
    GDOPThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(mat_tot.U_grid)));
    dataDens = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(mat_tot.U_grid)));
    velThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(mat_tot.U_grid)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TQC_err = 1;
end

%%

%% Prepare variables for QC tests

try
    % Variance Threshold QC test
    variance = ((mat_tot.U_grid.^2).*(mat_tot.U_acc).^2) + ((mat_tot.V_grid.^2).*(mat_tot.V_acc).^2);
    
    % Velocity Threshold QC test
    totVel = sqrt(((mat_tot.U_grid).^2) + ((mat_tot.V_grid).^2));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TQC_err = 1;
end

%%

%% Populate QC variables
try
    % Velocity Threshold quality flags
    velThr(totVel>Total_QC_params.VelThr) = 4;
    velThr(totVel<=Total_QC_params.VelThr) = 1;
    
    % GDOP Threshold quality flags
    if(sum(sum(isnan(mat_tot.GDOP)))<numel(mat_tot.GDOP))
        GDOPThr(mat_tot.GDOP>Total_QC_params.GDOPThr) = 4;
        GDOPThr(mat_tot.GDOP<=Total_QC_params.GDOPThr) = 1;
    else
        GDOPThr(~isnan(totVel)) = 0;
    end
    
    % Data Density Threshold quality flag
    dataDens(~isnan(totVel)) = 1;
    
    % Variance Threshold quality flags
    varThr(variance>Total_QC_params.VarThr) = 4;
    varThr(variance<=Total_QC_params.VarThr) = 1;
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TQC_err = 1;
end

%%

%% Populate the overall quality variable

% Current data file
try
    for ii=1:size(overall,1)
        for jj = 1:size(overall,2)
            if(~isnan(mat_tot.U_grid(ii,jj)))
                if((varThr(ii,jj) == 1) && (velThr(ii,jj) == 1) && (GDOPThr(ii,jj) == 1) && (dataDens(ii,jj) == 1))
                    overall(ii,jj) = 1;
                else
                    overall(ii,jj) = 4;
                end
            end
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TQC_err = 1;
end

%%

if(TQC_err==0)
    disp(['[' datestr(now) '] - - ' 'curAscTotalQCtests_v11.m successfully ecexuted.']);
end

return