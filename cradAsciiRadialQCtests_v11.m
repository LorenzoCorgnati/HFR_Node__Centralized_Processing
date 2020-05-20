%% cradAsciiRadialQCtests_v11.m
% This function performs the QC tests on radial velocity data. The tests are
% the ones defined in INCREASE project compliant to CMEMS needs.
% In particular, the following tests are performed:
%       - Velocity threshold
%       - Variance threshold
%       - Temporal derivative
%       - Median Filter
%       - Average radial bearing
%       - Radial count
%       - Vector Over Water

% INPUT:
%         radData: structure containing radial data
%         Radial_QC_params: structure containing parameters for radial QC tests

% OUTPUT:
%         overall: overall quality flag (good data value is assigned
%                         if and only if all QC tests are passed
%         varThr: Variance threshold quality flags
%         tempDer: Temporal derivative quality flags
%         velThr: Velocity threshold quality flags
%         overWater: Vector Over Water quality flags
%         medFilt: Median Filter quality flags
%         avgRadBear: Average Radial Bearing quality flag
%         radVelMedianFiltered: radial velocities after the application of the median filter
%         radCount: Radial Count quality flag


% Author: Lorenzo Corgnati
% Date: May 18, 2020

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, overWater, varThr, tempDer, velThr, medFilt, avgRadBear, radVelMedianFiltered, radCount] = cradAsciiRadialQCtests_v11(radData, Radial_QC_params)

display(['[' datestr(now) '] - - ' 'cradAsciiRadialQCtests_v11.m started.']);

RQC_err = 0;

%% Prepare QC flag variables

try
    overall = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    overWater = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    varThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    tempDer = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    velThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    medFilt = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(radData.rdva)));
    avgRadBear = netcdf.getConstant('NC_FILL_BYTE');
    radCount = netcdf.getConstant('NC_FILL_BYTE');
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    RQC_err = 1;
end

%%

%% Prepare QC tests

try
    % Velocity Threshold QC test
    maxspd_R = Radial_QC_params.VelThr;
    
    % Average Radial Bearing QC test
    avgBear = mean(radData.drva(~isnan(radData.drva)));
    
    % Median Filter QC test
    radVelMedianFiltered = radData.rdva;
    
    % Radial Count QC test
    radVectors = sum(sum(~isnan(radData.rdva)));

catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    RQC_err = 1;
end

%%

%% Populate QC variables
if (RQC_err == 0)
    try
        % Over Water quality flags
        overWater(~isnan(radData.rdva)) = 1;
        
        % Velocity Threshold quality flags
        velThr((abs(radData.rdva) <= maxspd_R)) = 1;
        velThr((abs(radData.rdva) > maxspd_R)) = 4;
        
        % Variance Threshold quality flags
        varThr((radData.hcss > Radial_QC_params.VarThr)) = 4;
        varThr((radData.hcss <= Radial_QC_params.VarThr)) = 1;
        
        % Average Radial Bearing quality flag
        if ((avgBear >= Radial_QC_params.AvgRadBear(1)) && (avgBear <= Radial_QC_params.AvgRadBear(2)))
            avgRadBear = 1;
        else
            avgRadBear = 4;
        end
        
        % Radial Count quality flag
        if (radVectors > Radial_QC_params.RadCnt)
            radCount = 1;
        else
            radCount = 4;
        end
        
        % Median Filter quality flags
        for i=1:size(radData.rdva,1)
            for j = 1:size(radData.rdva,2)
                if(~isnan(radData.rdva(i,j)))
                    % Check vector distances
                    jj = j + 1;
                    radiusFlag = 0; % flag saying if the cell is inside (0) or outside (1) the search radius
                    % Scan the grid horizontally toward right
                    while((jj <= size(radData.rdva,2)) && radiusFlag == 0)
                        [cellDist, d2km] = lldistkm([radData.Lat(i) radData.Lon(j)], [radData.Lat(i) radData.Lon(jj)]);
                        if(cellDist > Radial_QC_params.MedFilt(1))
                            radiusFlag = 1;
                        end
                        jj = jj + 1;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically downward
                        ii = i + 1;
                        while((ii <= size(radData.rdva,1)) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([radData.Lat(i) radData.Lon(j)], [radData.Lat(ii) radData.Lon(j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii + 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (jj -1) - j;
                        vel4bear = radData.rdva(max(1,i-windowSpan):min(i+windowSpan,size(radData.rdva,1)),max(1,j-windowSpan):min(j+windowSpan,size(radData.rdva,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid horizontally toward left
                        jj = j - 1;
                        while((jj > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([radData.Lat(i) radData.Lon(j)], [radData.Lat(i) radData.Lon(jj)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            jj = jj - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (ii -1) - i;
                        vel4bear = radData.rdva(max(1,i-windowSpan):min(i+windowSpan,size(radData.rdva,1)),max(1,j-windowSpan):min(j+windowSpan,size(radData.rdva,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically upward
                        ii = i - 1;
                        while((ii > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([radData.Lat(i) radData.Lon(j)], [radData.Lat(ii) radData.Lon(j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = j - (jj + 1);
                        vel4bear = radData.rdva(max(1,i-windowSpan):min(i+windowSpan,size(radData.rdva,1)),max(1,j-windowSpan):min(j+windowSpan,size(radData.rdva,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % all the grid lies inside the search radius
                        % Use the whole grid as window for median filter
                        vel4bear = radData.rdva;
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = i - (ii + 1);
                        vel4bear = radData.rdva(max(1,i-windowSpan):min(i+windowSpan,size(radData.rdva,1)),max(1,j-windowSpan):min(j+windowSpan,size(radData.rdva,2)));
                    end
                    
                    % Check bearing distances within vel4filt indices matrix
                    refBear = radData.drva(i,j);
                    bearDist = radData.drva(max(1,i-windowSpan):min(i+windowSpan,size(radData.rdva,1)),max(1,j-windowSpan):min(j+windowSpan,size(radData.rdva,2)));
                    vel4filt = vel4bear(((abs(bearDist-refBear)) <= Radial_QC_params.MedFilt(2)));
                    
                    % Evaluate quality flag
                    medVal = median(vel4filt(:));
                    if(abs(radData.rdva(i,j) - medVal) <= Radial_QC_params.MedFilt(3))
                        medFilt(i,j) = 1;
                    else
                        radVelMedianFiltered(i,j) = medVal;
                        %                         medFilt(i,j) = 8;
                        medFilt(i,j) = 4;
                    end
                end
            end
        end
        
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        RQC_err = 1;
    end
end

%%

%% Populate the overall quality variable

if(RQC_err==0)
    try
        for ii=1:size(overall,1)
            for jj = 1:size(overall,2)
                if(velThr(ii,jj) ~= netcdf.getConstant('NC_FILL_BYTE'))
                    if((velThr(ii,jj) == 1) && (varThr(ii,jj) == 1) && (overWater(ii,jj) == 1) && (medFilt(ii,jj) == 1) && (avgRadBear == 1) && (radCount == 1))
                        overall(ii,jj) = 1;
                    else
                        overall(ii,jj) = 4;
                    end
                end
            end
        end
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        RQC_err = 1;
    end
    
end

%%

if(RQC_err==0)
    display(['[' datestr(now) '] - - ' 'cradAsciiRadialQCtests_v11.m successfully executed.']);
end

return