%% ruvRadialQCtests_v11.m
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
%         bear: radial velocity bearing variable from radial data
%         lond: longitude coordinates of the grid for radial velocities
%         latd: latitude coordinates of the grid for radial velocities
%         owtr: Vector Over Water quality flags
%         etmp: temporal quality variable from radial data
%         head: radial velocity heading variable from radial data
%         radVel: radial velocities from radial data
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
% Date: May 30, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, overWater, varThr, tempDer, velThr, medFilt, avgRadBear, radVelMedianFiltered, radCount] = ruvRadialQCtests_v11(bear, lond, latd, owtr, etmp, head, radVel, Radial_QC_params)

display(['[' datestr(now) '] - - ' 'ruvRadialQCtests_v11.m started.']);

RQC_err = 0;

%% Prepare QC flag variables

try
    overall = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    overWater = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    varThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    tempDer = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    velThr = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    medFilt = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(owtr,1), size(owtr,2)));
    avgRadBear = netcdf.getConstant('NC_FILL_BYTE');
    radCount = netcdf.getConstant('NC_FILL_BYTE');
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    RQC_err = 1;
end

%%

%% Prepare QC tests

try
    % Variance Threshold QC test
    varVec = etmp.^2;
    
    % Velocity Threshold QC test
    maxspd_R = Radial_QC_params.VelThr;
    
    % Average Radial Bearing QC test
    avgBear_HEAD = mean(head(~isnan(head)));
    avgBear = mean(bear(~isnan(bear)));
    
    % Median Filter QC test
    radVelMedianFiltered = radVel;
    
    % Radial Count QC test
    radVectors = sum(sum(~isnan(radVel)));
    
    % Temporal Derivative QC test
    tempDer_Thr = Radial_QC_params.TempDerThr.threshold;
    
    % Check if the files of the previous two hours exist
    if ((exist(Radial_QC_params.TempDerThr.hour2) == 2) && (exist(Radial_QC_params.TempDerThr.hour1) == 2))
        tD_go = true;
        [filepath1h,name1h,ext1h] = fileparts(Radial_QC_params.TempDerThr.hour1);
    else
        tD_go = false;
    end
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    RQC_err = 1;
end

%%

%% Populate QC variables
if (RQC_err == 0)
    try
        % Over Water quality flags
        overWater(owtr==0) = 1;
        overWater(owtr==128) = 4;
        
        % Velocity Threshold quality flags
        velThr((abs(radVel) <= maxspd_R)) = 1;
        velThr((abs(radVel) > maxspd_R)) = 4;
        
        % Variance Threshold quality flags
        varThr((varVec > Radial_QC_params.VarThr)) = 4;
        varThr((varVec <= Radial_QC_params.VarThr)) = 1;
        
        % Temporal Derivative quality flags
        if (tD_go)
            tempDer1h = tempDer;
            % Extract the radial velocity fields from the previous two hours files
            radVel1h = ncread(Radial_QC_params.TempDerThr.hour1,'RDVA');
            radVel2h = ncread(Radial_QC_params.TempDerThr.hour2,'RDVA');
            for rVr=1:size(radVel,1)
                for rVc=1:size(radVel,2)
                    if (~isnan(radVel1h(rVr, rVc)))
                        if ((isnan(radVel(rVr,rVc))) || isnan(radVel2h(rVr,rVc)))
                            tempDer1h(rVr,rVc) = 0;
                        elseif ((abs(radVel(rVr,rVc) - radVel1h(rVr,rVc)) < tempDer_Thr) && (abs(radVel2h(rVr,rVc) - radVel1h(rVr,rVc)) < tempDer_Thr))
                            tempDer1h(rVr,rVc) = 1;
                        else
                            tempDer1h(rVr,rVc) = 4;
                        end
                    end
                end
            end
            
            % Modify the VART_QC variable of the nc file of the previous hour
            ncwrite(Radial_QC_params.TempDerThr.hour1,'VART_QC',tempDer1h);
            disp(['[' datestr(now) '] - - ' [name1h,ext1h] ' previous time step nc file successfully updated with the Temporal Derivative QC variable.']);
        end
        % Set the QC flag for the current hour to 0 (no QC performed)
        tempDer(~isnan(radVel)) = 0;
        
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
        for i=1:size(radVel,1)
            for j = 1:size(radVel,2)
                if(~isnan(radVel(i,j)))
                    % Check vector distances
                    jj = j + 1;
                    radiusFlag = 0; % flag saying if the cell is inside (0) or outside (1) the search radius
                    % Scan the grid horizontally toward right
                    while((jj <= size(radVel,2)) && radiusFlag == 0)
                        [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(i,jj) lond(i,jj)]);
                        if(cellDist > Radial_QC_params.MedFilt(1))
                            radiusFlag = 1;
                        end
                        jj = jj + 1;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically downward
                        ii = i + 1;
                        while((ii <= size(radVel,1)) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(ii,j) lond(ii,j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii + 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (jj -1) - j;
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid horizontally toward left
                        jj = j - 1;
                        while((jj > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(i,jj) lond(i,jj)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            jj = jj - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (ii -1) - i;
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically upward
                        ii = i - 1;
                        while((ii > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(ii,j) lond(ii,j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = j - (jj + 1);
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % all the grid lies inside the search radius
                        % Use the whole grid as window for median filter
                        vel4bear = radVel;
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = i - (ii + 1);
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                    end
                    
                    % Check bearing distances within vel4filt indices matrix
                    refBear = bear(i,j);
                    bearDist = bear(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                    vel4filt = vel4bear(((abs(bearDist-refBear)) <= Radial_QC_params.MedFilt(2)));
                    
                    % Evaluate quality flag
                    medVal = median(vel4filt(:));
                    if(abs(radVel(i,j) - medVal) <= Radial_QC_params.MedFilt(3))
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
                    if((velThr(ii,jj) == 1) && (overWater(ii,jj) == 1) && (medFilt(ii,jj) == 1) && (avgRadBear == 1) && (radCount == 1))
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
    
    % Previous hour data file
    try
        if (tD_go)
            % Extract the QC variables from the previous hour file
            overall1h = ncread(Radial_QC_params.TempDerThr.hour1,'QCflag');
            overWater1h = ncread(Radial_QC_params.TempDerThr.hour1,'OWTR_QC');
            medFilt1h = ncread(Radial_QC_params.TempDerThr.hour1,'MDFL_QC');
            velThr1h = ncread(Radial_QC_params.TempDerThr.hour1,'CSPD_QC');
            avgRadBear1h = ncread(Radial_QC_params.TempDerThr.hour1,'AVRB_QC');
            radCount1h = ncread(Radial_QC_params.TempDerThr.hour1,'RDCT_QC');
            % Fill the overall QC variable
            for ii=1:size(overall1h,1)
                for jj = 1:size(overall1h,2)
                    if(~isnan(overall1h(ii,jj)))
                        if((tempDer1h(ii,jj) == 1) && (velThr1h(ii,jj) == 1) && (overWater1h(ii,jj) == 1) && (medFilt1h(ii,jj) == 1) && (avgRadBear1h == 1) && (radCount1h == 1))
                            overall1h(ii,jj) = 1;
                        else
                            overall1h(ii,jj) = 4;
                        end
                    else
                        overall1h(ii,jj) = netcdf.getConstant('NC_FILL_BYTE');
                    end
                end
            end
            % Modify the overall QC variable and the date_update attribute of the nc file of the previous hour
            ncwrite(Radial_QC_params.TempDerThr.hour1,'QCflag',int8(overall1h));
            ncwriteatt(Radial_QC_params.TempDerThr.hour1,'/','date_update',char([datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z']));
            disp(['[' datestr(now) '] - - ' [name1h,ext1h] ' previous time step nc file successfully updated with the overall QC variable.']);
        end
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        RQC_err = 1;
    end
end

%%

if(RQC_err==0)
    display(['[' datestr(now) '] - - ' 'ruvRadialQCtests_v11.m successfully executed.']);
end

return