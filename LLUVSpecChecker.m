%% LLUVSpecChecker.m
% This function retrieves the LLUVSpec version from the header of Codar ruv and tuv files.

% INPUT:
%         header: header of the Codar ruv or tuv file.

% OUTPUT:
%         LSC_err: error flag (0 = correct, 1 = error)
%         LLUVSpec: string containing the LLUVSpec version


% Author: Lorenzo Corgnati
% Date: July 9, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [LSC_err,LLUVSpec] = LLUVSpecChecker(header)

disp(['[' datestr(now) '] - - ' 'LLUVSpecChecker.m started.']);

LSC_err = 0;

warning('off', 'all');

try
    % Retrieve information from header
    for header_idx=1:length(header)
        splitLine = regexp(header{header_idx}, ' ', 'split');
        
        % Retrieve site codes and coordinates
        if(strcmp(splitLine{1}, '%LLUVSpec:'))
            LLUVSpecStr = strrep(header{header_idx}(length('%LLUVSpec:')+2:length(header{header_idx})), '"', '');
            LLUVSpecCellArray = strsplit(LLUVSpecStr);
            LLUVSpec = LLUVSpecCellArray{1};
            break
        end
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    LSC_err = 1;
end

if(LSC_err==0)
    disp(['[' datestr(now) '] - - ' 'LLUVSpecChecker.m successfully executed.']);
end

return

