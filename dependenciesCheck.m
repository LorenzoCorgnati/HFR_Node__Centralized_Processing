[fList,pList] = matlab.codetools.requiredFilesAndProducts('ruv2netCDF_v31.m');

for i=1:length(fList)
    disp(fList{i});
end