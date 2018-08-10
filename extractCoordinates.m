
names = dir('*.fcsv');

TenPercPoints = {};

for i = 1:length(names)
    filename = names(i).name;
    species = filename(1:end-9);
    perc = filename(end-8:end-7);
    
    %Read In File
    %filename = 'Z:\Donatelli\Elongate Fish Vert Mechanics\Raw data - ...
    %            CT Files\0.1Measurements-Kylene\AnoplarchusInsignis001.20%.fcsv';
    delimiter = ',';
    startRow = 4;
    formatSpec = '%*s%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
                'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    XYZpts = [dataArray{1:end-1}];
    
    XYpts = [XYZpts(:,3),XYZpts(:,2)];
    XYptsTransform = [];
    for j = 1:6
        XYptsTransform = [XYptsTransform,XYpts(j,:)];
    end
    
    % Clear temporary variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    TenPercPoints = [TenPercPoints; [species, perc, num2cell(XYptsTransform)]];
end