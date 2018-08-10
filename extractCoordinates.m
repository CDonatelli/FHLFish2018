% To Run this code:
% 1) put all fcsv files in the same folder
% 2) Make that folder your working directory (navigate to it in the bar
%    above the script editor window)
% 3) Run the ode by clicking "Run" in the editor tab

names = dir('*.fcsv');
TenPercPoints = {};

for i = 1:length(names)
    filename = names(i).name;
    % edit this based on your file naming scheme
    
    % the end should be end-## where the ## is the number of characters
    % following your species name. For example:
    % Apurp001_20.fcsv <-- ## = 8
    species = filename(1:end-8);
    
    % the end should be end-n1:end-n2 where the n's are the number of ch
    % caracters surrounding your vert percent. For example:
    % Apurp001_20.fcsv <-- n1 = 7; n2 = 5
    perc = filename(end-6:end-5);
    
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