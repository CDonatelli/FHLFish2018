% To Run this code:
% 1) put all fcsv files in the same folder
% 2) Make that folder your working directory (navigate to it in the bar
%    above the script editor window)
% 3) Run the ode by clicking "Run" in the editor tab

names = dir('*.fcsv');
dataTable = {};

for i = 1:length(names)
    filename = names(i).name;
    % edit this based on your file naming scheme
    
    % the end should be end-## where the ## is the number of characters
    % following your species name. For example:
    % Apurp001_20.fcsv <-- ## = 8
%     species = filename(1:5);
    individual = filename(1:end-5);
    % the end should be end-n1:end-n2 where the n's are the number of ch
    % caracters surrounding your vert percent. For example:
    % Apurp001_20.fcsv <-- n1 = 7; n2 = 5
%     scale = filename(16:20);
    
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
    
%     XYpts = [XYZpts(:,3),XYZpts(:,2)];
%     XYptsTransform = [];
%     for j = 1:6
%         XYptsTransform = [XYptsTransform,XYpts(j,:)];
%     end
    snout = XYZpts(1,:);
    CP = XYZpts(2,:);
    skull = XYZpts(3,:);
    tailStart = XYZpts(4,:);
    NCleft = XYZpts(5,:);
    NCright = XYZpts(6,:);
    
    standardLength = norm(snout-CP);
    abdLength = norm(skull-tailStart);
    headLength = norm(snout-skull);
    tailLength = norm(tailStart-CP);
    headWidth = norm(NCleft-NCright);
    
    BTratio = abdLength/standardLength;
    Bratio = abdLength/standardLength;
    Tratio = tailLength/standardLength;
    Hratio = headLength/standardLength;
    headAR = headLength/headWidth;
    
%     height = norm(XYZpts(1,:)-XYZpts(2,:));
%     width = norm(XYZpts(3,:)-XYZpts(4,:));
%     overlap = norm(XYZpts(3,:)-XYZpts(6,:));
%     underlap = norm(XYZpts(4,:)-XYZpts(5,:));
    
    % Clear temporary variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    dataTable = [dataTable; [individual, ...
        num2cell(standardLength), num2cell(abdLength), num2cell(headLength), ...
        num2cell(tailLength), num2cell(headWidth), num2cell(BTratio), ...
        num2cell(Bratio), num2cell(Tratio), num2cell(Hratio), num2cell(headAR)]];
end

writecell(dataTable,"dataFile.csv")