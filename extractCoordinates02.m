function [] = extractCoordinates02(datafile)
% To Run this code:
% 1) put all fcsv files in the same folder
% 2) Make that folder your working directory (navigate to it in the bar
%    above the script editor window)
% 3) Run the ode by clicking "Run" in the editor tab

names = dir('*.fcsv');
dataTable = {};

for i = 60:66
    filename = names(i).name;
    % edit this based on your file naming scheme
    
    % the end should be end-## where the ## is the number of characters
    % following your species name. For example:
    % Apurp001_20.fcsv <-- ## = 8
    species = filename(1:end-9);
    individual = filename(end-8);
    % the end should be end-n1:end-n2 where the n's are the number of ch
    % caracters surrounding your vert percent. For example:
    % Apurp001_20.fcsv <-- n1 = 7; n2 = 5
    scale = filename(end-6:end-5);
    
    %Read In File
    opts = delimitedTextImportOptions("NumVariables", 14);

    % Specify range and delimiter
    opts.DataLines = [4, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["columnsid", "x", "y", "z", "ow", "ox", "oy", "oz", "vis", "sel", "lock", "label", "desc", "associatedNodeID"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["label", "desc"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["label", "desc"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["columnsid", "associatedNodeID"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["columnsid", "associatedNodeID"], "ThousandsSeparator", ",");

    A = readtable(filename, opts);
    
    XYZ = table2array(A(:,2:4));
    XYZpts = [];
    
    labels = table2array(A(:,12));
    stringList = ["center", "TL", "ML", "MR", "BL", "BR"];
    for j = 1:6
        rowOfInterest = contains(labels, stringList(j));
        XYZpts = [XYZpts; XYZ(rowOfInterest,:)];
    end
    
%     XYpts = [XYZpts(:,3),XYZpts(:,2)];
%     XYptsTransform = [];
%     for j = 1:6
%         XYptsTransform = [XYptsTransform,XYpts(j,:)];
%     end
    
%     p = polyfit([XYZpts(1,2);XYZpts(4:5,2)] , [XYZpts(1,3);XYZpts(4:5,3)],1);
%     yP = polyval(p, XYZpts(:,2));
%     XYZptsOG = XYZpts; XYZptsOG(:,3) = XYZptsOG(:,3)-XYZptsOG(1,3);
%     XYZpts(:,3) = XYZpts(:,3)-yP;

    Dpost = norm(XYZpts(5,:)-XYZpts(6,:));
    CBL = norm(XYZpts(2,:)-XYZpts(5,:));
    d = norm(XYZpts(3,:)-XYZpts(4,:));
        if XYZpts(2,3) > XYZpts(5,3)
            XYZpts = [XYZpts; XYZpts(2,1), XYZpts(2,2)-Dpost, XYZpts(6,3)+CBL];
        else
            XYZpts = [XYZpts; XYZpts(2,1), XYZpts(2,2)-Dpost, XYZpts(6,3)-CBL];
        end
    Dant = norm(XYZpts(2,:)-XYZpts(7,:));
    I = norm(XYZpts(2,:)-XYZpts(1,:));
    J = norm(XYZpts(7,:)-XYZpts(1,:));
    K = norm(XYZpts(5,:)-XYZpts(1,:));
    L = norm(XYZpts(6,:)-XYZpts(1,:));
    
    alphaAnt = acos((I^2+J^2-Dant^2)/(2*I*J))*180/pi;
    alphaPost = acos((K^2+L^2-Dpost^2)/(2*K*L))*180/pi;
    
    AsquishAnt = 0.5*norm(cross(XYZpts(2,:)-XYZpts(1,:),XYZpts(7,:)-XYZpts(1,:)));
    AsquishPost = 0.5*norm(cross(XYZpts(5,:)-XYZpts(1,:),XYZpts(6,:)-XYZpts(1,:)));
    AboneL = 0.5*norm(cross(XYZpts(5,:)-XYZpts(3,:),XYZpts(2,:)-XYZpts(3,:)));
    AboneR = 0.5*norm(cross(XYZpts(6,:)-XYZpts(4,:),XYZpts(7,:)-XYZpts(4,:)));
    Asquish = AsquishAnt+AsquishPost; Abone = AboneL+AboneR;
    PercSquish = Abone/Asquish;
    
    h = figure;
    plot(XYZpts(:,2), XYZpts(:,3),'bO')
%     hold on
%     plot(XYZptsOG(:,2), XYZptsOG(:,3),'rO')
    title(filename);
    axis equal
    uiwait(h)
    
    dataTable = [dataTable; [species, individual, scale, ...
        num2cell(Dpost), num2cell(Dant), num2cell(alphaPost), num2cell(alphaAnt),...
        num2cell(d), num2cell(CBL), num2cell(AsquishAnt), num2cell(AsquishPost),...
        num2cell(AboneL), num2cell(AboneR), num2cell(PercSquish)]];
end

writecell(dataTable,datafile)
end