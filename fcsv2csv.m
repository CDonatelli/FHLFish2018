% To Run this code:
% 1) put all fcsv files in the same folder
% 2) Make that folder your working directory (navigate to it in the bar
%    above the script editor window)
% 3) Run the ode by clicking "Run" in the editor tab

names = dir('*.fcsv');
dataTable = {};

for i = 1:length(names)
    filename = names(i).name;
    
    %Read In File
    opts = delimitedTextImportOptions("NumVariables", 14);

    % Specify range and delimiter
    opts.DataLines = [4, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["ColumnsId", "x", "y", "z", "ow", "ox", "oy", "oz", "vis", "sel", "lock", "label", "desc", "associatedNodeID"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "desc", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "desc", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["ColumnsId", "label", "associatedNodeID"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["ColumnsId", "label", "associatedNodeID"], "ThousandsSeparator", ",");

    % Import the data
    slicerDat = readtable(filename, opts);
    clear opts

    writetable(slicerDat,[filename(1:end-5),'_MAT.csv'])
end