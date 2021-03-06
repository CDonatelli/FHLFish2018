function [T, A] = kinematicsTableGenerator(csvName)
    list = dir('*.mat');

    SwimmingSpeed = []; StrideLength = [];
    Species = []; Trial = [];

    for i = 1:length(list)
        NameStr = list(i).name;
        Struct = load(NameStr);
        Struct = Struct.Struct;
        
        Species = [Species; NameStr(1:5)];
        Trial = [Trial; str2mat(NameStr(12))+1];
        SwimmingSpeed = [SwimmingSpeed; Struct.swimmingSpeed];
        StrideLength = [StrideLength; Struct.bendingStrideLength];

    end
    
    T = table(Species, Trial, SwimmingSpeed, StrideLength);
    A = table2array(T);
    writetable(T,csvName);
   
end