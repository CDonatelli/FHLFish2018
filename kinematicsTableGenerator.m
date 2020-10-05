function [T, A] = kinematicsTableGenerator(csvName)
    list = dir('*.mat');

    SwimmingSpeed = []; StrideLength = [];
    Species = []; Trial = [];

    for i = 1:length(list)
        NameStr = list(i).name;
        Struct = load(NameStr);     
        fns = fieldnames(Struct);
        Struct = Struct.(fns{1});
        
        Species = [Species; NameStr(1:5)];
        Trial = [Trial; str2double(NameStr(12))+1];
        SwimmingSpeed = [SwimmingSpeed; Struct.swimmingSpeed];
        StrideLength = [StrideLength; Struct.bendingStrideLength];

    end
    
    T = table(Species, Trial, SwimmingSpeed, StrideLength);
    A = table2array(T);
    writetable(T,csvName);
   
end