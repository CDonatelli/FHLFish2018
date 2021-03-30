
function [T, A] = kinematicsTableGenerator(csvName)
    list = dir('*bw.mat');

    SwimmingSpeed = []; StrideLength = [];
    TailFrequency = []; TailAmplitude = [];
    Fish = []; Trial = [];

    for i = 1:length(list)
        NameStr = list(i).name;
        Struct = load(NameStr);     
        fns = fieldnames(Struct);
        Struct = Struct.(fns{1});
        
        Fish = [Fish; Struct.fishID];
        Trial = [Trial; Struct.trialID];
        SwimmingSpeed = [SwimmingSpeed; Struct.swimmingSpeed];
        StrideLength = [StrideLength; Struct.strideLength];
        TailFrequency = [ TailFrequency; Struct.bendingFrequency];
        TailAmplitude = [TailAmplitude; Struct.Amplitudes];

    end
    
    T = table(Fish, Trial, SwimmingSpeed, StrideLength, TailFrequency, TailAmplitude);
    A = table2array(T);
    writetable(T,csvName);
   
end