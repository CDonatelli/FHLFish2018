function [T, A] = wobbleTableGenerator(csvName)
    list = dir('*.mat');

    SwimmingSpeed = [];
    F = [];
    T = [];
    SL = [];
    lambda = [];
    WS = [];
    slip = [];
    amplitude = [];
    wobF = [];
    wobT = [];
    wobSL = [];
    wobLambda = [];
    wobAmp = [];
    averageVerts = [];
    guageLength = [];
    RowNames = [];
    for i = 2:length(list)
        NameStr = list(i).name;
        Struct = load(NameStr);
        Struct = Struct.Struct;
        
        individual = NameStr(1:end-6);
        trial = NameStr(end-5);
        SwimmingSpeed = [SwimmingSpeed; Struct.swimmingSpeed];
        F = [F; Struct.bendingFrequency];
        T = [T; Struct.bendingPeriod];
        SL = [SL; Struct.bendingStrideLength];
        lambda = [lambda; Struct.wavelength];
        WS = [WS; Struct.bendingWS];
        slip = [slip; Struct.bendingSlip];
        amplitude = [amplitude; Struct.bendingAmp];
        wobF = [wobF; Struct.wobbleFrequency];
        wobT = [wobT; Struct.wobblePeriod];
        wobSL = [wobSL; Struct.wobbleStrideLength];
        wobLambda = [wobLambda; Struct.wobbleWavelength];
        wobAmp = [wobAmp; Struct.wobbleAmp];
        guageLength = [guageLength; Struct.guageLength];
        averageVerts = [averageVerts; Struct.averageVerts];
        RowNames = [RowNames; NameStr(1:8)];
    end
    RowNames = cellstr(RowNames);
    T = table(  SwimmingSpeed, ...
                F, ...
                T, ...
                SL, ...
                lambda, ...
                WS,...
                slip,...
                amplitude, ...
                wobF, ...
                wobT, ...
                wobSL, ...
                wobLambda, ...
                wobAmp, ...
                guageLength, ...
                averageVerts, ...
                'RowNames',RowNames);
   A = table2array(T);
   writetable(T,csvName);
   
   
   
end