function dataStruct = dltNoseTailKine(testData)
%Note, make sure that you don't input any NaN

dataStruct = struct;

    prompt = {'Fish ID:', 'Trial ID:', 'Frame Rate:', 'Fish Length (mm):','Scale (mm/pixel):'};
    dlgtitle = 'Fish Data: Leave defaults if unknown';
    dims = [1 35];
    definput = {'n/a','n/a', '120','0','0'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
        dataStruct.fishID = answer{1};
        dataStruct.trialID = answer{2};
        dataStruct.FrameRate = str2num(answer{3});
        dataStruct.fishLength = str2num(answer{4});
        dataStruct.scale = str2num(answer{5}); %mm/pixel
        
    length = dataStruct.fishLength;
    [m,~] = size(testData);

    dataStruct.tailPts = [smooth(testData(:,7)), smooth(testData(:,8))];
    dataStruct.nosePts = [smooth(testData(:,1)), smooth(testData(:,2))];

    fr = dataStruct.FrameRate;
    total = m/fr; 
    dataStruct.t = linspace(0,total,m);

    p = polyfit(dataStruct.tailPts(:,1), dataStruct.tailPts(:,2),2);     % fit line for the tail wave
    yT = polyval(p, dataStruct.tailPts(:,1)); % y values for that line

    tailY = smooth(dataStruct.tailPts(:,2)-yT);
    
    %%%%% Mini-Peak Finder
        [p1,k1] = findpeaks(tailY,dataStruct.t,'MinPeakProminence',5);
        [p2,k2] = findpeaks(-tailY,dataStruct.t,'MinPeakProminence',5);
        p = [p1;abs(p2)]; k = [k1';abs(k2')];
        peaks = [k,p]; peaks = sortrows(peaks);
        k = peaks(:,1); p = peaks(:,2);
        tailPeaks = [k,p];

    [tf,idx] = ismember(k,dataStruct.t);
    tailPeakX = dataStruct.tailPts(idx,2);
    
    tailAmps = abs(p*dataStruct.scale);
    
%     nose = [dataStruct.nosePts(1,1), dataStruct.nosePts(1,2); ...
%             dataStruct.nosePts(end,1), dataStruct.nosePts(end,2);];

    distance = arclength(dataStruct.nosePts(:,1), dataStruct.nosePts(:,2));

%     distance = pdist(nose, 'euclidean');
    distance = distance.*dataStruct.scale/1000;
    %speed in m/s
    dataStruct.swimmingSpeed = (distance/dataStruct.t(end));
    dataStruct.bendingFrequency = size(p,1)/2/dataStruct.t(end);
    dataStruct.bendingPeriod = 1/dataStruct.bendingFrequency;
    dataStruct.bendingStrideLength = distance/1000/(size(p,1)/2);
    dataStruct.bendingAmp = median(tailAmps)/1000;
    dataStruct.bendingAmps = tailAmps;
end