format long g
files = dir('*.csv');
TrialNames = {}; Indiv = {}; Trial = {};
Speed = []; Amps = []; Angs = []; Freq = []; Stride = []; Position = [];
npts = 20;
    
for i = 1:3:length(files)
    disp(['Trial ', num2str(i), ' out of ', num2str(length(files)), ' started.'])
    name = files(i).name(1:end-5);
    
    opts.DataLines = [1, Inf]; opts.ExtraColumnsRule = "ignore";
    opts.Delimiter = ","; opts.EmptyLineRule = "read";
    
    time = readtable([name,'T.csv'], 'ReadVariableNames',false); 
    endTime = time(1,end); timeEnd = table2array(endTime); 
    timeEnd = cell2mat(timeEnd);timeEnd = str2num(timeEnd);
    
    X = readtable([name,'X.csv']); X = table2array(X);
    Y = readtable([name,'Y.csv']); Y = table2array(Y);
    nfr = size(X, 2);
    time = linspace(0,timeEnd,nfr)';
    
    % For xiphister since there are only 20 pts instead of 21
%     x = []; y = [];
%     for j = 1:nfr
%         % Generate equation if the midline
%         [pts, deriv, funct] = interparc(npts, smooth(X(:,j)), smooth(Y(:,j)), 'spline');
%         % add those points to an array
%         x = [x,smooth(pts(:,1))]; y = [y,smooth(pts(:,2))];
%     end
%     X = x; Y = y;
    
    % Calculate scale using the fish length as the scale bar
        FishPixelLengths = [];
        % Loop through a few frames of the vide and calculate the
        % length of the midline at each frame
        for o = 1:round(nfr/15):nfr
            FishPixelLengths = [FishPixelLengths,...
                arclength(smooth(X(:,o)),smooth(Y(:,o)))];
        end
        % Use the known length of the fish and the median of the 
        % midline measurements to calculate the scale of the video
        fishLengthPx = median(FishPixelLengths);

%     polypPts = round(1:4.95:100);
%%%% 2D Wave Kinematics
    disp(['Calculating Amp + Ang for trial ', num2str(i), ' out of ', num2str(length(files)), '.'])
    Amplitudes = [];
    Angles = [];
    Waves = [];
    for j = 1:npts
        pointY = Y(j,:); pointX = X(j,:); % For Xmuco data
%         pointY = Y(polypPts(j),:); pointX = X(polypPts(j),:);

        p = polyfit(pointX, pointY,3);  % fit line for the tail wave
        yT = polyval(p, pointX);        % y values for that line
        pointY = pointY - yT;            % subtract y values to get amplitude
        pointY = smooth(pointY,15,'rlowess'); pointX = smooth(pointX,15, 'rlowess');
    %%%%% Peak Finder
        % For walking MPP = 1; For swimming MPP = 0.075
        [AmpPks,AmpLoc] = findpeaks(abs(pointY),time,'MinPeakProminence',0.075);
        Amplitudes = [Amplitudes; median(AmpPks)/fishLengthPx];
        Waves = [Waves; length(AmpPks)/2];
        
        % For walking MPP = 0.5; For swimming MPP = 0.05
        tAngles = smooth((atan2(pointY, pointX)*180/pi)');
        [AngPks,AngLoc] = findpeaks(abs(tAngles),time,'MinPeakProminence',0.05);
        Angles = [Angles; median(AngPks)];
        
%         subplot(1,2,1);
%         plot(time, abs(pointY))
%         hold on
%         plot(AmpLoc, AmpPks, 'ro')
%         subplot(1,2,2);
%         plot(time, abs(tAngles))
%         hold on
%         plot(AngLoc, AngPks, 'ro')
%         pause
%         close all
        
    end
    nose = [X(1,1),Y(1,1);X(1,end),Y(1,end)];
    distance = pdist(nose, 'euclidean');
    distance = distance./fishLengthPx;                  
    swimmingSpeed = (distance/time(end));
    strideLength = distance/Waves(end);
    bendingFrequency = Waves/time(end);
    
%     %PolypNaming
%     TrialNames(end+1:end+npts,1) = {name(1:end-4)};
%     Indiv(end+1:end+npts,1) = {name(1:8)}; 
%     Trial(end+1:end+npts,1) = {name(10:13)};

    %XmucoNaming
%     TrialNames(end+1:end+npts,1) = {name(1:end-4)};
%     Indiv(end+1:end+npts,1) = {name(1:6)}; 
%     Trial(end+1:end+npts,1) = {name(8)};
    TrialNames(end+1,1) = {name(1:end-4)};
    Indiv(end+1,1) = {name(1:6)}; 
    Trial(end+1,1) = {name(8)};
    
%     speedVec = zeros(npts,1); speedVec(:,1) = swimmingSpeed;
    Speed = [Speed;swimmingSpeed]; 
%     slVec = zeros(npts,1); slVec(:,1) = strideLength;
   	Stride = [Stride;strideLength];
%     Position = [Position; round(linspace(0,npts,npts)./npts.*100)'];
    Amps = [Amps; Amplitudes(end)]; Angs = [Angs; Angles(end)];
    Freq = [Freq; bendingFrequency(end)];
        
end

T = table(TrialNames, Indiv, Trial, Position, Speed, Stride, Amps, Angs, Freq);
writetable(T,csvName,'WriteRowNames',true);