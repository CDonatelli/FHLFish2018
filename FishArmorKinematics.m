
files = dir('*.mat');   

for i = 1:length(files)
    
    Struct = load(files(i).name);
    Struct = Struct.midlineStruct;
    
    if isfield(Struct, 'fishLength') == 1
      lengthMM = Struct.fishLength;
    else
        prompt = {'Enter fish length in mm (leave as 0 if unknown)'};
        dlgtitle = 'Fish Length';
        dims = [1 35];
        definput = {'0'};
        lengthMM = inputdlg(prompt,dlgtitle,dims,definput);
        lengthMM = str2double(lengthMM);
    end
    
    
    mids = Struct.midLines;
    nfr = size(mids,2);
        
    % points to generate data for
%     npts = 21;
%     % initiating new variables
%     x = []; y = []; tailPtCordsY = []; tailPtCordsX = [];
%     Curvature = [];
%     for i = 1:nfr
%         xPts = sgolayfilt(mids(i).MidLine(:,1),1,15); 
%         yPts = sgolayfilt(mids(i).MidLine(:,2),1,15);
%         randPts = rand(1,length(xPts))/1000; xPts = xPts+randPts';
%         % Generate equation if the midline
%         [pts, deriv, funct] = interparc(npts, xPts, yPts, 'spline');
%         % add those points to an array
%         x = [x,pts(:,1)]; y = [y,pts(:,2)];
%         % L: Arc length R: Curvature radius K: Curvature vector
% %         [L, R, K] = curvature(x,y);
% %         Curvature = [Curvature, abs(R)./fishPixels];
%     end
    X = Struct.X; Y = Struct.Y;
    % figure out time for each frame and make a vector of times
    
        % Calculate the scale of the video using the fish length
    % as the scale bar
        FishPixelLengths = [];
        % Loop through a few frames of the vide and calculate the
        % length of the midline at each frame
        for j = 1:round(nfr/15):nfr
            FishPixelLengths = [FishPixelLengths,...
                arclength(Y(:,j),X(:,j))];
        end
        % Use the known length of the fish and the median of the 
        % midline measurements to calculate the scale of the video
        fishPixels = median(FishPixelLengths);
        
        Struct.fishLength = lengthMM;
        Struct.fishPixels = fishPixels;
        VidScaleCalc = lengthMM/fishPixels;
        VidScale = 5/11.67; %mm/pixel
        Struct.VidScale = VidScale;
        Struct.VidScaleCalc = VidScaleCalc;
    
%     prompt = {'Enter video frame rate.)'};
%     dlgtitle = 'Frame Rate';
%     dims = [1 35];
%     definput = {'0'};
%     fr = inputdlg(prompt,dlgtitle,dims,definput);
%     fr = str2double(fr);
    fr = Struct.originalFrameRate/Struct.skipRate;
    
    total = nfr/fr; 
    Struct.t = linspace(0,total,nfr)';
    
    Amplitudes = []; Waves = [];
    
    pkVal = 0;
    pkProm = 2.5;
    
%     while pkVal == 0
        pointY = X(21,:); pointX = Y(21,:);
        p = polyfit(pointX, pointY,3);  % fit line for the tail wave
        yT = polyval(p, pointX);        % y values for that line
        pointY = pointY - yT;           % subtract y values to get amplitude
        pointY = smooth(pointY,15,'rlowess'); pointX = smooth(pointX,15, 'rlowess');
        [AmpPks,AmpLoc] = findpeaks(pointY,Struct.t,'MinPeakProminence',pkProm);
        plot(Struct.t, pointY) 
        hold on
        plot(AmpLoc, AmpPks, 'bo')
        hold off
        
%         answer = questdlg('Are you getting the correct peaks?', ...
%             ['Current value is' num2str(pkProm)], 'Yes','No','No');
%         % Handle response
%         switch answer
%             case 'Yes'
%                 pkVal = 1;
%             case 'No'
%                 pkVal = 0;
%                 prompt = {'Enter a new prominence value.)'};
%                 dlgtitle = ['Current value is' num2str(pkProm)];
%                 dims = [1 35];  definput = {num2str(pkProm)};
%                 pkProm = inputdlg(prompt,dlgtitle,dims,definput);
%                 pkProm = str2double(pkProm);
%         end
%     end

%     for j = 1:npts
%         pointY = y(j,:); pointX = x(j,:);
%         p = polyfit(pointX, pointY,3);  % fit line for the tail wave
%         yT = polyval(p, pointX);        % y values for that line
%         pointY = pointY - yT;           % subtract y values to get amplitude
%         pointY = smooth(pointY,15,'rlowess'); pointX = smooth(pointX,15, 'rlowess');
%         
%         %%%%% Peak Finder
%         [AmpPks,AmpLoc] = findpeaks(abs(pointY),Struct.t,'MinPeakProminence',pkProm);
        Amplitudes = [Amplitudes; prctile(AmpPks,95)*VidScale/Struct.fishLength];
        Waves = [Waves; length(AmpPks)/2];
% %         MaxCurvature = [MaxCurvature; max(Curvature(:,j))];
%         
%         % For walking MPP = 0.5; For swimming MPP = 0.05
% %         tAngles = smooth((atan2(pointY, pointX)*180/pi)');
% %         [AngPks,AngLoc] = findpeaks(abs(tAngles),struct.t,'MinPeakProminence',0.05);
% %         Angles = [Angles; median(AngPks)];
%         
% %                 plot(struct.t, abs(pointY))
% %                 hold on
% %                 plot(AmpLoc, AmpPks, 'ro')
% %                 pause
% %                 close all
%         
%     end
    
    t0 = find(Struct.t == min(AmpLoc)); t1 = find(Struct.t == max(AmpLoc));
    noseAtPeaks = [Y(1,t0),X(1,t0);Y(1,t1),X(1,t1)];
    peak2peakDist = pdist(noseAtPeaks, 'euclidean');
    peak2peakDist = peak2peakDist*VidScale/Struct.fishLength;
    
    nose = [Y(1,1),X(1,1);Y(1,end),X(1,end)];
    distance = pdist(nose, 'euclidean');
    Struct.distanceTotal = distance*VidScale/Struct.fishLength;
    Struct.distanceP2P = peak2peakDist;
    Struct.swimmingSpeed = peak2peakDist/(Struct.t(t1)-Struct.t(t0));        
    Struct.strideLength = peak2peakDist/Waves(end);
    Struct.bendingFrequency = Waves./(Struct.t(t1)-Struct.t(t0));
    Struct.Amplitudes = Amplitudes; 
    Struct.Waves = Waves; 
%     struct.Angles = Angles;
%     struct.Curvature = Curvature;
%     struct.maxCurvature = MaxCurvature;
    
    save([files(i).name(1:end-17),'.mat'], 'Struct');
end
