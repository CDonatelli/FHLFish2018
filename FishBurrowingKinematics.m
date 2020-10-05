function struct = FishBurrowingKinematics(structfile)
    load(structfile);
    clear 'structfile'
    vars = who;
    eval(['struct=',cell2mat(vars(1))]);
    %This imports the data

    prompt = {'Enter fish length in mm (leave as 0 if unknown)'};
    dlgtitle = 'Fish Length';
    dims = [1 35];
    definput = {'0'};
    lengthMM = inputdlg(prompt,dlgtitle,dims,definput);
    lengthMM = str2double(lengthMM);
    
    mids = struct.midLines;
    nfr = size(mids,2);
    % Calculate the scale of the video using the fish length
    % as the scale bar
        FishPixelLengths = [];
        % Loop through a few frames of the vide and calculate the
        % length of the midline at each frame
        for i = 1:round(nfr/15):nfr
            FishPixelLengths = [FishPixelLengths,...
                arclength(mids(i).MidLine(:,1),mids(i).MidLine(:,2))];
        end
        % Use the known length of the fish and the median of the 
        % midline measurements to calculate the scale of the video
        fishPixels = median(FishPixelLengths);
        
        struct.fishLength = lengthMM;
        VidScale = lengthMM/fishPixels;
        struct.VidScale = VidScale;
        
    % points to generate data for
    npts = 21;
    % initiating new variables
    x = []; y = []; tailPtCordsY = []; tailPtCordsX = [];
    for i = 1:nfr
        xPts = mids(i).MidLine(:,1); yPts = mids(i).MidLine(:,2);
        randPts = rand(1,length(xPts))/1000; xPts = xPts+randPts';
        % Generate equation if the midline
        [pts, deriv, funct] = interparc(npts, xPts, yPts, 'spline');
        % add those points to an array
        x = [x,pts(:,1)]; y = [y,pts(:,2)];
        
    end
    struct.X = x; struct.Y = y;
    % figure out time for each frame and make a vector of times
    
    prompt = {'Enter video frame rate.)'};
    dlgtitle = 'Frame Rate';
    dims = [1 35];
    definput = {'0'};
    fr = inputdlg(prompt,dlgtitle,dims,definput);
    fr = str2double(fr);
    
    total = nfr/fr; 
    struct.t = linspace(0,total,nfr)';
    s = linspace(1,lengthMM,npts);
    struct.s = s';
    
    Amplitudes = []; Waves = []; Angles = [];

    for j = 1:npts
        pointY = y(j,:); pointX = x(j,:);
        p = polyfit(pointX, pointY,3);  % fit line for the tail wave
        yT = polyval(p, pointX);        % y values for that line
        pointY = pointY - yT;           % subtract y values to get amplitude
        pointY = smooth(pointY,15,'rlowess'); pointX = smooth(pointX,15, 'rlowess');
        
        %%%%% Peak Finder
        [AmpPks,AmpLoc] = findpeaks(abs(pointY),struct.t,'MinPeakProminence',0.075);
        Amplitudes = [Amplitudes; median(AmpPks)/fishPixels];
        Waves = [Waves; length(AmpPks)/2];
        
        % For walking MPP = 0.5; For swimming MPP = 0.05
        tAngles = smooth((atan2(pointY, pointX)*180/pi)');
        [AngPks,AngLoc] = findpeaks(abs(tAngles),struct.t,'MinPeakProminence',0.05);
        Angles = [Angles; median(AngPks)];
        
%                 subplot(1,2,1);
%                 plot(struct.t, abs(pointY))
%                 hold on
%                 plot(AmpLoc, AmpPks, 'ro')
%                 subplot(1,2,2);
%                 plot(struct.t, abs(tAngles))
%                 hold on
%                 plot(AngLoc, AngPks, 'ro')
%                 pause
%                 close all
        
    end
    nose = [x(1,1),y(1,1);x(1,end),y(1,end)];
    distance = pdist(nose, 'euclidean');
    struct.distance = distance./fishPixels;                  
    struct.swimmingSpeed = (struct.distance/struct.t(end));
    struct.strideLength = distance/Waves(end);
    struct.bendingFrequency = Waves/struct.t(end);
    struct.Amplitudes = Amplitudes; 
    struct.Waves = Waves; 
    struct.Angles = Angles;
    
    eval([cell2mat(vars(1)), '= struct'])
    save(cell2mat(vars(1)), cell2mat(vars(1)));
end
