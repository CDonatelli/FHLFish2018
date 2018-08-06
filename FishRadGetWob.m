function struct = FishRadGetWob(structfile)

load(structfile);
clear 'structfile'
vars = who;
eval(['struct=',cell2mat(vars(1))])

%disp('Select the directory containing your BW.tif files');
%direct = uigetdir;
%cd(direct);
%DataFile = input('What is the name of the structure file? :');
%ImageFile = input('What is the prefix of the image files? :');
%Format = input('What is the format of the image? :');
%FileNameList = dir([ImageFile,'*.tif']);
%PicNames = cell(length(FileNameList),1); % create a cell array                        
% for t = 1:length(FileNameList);          % put nnames into cell array
%     PicNames{t}= FileNameList(t).name;
% end
X = struct.X;
Y = struct.Y;
nfr = size(X,2);

for i = 1:nfr
    Thickness(i).Frame = DataFile(i).Frame; % create a structure file to 
                                            % add thicknesses to
    % Read in midline data from MidlineCust
    XF = DataFile(i).MidLine(:,1); YF = DataFile(i).MidLine(:,2);
%     if i <10 % read in ith image file
%         FI = imread([ImageFile,'0',num2str(i),'.tif']);
%     else
%         FI = imread([ImageFile,num2str(i),'.tif']);
%     end
    FI = imread(char(PicNames(i)));
    DFi = double(FI);    % convert image file to double
    imshow(DFi); hold on % show the image
    plot(XF,YF,'b');     % plot the midline
    for j = 2:length(XF)
        % create a new array of X's in order to create a line perpendicular
        % to the midline at the jth point
        Xup = XF(j)+25; Xend = XF(j)-25; Xs = Xend:0.25:Xup;
        % calculate the slope of the midline from j-1 to j
        M =(YF(j-1)-YF(j))/(XF(j-1)-XF(j));
        B = YF(j) - XF(j)*(1/(-M)); % calculate B of the perpendicular line
        if abs(M) <= (1e-3)         % if slope is close to 0 make it 0
            M = 0;
        end
        if M == 0            % if the slope is zero, create a vertical line
            perpY = (YF(j)-30):0.25:(YF(j)+29);
            x(1:0.5:length(perpY/4)) = XF(j); Xs = x;
        else                         % if the slope is not zero, calculate
            perpY = (1/(-M))*Xs + B; % set of Y's for perpendicular line
        end
        plot(Xs,perpY,'r');          % plot the perpendicular line
% -->   plot(Xs,perpY,'r*');   % debug, make sure it's getting enough pts
        big = [];                    % create an empty array to 
        [xx,yy] = size(DFi);         % get size of the image                 
        for p = 1:length(Xs)                   
            if perpY(p) > xx || perpY(p) <= 0.5 % fill with values of the
                big = [big, p];   
            end                                  % points where y is out of
            if Xs(p) > yy || Xs(p) <=0.5         %  the range in the tif
                big = [big, p];
            end        
        end
        Xs(big) = []; perpY(big) = [];        % subtract those points
        Black = [];             
        for k = 1:length(Xs)                  % create an array of points
            % CHANGE DEPENDING ON NATURE OF TIF FILES
            % 0 if Fish are black, 1 if fish are white
            if DFi(round(perpY(k)), round(Xs(k))) == 0 % where the perp line
                add = [Xs(k),perpY(k)];        % passes through a black spot
                Black = [Black; add];          % presumably the fish
            end
        end
        % plot(Black, 'co');
        % add the thicknesses at each point along the midline to a
        % structure array
        Thickness(i).Thick(j) = sqrt((Black(end,1)-Black(1,1))^2 + ...
            (Black(end,2)-Black(1,2))^2);
    end
end
eval(sprintf('%s=%s',[ImageFile,'Thk'],'Thickness'));       % rename output
save([ImageFile,'Thk.mat'],[ImageFile, 'Thk']);         % save midline data

end

function FrameOut = ProcessImage(Frame,Level)
h = ones(5,5) / 25;     %blur the image to kill line artifacts
BlurredImage = imfilter(Frame,h);
%CroppedImage = double(BlurredImage);%(65:405,80:660); %this removes the borders from the images
                                             % assumes motionscope frame
%Level = graythresh(BlurredImage)*.8;         %set threshold a little darker than the auto computed one
FrameOut = ~im2bw(BlurredImage,Level);       %make image binary and invert it so fish is white
FrameOut = bwareaopen(FrameOut, 75);
FrameOut(1:3,:) = [];
FrameOut(end-2:end,:) = [];
FrameOut(:,1:3) = [];
FrameOut(:,end-2:end) = [];
