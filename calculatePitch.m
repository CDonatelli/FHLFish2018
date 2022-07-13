nose = [smooth(pitchDat(:,1)),smooth(pitchDat(:,2))];
tail = [smooth(pitchDat(:,3)),smooth(pitchDat(:,4))];
 
noseMin = min(nose(:,2)); tailMin = min(tail(:,2));
actualMin = min([noseMin, tailMin]);

[m,~] = size(nose);
ground = [nose(:,1), repmat(actualMin-10, [m, 1])];

angle = [];
for i = 1:m
    n1 = (tail(i,:) - nose(i,:)) / norm(tail(i,:) - nose(i,:));  % Normalized vectors
    n2 = (ground(i,:) - nose(i,:)) / norm(ground(i,:) - nose(i,:));
    angle1 = acos(dot(n1, n2)); 
    angle = [angle;angle1];
end
angle = angle*(180/pi);
medianAng = nanmedian(angle);

%https://www.mathworks.com/matlabcentral/answers/331017-calculating-angle-between-three-points
