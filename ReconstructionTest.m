
imageList = dir('*.jpg') ;
[m,n] = size(imageList);
theta = m/180;

rowVals = [];
for i = 1:1:m
    currentImage = imread(imageList(i).name);
    blueImage = currentImage(:,:,3);
    [imageM, imageN] = size(blueImage);
    rowVals = [rowVals; blueImage(1500,:)];
end
imageExample = imread(imageList(1).name);
imshow(rowVals);
output_size = max(size(imageExample));

imagesc(rowVals)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)'); 
ylabel('Parallel Sensor Position - x\prime (pixels)');


%[projections, Xp] = radon(rowVals,theta);

recondImage = iradon(rowVals',theta, output_size);

imshow(recondImage)