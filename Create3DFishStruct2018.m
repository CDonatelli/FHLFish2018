function [ name ] = Create3DFishStruct2018( nameString )

name = struct;
videos = input('How many videos?: ');
% currentFolder = pwd;
% directory = uigetdir();
% cd('E:\FHL2015\Pictures\Best Trial Pictures\BWpics');
disp('Select the Dorsal Image');
dImName = uigetfile('','Select the Dorsal Image');
name.dorsalIm = imread(dImName);
disp('Select the Lateral Image');
lImName = uigetfile('','Select the Lateral Image');
name.lateralIm = imread(lImName);
name.fishLength = input('What is the length of the fish?: ');
% name.twistPts = input('Enter the twisting points ([paste]): ');
% cd(currentFolder);

%%%%%%%%% FOR 3D Analysis
% name = imageInfo(name);
%%%%%%%%%%

save(nameString, 'name');

    for i = 1:videos
        Struct = name;
        Struct.vid = input(['Name of video for struct ',num2str(i),': ']);
        eval([[nameString,'_',num2str(i)], '= Struct'])
        save([nameString,'_',num2str(i)], [nameString,'_',num2str(i)]);
    end

end
