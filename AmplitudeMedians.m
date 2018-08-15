%Change your working directory to the folder containing all of the struct
%files you whish to extract data from

%update for the prefixes of the individuals in your data folder
prefixes = ['Aflav'; 'Ainsi'; 'Lsagi';'Plaet';'Rjord';'Xmuco'];
medians = [];
for i = 1:size(prefixes,1)
    % finds all the files with the ith prefix you listed
    files = dir([prefixes(i,:),'*']);
    BodyAmps = [];
    BodyAmpsRaw = [];
    for j = 1:size(files)
        % loads ethe jth file
        Data = load(files(j).name);
        Data = Data.Struct;
        BodyAmps = [BodyAmps, Data.BodyAmps'./Data.fishLength];
        BodyAmpsRaw = [BodyAmpsRaw, Data.BodyAmps'];
    end
    medians = median(BodyAmps,2);
    eval([[prefixes(i,:),'AmpsRaw'],'=BodyAmpsRaw'])
    eval([[prefixes(i,:),'Amps'],'=medians'])
end