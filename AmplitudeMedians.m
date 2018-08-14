
prefixes = ['Aflav'; 'Ainsi'; 'Lsagi';'Plaet';'Rjord';'Xmuco'];
medians = [];
for i = 1:6
    files = dir([prefixes(i,:),'*']);
    BodyAmps = [];
    BodyAmpsRaw = [];
    for j = 1:size(files)
        Data = load(files(j).name);
        Data = Data.Struct;
        BodyAmps = [BodyAmps, Data.BodyAmps'./Data.fishLength];
        BodyAmpsRaw = [BodyAmpsRaw, Data.BodyAmps'];
    end
    medians = median(BodyAmps,2);
    eval([[prefixes(i,:),'AmpsRaw'],'=BodyAmpsRaw'])
    eval([[prefixes(i,:),'Amps'],'=medians'])
end