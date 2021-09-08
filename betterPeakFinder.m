function [k,p] = betterPeakFinder(X,Y)

[p1,k1] = findpeaks(Y,X,'MinPeakProminence',2);
[p2,k2] = findpeaks(-Y,X,'MinPeakProminence',2);

p = [p1'; -p2']; k = [k1'; k2'];
peaks = [k,p]; 
peaks = sortrows(peaks);
k = peaks(:,1); p = peaks(:,2);
        
% check if peakfinder did it's damn job
    activeFigure = figure;
    plot(X,Y);
    hold on
    plot(k,p,'r*');
    
    prompt = {'How Many False Peaks?', 'How Many Missing Peaks?'};
    BoxName = 'FindPeak Error correction';
    default = {'0','0'};
    answer = inputdlg(prompt, BoxName,1,default);
    answer = str2double(answer);
    
    % if the user needs to eliminate points
    if answer(1) ~= 0
        uiwait(msgbox({ 'Click and drag to draw a', ...
                        'square under incorrect peaks.'}));
        % to eliminate peaks
        for i = 1:answer(1)
            rect = getrect();
            elim = find(k>rect(1) & k<(rect(1)+rect(3)));
            k(elim) = []; p(elim) = [];
            plot(k,p,'o');
        end
    else
    end
    
    % if the user needs to specify points
    if answer(2) ~= 0
        uiwait(msgbox({ 'Click the point on the graph', ...
                        'to specify missing peaks.', ...
                        'Then press ENTER to continue'}));
        [x, y] = getpts();
        p = [p;y];
        k = [k;x];
        plot(k,p,'o');
    else
    end
    close(activeFigure)

end