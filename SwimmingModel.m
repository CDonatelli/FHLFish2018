function SwimmingModel(nameString)

    % Add the amp factor from R here
    AmpFactor = 0.083;
    Intercept = -0.004;
    t = 0:0.05:2;
    axisLimits = [-1,3,-1,1];

    h = figure;
    axis(axisLimits)
    pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])

    % This is the code to generate the background image
    for i = 0:0.75:2*pi
        carrier = sin(pi*t+i).*(AmpFactor*t);
        plot(t,carrier,'LineWidth',1,'Color','k');
        axis(axisLimits)
        pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])
        hold on
    end
    print(h, '-djpeg', [nameString,'background']);
    img = imread('background.jpg');

    % This is the code to generate the gif
    h = figure;
    carrier = sin(pi*t+0).*(AmpFactor*t);
    plot(t,carrier,'LineWidth',2,'Color','b');
    axis(axisLimits)
    pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])
    gif([nameString,'.gif'])
    hold off
    for i = 0.5:0.5:30
        carrier = sin(pi*t+i).*(AmpFactor*t);
        plot(t,carrier,'LineWidth',2,'Color','b');
        axis(axisLimits)
        pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])
        pause(0.1)
        gif
    end

end