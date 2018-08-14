
f_c = 0.05;
t = 0:0.05:2;
carrier = sin(pi*t+0).*(f_c*t);

h = figure;
% plot(t,carrier,'LineWidth',2,'Color','b');
axisLimits = [-1,3,-1,1];
axis(axisLimits)
pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])
gif('myfile.gif')

for i = 0:0.5:30
    carrier = sin(pi*t+i).*(f_c*t);
    plot(t,carrier,'LineWidth',2,'Color','b');
    axis(axisLimits)
    pbaspect([(axisLimits(2)-axisLimits(1)) (axisLimits(4)-axisLimits(3)) 1])
    pause(0.1)
    gif
end