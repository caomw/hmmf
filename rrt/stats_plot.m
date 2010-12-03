clc;
clear all;
close all;

load stats_rrtsnobb_20000_obs2.dat;
a = stats_rrtsnobb_20000_obs2;
bins = 100;

times = a(:,3);
costs = a(:,2);

plot_time_hist = 1;
if plot_time_hist,
    lowt = min(times);
    hight = max(times);
    delt = (hight-lowt)/bins;
    
    xt = lowt:delt:hight;
    avg = mean(times);
    stddev = std(times);
    hist(times, xt);
    show = sprintf('Avg: %s [ms]\nStd.dev: %s [ms]', num2str(avg), num2str(stddev));
    legend(show);
    grid on;
    print -dpng -r300 stats_rrtsnobb_20000_obs2.png
end

plot_cost_hist = 0;
if plot_cost_hist,
    figure(2);
    lowc = min(costs);
    highc = max(costs);
    delc = (highc-lowc)/bins;
    
    xc = lowc:delc:highc;
    avg = mean(costs);
    stddev = std(costs);
    hist(costs, xc);
    show = sprintf('Avg: %s [c]\nStd.dev: %s [c]', num2str(avg), num2str(stddev));
    legend(show);
    grid on;
    print -dpng -r300 stats_rrtsnobb_20000_ob3_cost.png
end

plot_avg_cost = 0;
if plot_avg_cost,
    b = load('rrtsbb_avgcost.dat');
    c = load('rrtbb_avgcost.dat');
    plot(b(1245:end,1), b(1245:end,2), 'r-', 'LineWidth', 2);
    hold on;
    plot(c(600:end,1), c(600:end,2), 'b-', 'LineWidth', 2);
    grid on;
    axis([0 40000 55 65]);
    xlabel('Iterations');
    ylabel('Unsmoothened path cost');
    legend('RRT*', 'RRT')
end