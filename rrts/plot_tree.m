load nodes.txt
load edges.txt
load optpath.txt

figure
daspect ([1 1 1]);


% Draw the nodes
max_a = max(nodes(:,1));
hold on, plot (nodes (:,1), nodes(:,2), 'g.')


% Draw the edges 
for i = 1 : length (edges) 
    hold on, plot ( [edges(i,1), edges(i,3)], [edges(i,2), edges(i,4)], 'Color', [0.1 0.5 0.7], 'LineWidth', 0.5);
end

%{
% Draw the optimal path
for i = 1:length (optpath)-1
    hold on, plot ([optpath(i,1), optpath(i+1,1)], [optpath(i,2), optpath(i+1,2)], 'r-', 'LineWidth', 4.0);    
end
%}

grid on;
axis([-10 10 -10 10]);
hold on;
alpha(0.5);
rectangle('Position', [-5.5, -7.5, 5, 5], 'Facecolor', 'r');
rectangle('Position', [0.5, -7.5, 5, 5], 'Facecolor', 'r');