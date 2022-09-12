%% Averaged edge mode dynamics for different initial times
clear all
close all
home 

cd upper_edge_data                                                          % points to th data directory
files = dir('*.mat');                                                       % loads all .mat files
for j=1:length(files)                                                       % loop through the files
    datafile = load(files(j).name, 'evolution_matx');                       % loads the data from each file  
    t0_dynamics = datafile.evolution_matx;                                  % selects the variable from the file
    if j==1
        av_dynamics = t0_dynamics;
    else
    av_dynamics = av_dynamics + t0_dynamics;
    end
end

N =length(files);
av_dynamics = av_dynamics/ N;                                               % Number average

figure(1)
s=surf(av_dynamics);
set(gcf,'Position',[700 200 1000 700])
box on
view(2)
c=jet;
colormap(c);
caxis([0 0.06])
colorbar
xlabel('Edge Site','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel('Time [T]','FontSize',25,'FontWeight','bold','Interpret','latex')
title('Upper edge dynamics, $n_{t_0}=20$', 'FontSize',25,'FontWeight','bold','Interpret','latex')
h = colorbar;
set(get(h,'label'),'string','$\vert \Psi \vert^2$ [arbitrary units]','FontSize',20,'FontWeight','bold','Interpret','latex');
xlim([1 length(av_dynamics(1, :))])
ylim([1 length(av_dynamics(:, 1))])
s.EdgeColor = 'flat';
ax = gca;
ax.FontSize =25; 
ax.XTick = [0:10:length(av_dynamics(1, :))];