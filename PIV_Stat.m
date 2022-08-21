% PIV Statistics Tool (FMT Lab Exercise)
% (c) Aaron Sequeira, Aniruddha Paranjape, Nikhil Joseph Jose, Jahnavi Raghavaraju, 2021
% (Group 7)

clc
clearvars
tic

%% Read and Process Data

win = 32;   % Result Window Size (px)

% Load Results
load(['Results/Res_' num2str(win) '.mat']);     % Load Results

% Calculate Statistics
fprintf('Calculating Statistics...\n')
U_mean = mean(abs(U),3);
U_std  = std(abs(U),0,3);
V_mean = mean(abs(V),3);
V_std  = std(abs(V),0,3); 

Vel_mean = sqrt(U_mean.^2 + V_mean.^2);
Vel_std  = sqrt(U_std.^2 + V_std.^2);

% Calculate Profile at x = 1.2c
y_wk = linspace(-40,40,100);
x_wk = 125*ones(size(y_wk));

Velm_wk   = interp2(xf,yf,U_mean,x_wk,y_wk);
Velstd_wk = interp2(xf,yf,U_std,x_wk,y_wk);

fprintf('Done!\n')

%% Plotting Routine

% Velocity Contours
figure('Name','Velocity')
clf
set(gcf, 'Position',  [240, 400, 1050, 320])
scale = 0.05;

subplot(1,2,1)
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

contourf(xf,yf,U_mean,25, 'LineStyle','none')
hold on
plot(xm,ym,'k','Linewidth',1.5)
xlabel('x (mm)')
ylabel('y (mm)')
title('Mean Velocity Magnitude (m/s)')
colorbar
colormap('jet')

subplot(1,2,2)
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

contourf(xf,yf,Vel_std,25, 'LineStyle','none')
hold on
plot(xm,ym,'k','Linewidth',1.5)
xlabel('x (mm)')
ylabel('y (mm)')
title('RMS Fluctuation Magnitude (m/s)')
colorbar
colormap('jet')

% Velocity Profiles at x = 1.2c
figure('Name','Velocity')
set(gcf, 'Position',  [380, 100, 750, 280])
clf

subplot(1,2,1)
plot(Velm_wk,y_wk,'k','Linewidth',1.5)
grid on
xlabel('Velocity (m/s)')
ylabel('y (mm)')
title('Mean Velocity Profile at x = 1.2c')

subplot(1,2,2)
plot(Velstd_wk,y_wk,'r','Linewidth',1.5)
grid on
xlabel('Velocity (m/s)')
ylabel('y (mm)')
title('Fluctuation Profile at x = 1.2c')


%% End
toc

