% Custom PIV Image Processing Tool (FMT Lab Exercise)
% (c) Aaron Sequeira, Aniruddha Paranjape, Nikhil Joseph Jose, Jahnavi Raghavaraju, 2021
% (Group 7)

clc
clearvars
tic

%% Processor Inputs

win = 32;       % Window size in px (16,32,64)
pxs = 4.4;      % Pixel Size (Microns)
M   = 0.04348;  % Magnification Factor
dt  = 100;      % Pulse Separation Time (micro-s)

x0 = 125;    % Leading Edge px (x) 
y0 = 530;    % Leading Edge px (y)

f_beg = 10;    % First File
f_end = 11;    % Last File

plt = 1;    % Plot Results? (Yes => 1, No => 0)


%% Read and Process Images

% Image Folder Location
FoldRead = 'Images/';                % Image Folder
load('Masks/Mask_Alpha_15.mat');     % Load Mask

fprintf('Reading Files...\n')
% Read Images
k = 1;
for i = f_beg:f_end
    
    FileRead = ['B' num2str(i,'%05.f') '.tif'];
    
    im  = imread([FoldRead FileRead]);   % Read
    [rlen,~] = size(im);                 % Image Size
    
    % Split and Store Images
    t1 = im(1:rlen/2,:); t2 = im(rlen/2+1:rlen,:);
    [rlen,clen] = size(t1);
    
    % Mask Airfoil (Image Mask)
    bw = roipoly(t1,xmask,ymask);
    t1(bw) = nan; t2(bw) = nan;
    
    I1(:,:,k) = t1; I2(:,:,k) = t2;
    k = k + 1;
    
end
fprintf('Done!\n\n')


%% PIV Processing Function

% Pre-Allocate
U = zeros(fix(rlen/win),fix(clen/win),f_end-f_beg+1);
V = zeros(fix(rlen/win),fix(clen/win),f_end-f_beg+1);

if plt == 1
    figure('Name','RAW')
    figure('Name','Velocity')
end
    
for i = 1:(f_end-f_beg+1)
    fprintf(['Processing Image Pair ' num2str(i) '\n']) 
    
    % Function Call
    [u,v,Vel,xf,yf,xm,ym] = PIV_func(I1(:,:,i),I2(:,:,i),win,pxs,M,dt,x0,y0,rlen,clen,xmask,ymask);
    
    % Plotting
    if plt == 1
        
        % Plot Raw Image
        figure(1)
        set(gcf, 'Position',  [320, 250, 400, 300])
        imagesc(250.*I1(:,:,i))
        colormap('gray')
        hold on
        plot(xmask,ymask,'r','Linewidth',1)
        xlabel('x (px)')
        ylabel('y (px)')
        title('Raw Image')
        
        % Plot Velocity
        figure(2)
        clf
        set(gcf, 'Position',  [750, 250, 450, 300])
        contourf(xf,yf,Vel,25, 'LineStyle','none')
        hold on
        plot(xm,ym,'k','Linewidth',1.5)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title('Velocity (m/s)')
        colorbar
        colormap('jet')
        pause(0.1)
        
    end

    % Saving
    U(:,:,i) = u;
    V(:,:,i) = v;
    save(['Results/Res_' num2str(win) '.mat'],'U','V','xf','yf','xm','ym');

end


%% End
toc

