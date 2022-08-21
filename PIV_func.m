function [u,v,Vel,xf,yf,xm,ym] = PIV_func(I1,I2,win,pxs,M,dt,x0,y0,rlen,clen,xmask,ymask)
% PIV Image Processing Function
% (c) Aaron Sequeira, Aniruddha Paranjape, Nikhil Joseph Jose, Jahnavi Raghavaraju, 2021
% (Group 7)

%-------------------------------------------------------------------------
% Function uses PIV image pair as input and executes a cross-correlation
% algorithm. Output is domain Velocity magnitude and x,y components.
%-------------------------------------------------------------------------


%% Load and Prep Images

%Crop Images
I1 = I1(1:end-mod(rlen,win),1:end-mod(clen,win));
I2 = I2(1:end-mod(rlen,win),1:end-mod(clen,win));
[rlen,clen] = size(I1);

% Create Interrogation Windows
n   = rlen/win; m = clen/win;
iw1 = zeros(win,win,n*m);
iw2 = zeros(win,win,n*m);
xf = zeros(n,m);
yf = zeros(n,m);

k = 1;
for i = 1:n
    for j = 1:m
        iw1(:,:,k) = I1((i-1)*win+1:i*win,(j-1)*win+1:j*win);
        iw2(:,:,k) = I2((i-1)*win+1:i*win,(j-1)*win+1:j*win);
        
        xf(i,j) = (j-1)*win + 1 + win/2 - x0;
        yf(i,j) = (i-1)*win + 1 + win/2 - y0;
        
        k=k+1;
    end
end

% Field of View Coordinates
xf = xf*pxs*1e-03/M; yf = -yf*pxs*1e-03/M;

% Airfoil Masking (Physical Domain)
xm = (xmask-x0)*pxs*1e-03/M; ym = -(ymask-y0)*pxs*1e-03/M;  % Mask Coordinates

[in,on] = inpolygon(xf,yf, xm,ym);        % Logical Matrix
inon = in | on;                           % Combine ‘in’ And ‘on’
idx = find(inon(:));                      % Linear Indices Of ‘inon’ Points


%% Correlation and Velocities

% Pre-allocate
u = zeros(1,m*n);
v = zeros(1,m*n);

% Loop Through All Interrogation Windows
for i = 1:m*n

    % Store Current Interrogation Windows in Temp Variables
    t1  = iw1(:,:,i); t2 = iw2(:,:,i);
    
    % Create Lag Mesh
    xy    = -(win-1):(win-1);
    [X,Y] = meshgrid(flip(xy),xy);
    
    % Cross-Correlate
    cor = xcorr2(t1-mean(t1(:)),t2-mean(t2(:)));

    % Find Correlation Peak
    [~, imax] = max(abs(cor(:)));
    [ypeak, xpeak] = ind2sub(size(cor),imax(1));

    % Calculate Pixel Displacements through Gaussian Interpolation
    dx = X(ypeak,xpeak); dy = Y(ypeak,xpeak);

    if (xpeak>1) && (xpeak<2*win-1) && (ypeak>1) && (ypeak<2*win-1)
        dx = dx + (log(cor(ypeak-1,xpeak)) - log(cor(ypeak+1,xpeak)))/...
            (2*log(cor(ypeak-1,xpeak)) - 4*log(cor(ypeak,xpeak)) + 2*log(cor(ypeak+1,xpeak)));

        dy = dy + (log(cor(ypeak,xpeak-1)) - log(cor(ypeak,xpeak+1)))/...
            (2*log(cor(ypeak,xpeak-1)) - 4*log(cor(ypeak,xpeak)) + 2*log(cor(ypeak,xpeak+1)));
    else
        dx = 0; dy = 0;     % Peak on Boundaries 
    end

    % Compute Velocities
    u(i) = (dx*pxs)/(M*dt);
    v(i) = (dy*pxs)/(M*dt);

end


%% Matrix Assembly

% Velocities
u   = reshape(real(u),m,n).';   % Assemble x-Velocity Matrix
v   = reshape(real(v),m,n).';   % Assemble y-Velocity Matrix
u(idx) = nan; v(idx) = nan;     % Mask Airfoil

Vel = sqrt(u.^2 + v.^2);    % Velocity Magnitude


end

