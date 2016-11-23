function phaseplot(xx,yy,u,umag)
% PHASEPLOT plot 2d vector field via color coding and scaled arrow grid
% PHASEPLOT(XX,YY,U,UMAG) plots the 2d vector field U defined on the mesh 
% grid [XX,YY] as a color-coded plot (hue indicates angle, saturation 
% indicates magnitude) and as a superimposed arrow grid (arrow direction 
% indicates angle, arrow length indicates magnitude relative to reference 
% magnitude UMAG). A color wheel is also drawn as legend.
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/29487
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

N = size(xx,1);
u1 = reshape(u(:,1),N,N);  u2 = reshape(u(:,2),N,N);

% downsampling and scaling for arrow plots
dsamp = 16; ind = dsamp/2:dsamp:N;
dscal = 10*max(umag);

figure(1);  subplot(2,2,[1,3]); 
% color-coded plot
vfcolor(u1,u2,max(umag),xx(1,[1 end]),yy([1 end],1)'); 
hold on;
% quiver plot 
quiver(xx(ind,ind),yy(ind,ind),u1(ind,ind)/dscal,u2(ind,ind)/dscal,'k','AutoScale','off');
axis equal; axis off; hold off;
title('control');

% plot color wheel
h = subplot(2,2,2);
side = 256;
[xf,yf] = meshgrid(ceil(-side/2:side/2),ceil(-side/2:side/2));
nf = sqrt(xf.^2+yf.^2);
corners = find(nf >= side/2);
xf(corners) = 0;   yf(corners) = 0;
vfcolor(xf,yf,max(nf(:)),xx(1,[1 end]),yy([1 end],1)'); axis square
set(h,'XTick',0.5,'XTickLabel','-\pi/2','YTick',1,'YTickLabel','\pi');
title({'hue = angle','saturation = magnitude'})
end

%%
function img = vfcolor(xf,yf,maxs,xx,yy)
% VFCOLOR convert 2d vector field to hue-saturation color image and plot
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/29487

% rotate 45 degrees
xt = (xf+yf)/sqrt(2);
yt = (xf-yf)/sqrt(2);
xf = xt; yf = yt;

% convert vector field to color plot
[nrow,ncol] = size(xf);
h = mod(8*(atan2(yf,xf)/(2*pi)),8);
h2 = mod(ceil(h),8)+1;
h1 = mod(h2-2,8)+1;
hr = ceil(h)-h;
s = sqrt(xf.^2+yf.^2);
s = s./(maxs+eps);
% husl color palette (generated with seaborn http://seaborn.pydata.org)
cm = [0.9677975592919913, 0.44127456009157356, 0.5358103155058701;
      0.8087954113106306, 0.5634700050056693, 0.19502642696727285;
      0.5920891529639701, 0.6418467016378244, 0.1935069134991043;
      0.19783576093349015, 0.6955516966063037, 0.3995301037444499;
      0.21044753832183283, 0.6773105080456748, 0.6433941168468681;
      0.22335772267769388, 0.6565792317435265, 0.8171355503265633;
      0.6423044349219739, 0.5497680051256467, 0.9582651433656727;
      0.9603888539940703, 0.3814317878772117, 0.8683117650835491];
c = reshape(cm,8,1,3);
img1 = zeros(nrow,ncol,3);
img1(:) = c(h1,1,:);
img2 = zeros(nrow,ncol,3);
img2(:) = c(h2,1,:);
rhr = repmat(hr,[1,1,3]);
img = rhr.*img1+(1-rhr).*img2;
rs = repmat(s,[1,1,3]);
img = rs.*img+(1-rs);

% plot
imagesc(xx,yy,img);
axis xy equal tight
end