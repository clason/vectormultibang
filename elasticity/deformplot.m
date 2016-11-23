function deformplot(xx,yy,y,z)
% DEFORMPLOT plot 2d vector field as grid deformation
% DEFORMPLOT(XX,YY,Y,Z) plots the grid specified by its vertex coordinates 
% [XX,YY] deformed according to the two vector fields Y (the state, in red) 
% and Z (the target, in black)
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

% downsample mesh for plotting
N = size(xx,1);
dsamp = 16;
ind = 1:dsamp:N;
z = reshape(z,[N N 2])+cat(3,xx,yy);
y = reshape(y,[N N 2])+cat(3,xx,yy);

figure(2);  
% plot deformed mesh for target z (gray)
plot(z(:,ind,1),z(:,ind,2),'k',z(ind,:,1)',z(ind,:,2)','k');
hold on;
% plot deformed mesh for state y (red)
plot(y(1:end,ind,1),y(:,ind,2),'r',y(ind,:,1)',y(ind,:,2)','r');
% plot clamped boundary
aux = linspace(xx(1),xx(end),20);
plot([aux;aux-yy(1)-.1],[0*aux;0*aux-yy(1)-.1],'k');
axis equal; axis off; hold off;
title('achieved deformation (red) vs. target deformation (black)');