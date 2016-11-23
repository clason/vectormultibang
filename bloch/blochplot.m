function blochplot(d,x)
% BLOCHPLOT plot control and state for Bloch control problem
% BLOCHPLOT(D,X) plots the control X=(u,v) and the corresponding state M. 
% The structure D contains the problem parameters.
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% plot control
u = x(1:d.Nu);  v = x(1+d.Nu:end); 

figure(1);
plot3(d.tdis(1:end-1),u,v); grid on; title('control');
hold on
plot3(d.tdis(1:end-1),repmat(d.ub(1,:),d.Nu,1)',repmat(d.ub(2,:),d.Nu,1)','--');
hold off
xlabel('t');ylabel('u_1'),zlabel('u_2');

%% plot state
% compute state
M = cn_bloch(d,d.M0,u,v,d.w);

% plot
figure(2);
for z = 1:d.Nw
    Mz = squeeze(M(:,z,:));
    plot3(Mz(1,:),Mz(2,:),Mz(3,:)); grid on; title('state');
    hold on;
end
plot3(d.M0(1),d.M0(2),d.M0(3),'o'); text(d.M0(1),d.M0(2),d.M0(3),'M_0');
for z = 1:d.Nw
    Mz = squeeze(M(:,z,end));
    plot3(Mz(1),Mz(2),Mz(3),'x');
    text(Mz(1),Mz(2),Mz(3),'M(T)','HorizontalAlignment','Right');
    plot3(d.Md(1,z),d.Md(2,z),d.Md(3,z),'s');
    text(d.Md(1,z),d.Md(2,z),d.Md(3,z),'M_d');
end
hold off;

xlabel('M_x');ylabel('M_y');zlabel('M_z');
ax = gca; ax.YLim = [-1,1]; ax.XLim = [-1,1]; ax.ZLim = [0,1];
ax.View = [90-37.5,30];

