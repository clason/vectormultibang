function plot_flow(u, param)
% PLOT_FLOW plots multimaterial branched transport flow
% PLOT_FLOW(U,PARAM) plots the material flux U on the branched network described by PARAMs.
%
% April 12, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% show network
for k = 1:length(param.edges)
    plot(param.nodes(param.edges(k,:),1),param.nodes(param.edges(k,:),2),'color',.9*[1 1 1]);
    hold on;
end
axis equal;
axis off

%% show fluxes
maxFlow = max(abs(u(:)));
color = [.8;.4;0;.2]*[1 1 1];
thickness = 2*[3 2 1 4];
if param.M == 3
    order = [1 2 3];
else
    order = [4 1 2 3];
end
for j = order
    for k = 1:length(param.edges)
        if abs(u(j,k)) > .1 * maxFlow
            plot(param.nodes(param.edges(k,:),1),param.nodes(param.edges(k,:),2),'-o',...
                'color',(1-abs(u(j,k))/maxFlow)*[1 1 1]+abs(u(j,k))/maxFlow*color(j,:),...
                'MarkerFaceColor',(1-abs(u(j,k))/maxFlow)*[1 1 1]+abs(u(j,k))/maxFlow*color(j,:),...
                'LineWidth',thickness(j),...
                'MarkerSize',thickness(j)/4);
        end
    end
end
hold off;
end