function plot_FE_surf(FE, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot FE surf with options %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'all' = 'plot_physical_shape', 'plot_physical_surf', 'plot_physical_edge', 'plot_color_bar'
% 'physical_shape' = 'plot_physical_shape'
% 'physical_surf' = 'plot_physical_shape', 'plot_physical_surf'
% 'physical_edge' = 'plot_physical_edge'
% 'color_bar' = 'plot_color_bar'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from FE
FE_grid = FE.FE_grid; FE_feat = FE.FE_feat;

%% Check plot options
plot_physical_shape = false; plot_physical_surf = false;
plot_physical_edge = false; plot_color_bar = false;

for i = 1:length(varargin)
    switch lower(varargin{i})
        case 'all'
            plot_physical_shape = true;
            plot_physical_surf = true;
            plot_physical_edge = true;
            plot_color_bar = true;
            
        case 'physical_shape'
            plot_physical_shape = true;
        case 'physical_surf'
            plot_physical_shape = true;
            plot_physical_surf = true;
        case 'physical_edge'
            plot_physical_edge = true;
        case 'color_bar'
            plot_color_bar = true;
        otherwise
            error('Unknown keyword: %s', varargin{i});
    end
end
    
%% ====== Plot FE mesh ======
% --- Initializate figure ---
figure('Color','w', 'Units','normalized', 'Outerposition',[0 0 0.6 0.8]);
set(axes, 'Position', [0.1, 0.1, 0.7, 0.8]);
hold on
 
title('\bf{FE Surf Visualization}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
xlabel('x', 'Interpreter','latex', 'fontsize', 14);
ylabel('y', 'Interpreter','latex', 'fontsize', 14);
zlabel('Feature', 'Interpreter','latex', 'fontsize', 14);

set(gca, 'TickLabelInterpreter','latex', 'TickLength', [0.02 0.02], 'TickDir','in', 'fontsize', 12)
set(gca, 'box','on', 'GridLineStyle', 'none', 'GridAlpha', 0.2)
set(gca, 'Layer', 'Top')
view(3)
axis equal

set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 29.7]);
set(gcf, 'PaperOrientation','landscape');

% --- Physical shape ---
if plot_physical_shape
    plot3(FE_grid(:,1,1), FE_grid(:,1,2), FE_feat(:,1), 'Color','k', 'LineWidth',3);  % edge_x_lower
    plot3(FE_grid(:,end,1), FE_grid(:,end,2), FE_feat(:,end), 'Color','k', 'LineWidth',3);  % edge_x_upper
    plot3(FE_grid(1,:,1), FE_grid(1,:,2), FE_feat(1,:), 'Color','k', 'LineWidth',3);  % edge_y_lower
    plot3(FE_grid(end,:,1), FE_grid(end,:,2), FE_feat(end,:), 'Color','k', 'LineWidth',3);  % edge_y_upper
    
    uistack(gca, 'top');
end

% --- Physical surf ---
if plot_physical_surf
    surf(FE_grid(:,:,1), FE_grid(:,:,2), FE_feat(:,:), 'FaceColor','interp', 'FaceAlpha',0.8, ...
                                                         'EdgeColor', 'interp', 'EdgeAlpha', 0.5, 'LineStyle', 'none', 'LineWidth', 1, ...
                                                         'Marker', 'none')
end

% --- Physical edge ---
if plot_physical_edge
    idx_x = round(linspace(1, size(FE_grid,1), 5));
    for i_x = 2:length(idx_x)-1
        plot3(FE_grid(idx_x(i_x),:,1), FE_grid(idx_x(i_x),:,2), FE_feat(idx_x(i_x),:), 'LineStyle', '--', 'Color','r', 'LineWidth',1.5);  % edge_y
    end

    idx_y = round(linspace(1, size(FE_grid,2), 5));
    for i_y = 2:length(idx_y)-1
        plot3(FE_grid(:,idx_y(i_y),1), FE_grid(:,idx_y(i_y),2), FE_feat(:,idx_y(i_y)), 'LineStyle', '--', 'Color','r', 'LineWidth',1.5);  % edge_x
    end
end

% --- Color bar ---
if plot_color_bar
    cbar = colorbar('Position', [0.85 0.2 0.03 0.6], 'Ticks', linspace(min(FE_feat(:)), max(FE_feat(:)), 7), 'TickLabelInterpreter', 'latex');
    cbar.Ruler.TickLabelFormat = char(["%0.2f"]);
    colormap(parula(1024));
    clim([min(FE_feat(:)), max(FE_feat(:))]);

    cbar.FontSize = 14;
    cbar.Label.String = 'Color indication';
    cbar.Label.Interpreter = "latex";
    cbar.Label.FontSize = 14;
end

end
