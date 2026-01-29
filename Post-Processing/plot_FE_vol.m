function plot_FE_vol(FE, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot FE volume with options %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'all' = 'plot_physical_shape', 'plot_physical_vol', 'plot_color_bar'
% 'physical_shape' = 'plot_physical_shape'
% 'physical_vol' = 'plot_physical_shape', 'plot_physical_vol'
% 'color_bar' = 'plot_color_bar'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from FE
FE_2Dgrid = FE.FE_2Dgrid; FE_3Dgrid = FE.FE_3Dgrid; 
FE_feat = FE.FE_feat;

%% Check plot options
plot_physical_shape = false; plot_physical_vol = false; plot_color_bar = false;

for i = 1:length(varargin)
    switch lower(varargin{i})
        case 'all'
            plot_physical_shape = true;
            plot_physical_vol = true;
            plot_color_bar = true;
            
        case 'physical_shape'
            plot_physical_shape = true;
        case 'physical_surf'
            plot_physical_shape = true;
            plot_physical_vol = true;
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
 
title('\bf{FE Volume Visualization}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
xlabel('x', 'Interpreter','latex', 'fontsize', 14);
ylabel('y', 'Interpreter','latex', 'fontsize', 14);
zlabel('z', 'Interpreter','latex', 'fontsize', 14);

set(gca, 'TickLabelInterpreter','latex', 'TickLength', [0.02 0.02], 'TickDir','in', 'fontsize', 12)
set(gca, 'box','on', 'GridLineStyle', 'none', 'GridAlpha', 0.2)
set(gca, 'Layer', 'Top')
view(3)
axis equal

set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 29.7]);
set(gcf, 'PaperOrientation','landscape');

% --- Physical shape ---
if plot_physical_shape
    plot3(reshape(FE_3Dgrid.X(:,1,1),1,[]), reshape(FE_3Dgrid.Y(:,1,1),1,[]), reshape(FE_3Dgrid.Z(:,1,1),1,[]), 'Color','k', 'LineWidth',2);  % edge_x
    plot3(reshape(FE_3Dgrid.X(:,1,end),1,[]), reshape(FE_3Dgrid.Y(:,1,end),1,[]), reshape(FE_3Dgrid.Z(:,1,end),1,[]), 'Color','k', 'LineWidth',2);  % edge_x
    plot3(reshape(FE_3Dgrid.X(:,end,1),1,[]), reshape(FE_3Dgrid.Y(:,end,1),1,[]), reshape(FE_3Dgrid.Z(:,end,1),1,[]), 'Color','k', 'LineWidth',2);  % edge_x
    plot3(reshape(FE_3Dgrid.X(:,end,end),1,[]), reshape(FE_3Dgrid.Y(:,end,end),1,[]), reshape(FE_3Dgrid.Z(:,end,end),1,[]), 'Color','k', 'LineWidth',2);  % edge_x

    plot3(reshape(FE_3Dgrid.X(1,:,1),1,[]), reshape(FE_3Dgrid.Y(1,:,1),1,[]), reshape(FE_3Dgrid.Z(1,:,1),1,[]), 'Color','k', 'LineWidth',2);  % edge_y
    plot3(reshape(FE_3Dgrid.X(1,:,end),1,[]), reshape(FE_3Dgrid.Y(1,:,end),1,[]), reshape(FE_3Dgrid.Z(1,:,end),1,[]), 'Color','k', 'LineWidth',2);  % edge_y
    plot3(reshape(FE_3Dgrid.X(end,:,1),1,[]), reshape(FE_3Dgrid.Y(end,:,1),1,[]), reshape(FE_3Dgrid.Z(end,:,1),1,[]), 'Color','k', 'LineWidth',2);  % edge_y
    plot3(reshape(FE_3Dgrid.X(end,:,end),1,[]), reshape(FE_3Dgrid.Y(end,:,end),1,[]), reshape(FE_3Dgrid.Z(end,:,end),1,[]), 'Color','k', 'LineWidth',2);  % edge_y
    
    plot3(reshape(FE_3Dgrid.X(1,1,:),1,[]), reshape(FE_3Dgrid.Y(1,1,:),1,[]), reshape(FE_3Dgrid.Z(1,1,:),1,[]), 'Color','k', 'LineWidth',2);  % edge_z
    plot3(reshape(FE_3Dgrid.X(1,end,:),1,[]), reshape(FE_3Dgrid.Y(1,end,:),1,[]), reshape(FE_3Dgrid.Z(1,end,:),1,[]), 'Color','k', 'LineWidth',2);  % edge_z
    plot3(reshape(FE_3Dgrid.X(end,1,:),1,[]), reshape(FE_3Dgrid.Y(end,1,:),1,[]), reshape(FE_3Dgrid.Z(end,1,:),1,[]), 'Color','k', 'LineWidth',2);  % edge_z
    plot3(reshape(FE_3Dgrid.X(end,end,:),1,[]), reshape(FE_3Dgrid.Y(end,end,:),1,[]), reshape(FE_3Dgrid.Z(end,end,:),1,[]), 'Color','k', 'LineWidth',2);  % edge_z

    uistack(gca, 'top');
end

% --- Physical vol ---
if plot_physical_vol
    scatter3(FE_3Dgrid.X(:), FE_3Dgrid.Y(:), FE_3Dgrid.Z(:), 12, FE_feat(:), ...
                 'Marker','o', 'MarkerEdgeColor','flat', 'MarkerFaceColor','flat', 'MarkerEdgeAlpha',1, 'MarkerFaceAlpha', 1)

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
