function plot_FE_curve(FE, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot FE curve with options %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'all' = 'plot_physical_curve', 'plot_color_bar'
% 'physical_curve' = 'plot_physical_curve'
% 'color_bar' = 'plot_color_bar'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from FE
FE_grid = FE.FE_grid;
FE_feat = FE.FE_feat;

%% Check plot options
plot_physical_curve = false; plot_color_bar = false;

for i = 1:length(varargin)
    switch lower(varargin{i})
        case 'all'
            plot_physical_curve = true;
            plot_color_bar = true;
            
        case 'physical_curve'
            plot_physical_curve = true;
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
 
title('\bf{FE Curve Visualization}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
xlabel('Feature', 'Interpreter','latex', 'fontsize', 14);
ylabel('z', 'Interpreter','latex', 'fontsize', 14);

set(gca, 'TickLabelInterpreter','latex', 'TickLength', [0.02 0.02], 'TickDir','in', 'fontsize', 12)
set(gca, 'box','on', 'GridLineStyle', 'none', 'GridAlpha', 0.2)
set(gca, 'Layer', 'Top')
view(2)
% axis equal

set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 29.7]);
set(gcf, 'PaperOrientation','landscape');

% --- Physical curve ---
if plot_physical_curve
    if plot_color_bar
        hold on
        surface([FE_feat(:), FE_feat(:)]', [FE_grid(:), FE_grid(:)]', zeros(2, numel(FE_grid(:))), [FE_feat(:), FE_feat(:)]', ...
                'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 2);
        scatter(FE_feat(:), FE_grid(:), 50, FE_feat(:), 'filled', 'MarkerEdgeColor', 'k');
    else
        plot(FE_feat(:), FE_grid(:), 'LineStyle','-', 'LineWidth',1.5, 'Color','b', ...
                'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
    end

    uistack(gca, 'top');
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
