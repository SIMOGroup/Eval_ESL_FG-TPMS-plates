function [output] = cal_plot_stress_2D_GPT(IGA, Plate, plot_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and plot stress field of rectangular plate with 5-Dof 2D GPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'point' = Stress at a specific point in 3D
% 'point_dist' = Stress distribution at a specific point in 2D
% 'FE_surf' = Stress FE surface
% 'FE_vol' = Stress FE volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
norm_method = IGA.result.norm_method(3);

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h; 
stress_vector = ["xx"; "yy"; "xy"; "yz"; "zx"; "zz"];

%% Check plot options
point_plot = false; point_dist_plot = false;
FE_surf_plot = false; FE_vol_plot = false;
switch plot_type
    case 'point'
        point_plot = true;
        point_coord = varargin{1} .* [L, W, h];
    case 'point_dist'
        point_dist_plot = true;
        point_coord = varargin{1} .* [L, W];
        feature_idx = find(stress_vector == varargin{2});
    case 'FE_surf'
        FE_surf_plot = true;
        feature_idx = find(stress_vector == varargin{1});
    case 'FE_vol'
        FE_vol_plot = true;
        feature_idx = find(stress_vector == varargin{1});
    otherwise
        error('Unknown keyword: %s', plot_type);
end
%% ====== Calculation and plot ======
% --- Stress at a specific point in 3D ---
if point_plot
    point_stress = cal_point_stress_2D_GPT(IGA,Plate,point_coord);
    point_stress_norm = compute_stress_normalization(point_stress,norm_method,Plate);
    disp("> Normalized point stress field (" + strjoin(cellstr(stress_vector.'), ', ') + ") = (" + sprintf('%.4f, ', point_stress_norm) + ")")
end

% --- Stress distribution at a specific point in 2D ---
if point_dist_plot
    FE_stress_curve = {};
    FE_stress_curve.disp_idx = feature_idx;
    FE_stress_curve.npoint = 51;
    FE_stress_curve.z_coords = linspace(-h/2,h/2,FE_stress_curve.npoint);

    point_stress_dist = cal_point_stress_distribution_2D_GPT(IGA,Plate,point_coord,FE_stress_curve.z_coords);
    point_stress_dist_norm = compute_stress_normalization(point_stress_dist,norm_method,Plate);

    FE_stress_curve.FE_grid = FE_stress_curve.z_coords;
    FE_stress_curve.FE_feat = point_stress_dist_norm(:,FE_stress_curve.disp_idx);
    FE_stress_curve.data = point_stress_dist_norm;
    plot_FE_curve(FE_stress_curve, 'all')
end

% --- Stress FE surf ---
if FE_surf_plot
    FE_stress_surf = {};
    FE_stress_surf.disp_idx = feature_idx; 
    FE_stress_surf.z_coords = [0];

    FE_stress = cal_stress_field_2D_GPT(IGA,Plate,FE_stress_surf.z_coords);
    FE_stress_norm = compute_stress_normalization(FE_stress,norm_method,Plate);

    FE_stress_surf.FE_grid = IGA.NURBS.FE_grid; 
    FE_stress_surf.FE_feat = FE_stress_norm(:,:,1,FE_stress_surf.disp_idx);
    FE_stress_surf.data = FE_stress_norm;
    plot_FE_surf(FE_stress_surf, 'all')
end

% --- Stress FE vol ---
if FE_vol_plot
    FE_stress_vol = {};
    FE_stress_vol.disp_idx = feature_idx; 
    FE_stress_vol.npoint = 7;
    FE_stress_vol.z_coords = linspace(-h/2,h/2,FE_stress_vol.npoint);

    FE_stress = cal_stress_field_2D_GPT(IGA,Plate,FE_stress_vol.z_coords);
    FE_stress_norm = compute_stress_normalization(FE_stress,norm_method,Plate);
    
    idx_x = round(linspace(1, size(IGA.NURBS.FE_grid,1), 2*FE_stress_vol.npoint+1));
    idx_y = round(linspace(1, size(IGA.NURBS.FE_grid,2), 2*FE_stress_vol.npoint+1));
    FE_stress_vol.FE_2Dgrid = IGA.NURBS.FE_grid(idx_x,idx_y,:);
    FE_stress_vol.FE_3Dgrid.X = repmat(FE_stress_vol.FE_2Dgrid(:,:,1), 1, 1, FE_stress_vol.npoint);
    FE_stress_vol.FE_3Dgrid.Y = repmat(FE_stress_vol.FE_2Dgrid(:,:,2), 1, 1, FE_stress_vol.npoint);
    FE_stress_vol.FE_3Dgrid.Z = repmat(reshape(linspace(-Plate.geo.h/2, Plate.geo.h/2, FE_stress_vol.npoint), 1, 1, []), length(idx_x), length(idx_y), 1);
    FE_stress_vol.FE_feat = FE_stress_norm(idx_x,idx_y,:,FE_stress_vol.disp_idx);
    FE_stress_vol.data = FE_stress_norm;
    plot_FE_vol(FE_stress_vol, 'all')
end

switch plot_type
    case 'point'
        output = point_stress_norm;
    case 'point_dist'
        output = FE_stress_curve;
    case 'FE_surf'
        output = FE_stress_surf;
    case 'FE_vol'
        output = FE_stress_vol;
end
end
