function [output] = cal_plot_disp_2D_TSDPT(IGA, Plate, plot_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and plot displacement field of rectangular plate with 5-Dof 2D Reddy's TSDPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'point' = Displacement at a specific point in 3D
% 'point_dist' = Displacement distribution at a specific point in 2D
% 'FE_surf' = Displacement FE surface
% 'FE_vol' = Displacement FE volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
norm_method = IGA.result.norm_method(1);

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h; 
disp_vector = ["u"; "v"; "w"];

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
        feature_idx = find(disp_vector == varargin{2});
    case 'FE_surf'
        FE_surf_plot = true;
        feature_idx = find(disp_vector == varargin{1});
    case 'FE_vol'
        FE_vol_plot = true;
        feature_idx = find(disp_vector == varargin{1});
    otherwise
        error('Unknown keyword: %s', plot_type);
end
%% ====== Calculation and plot ======
% --- Displacement at a specific point in 3D ---
if point_plot
    point_disp = cal_point_displacement_2D_TSDPT(IGA,Plate,point_coord);
    point_disp_norm = compute_displacement_normalization(point_disp,norm_method,Plate);
    disp("> Normalized point displacement field (" + strjoin(cellstr(disp_vector.'), ', ') + ") = (" + sprintf('%.4f, ', point_disp_norm) + ")")
end

% --- Displacement distribution at a specific point in 2D ---
if point_dist_plot
    FE_disp_curve = {};
    FE_disp_curve.disp_idx = feature_idx;
    FE_disp_curve.npoint = 51;
    FE_disp_curve.z_coords = linspace(-h/2,h/2,FE_disp_curve.npoint);

    point_disp_dist = cal_point_displacement_distribution_2D_TSDPT(IGA,Plate,point_coord,FE_disp_curve.z_coords);
    point_disp_dist_norm = compute_displacement_normalization(point_disp_dist,norm_method,Plate);
    
    FE_disp_curve.FE_grid = FE_disp_curve.z_coords;
    FE_disp_curve.FE_feat = point_disp_dist_norm(:,FE_disp_curve.disp_idx);
    FE_disp_curve.data = point_disp_dist_norm;
    plot_FE_curve(FE_disp_curve, 'all')
end

% --- Displacement FE surf ---
if FE_surf_plot
    FE_disp_surf = {};
    FE_disp_surf.disp_idx = feature_idx; 
    FE_disp_surf.z_coords = [0];

    FE_disp = cal_displacement_field_2D_TSDPT(IGA,Plate,FE_disp_surf.z_coords);
    FE_disp_norm = compute_displacement_normalization(FE_disp,norm_method,Plate);

    FE_disp_surf.FE_grid = IGA.NURBS.FE_grid; 
    FE_disp_surf.FE_feat = FE_disp_norm(:,:,1,FE_disp_surf.disp_idx);
    FE_disp_surf.data = FE_disp_norm;
    plot_FE_surf(FE_disp_surf, 'all')
end

% --- Displacement FE vol ---
if FE_vol_plot
    FE_disp_vol = {};
    FE_disp_vol.disp_idx = feature_idx; 
    FE_disp_vol.npoint = 7;
    FE_disp_vol.z_coords = linspace(-h/2,h/2,FE_disp_vol.npoint);

    FE_disp = cal_displacement_field_2D_TSDPT(IGA,Plate,FE_disp_vol.z_coords);
    FE_disp_norm = compute_displacement_normalization(FE_disp,norm_method,Plate);

    idx_x = round(linspace(1, size(IGA.NURBS.FE_grid,1), 2*FE_disp_vol.npoint+1));
    idx_y = round(linspace(1, size(IGA.NURBS.FE_grid,2), 2*FE_disp_vol.npoint+1));
    FE_disp_vol.FE_2Dgrid = IGA.NURBS.FE_grid(idx_x,idx_y,:);
    FE_disp_vol.FE_3Dgrid.X = repmat(FE_disp_vol.FE_2Dgrid(:,:,1), 1, 1, FE_disp_vol.npoint);
    FE_disp_vol.FE_3Dgrid.Y = repmat(FE_disp_vol.FE_2Dgrid(:,:,2), 1, 1, FE_disp_vol.npoint);
    FE_disp_vol.FE_3Dgrid.Z = repmat(reshape(linspace(-Plate.geo.h/2, Plate.geo.h/2, FE_disp_vol.npoint), 1, 1, []), length(idx_x), length(idx_y), 1);
    FE_disp_vol.FE_feat = FE_disp_norm(idx_x,idx_y,:,FE_disp_vol.disp_idx);
    FE_disp_vol.data = FE_disp_norm;
    plot_FE_vol(FE_disp_vol, 'all')
end

switch plot_type
    case 'point'
        output = point_disp_norm;
    case 'point_dist'
        output = FE_disp_curve;
    case 'FE_surf'
        output = FE_disp_surf;
    case 'FE_vol'
        output = FE_disp_vol;
end
end
