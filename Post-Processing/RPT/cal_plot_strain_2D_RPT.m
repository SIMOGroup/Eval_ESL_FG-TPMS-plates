function [output] = cal_plot_strain_2D_RPT(IGA, Plate, plot_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and plot strain field of rectangular plate with 4-Dof 2D RPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'point' = Strain at a specific point in 3D
% 'point_dist' = Strain distribution at a specific point in 2D
% 'FE_surf' = Strain FE surface
% 'FE_vol' = Strain FE volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
norm_method = IGA.result.norm_method(2);

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h; 
strain_vector = ["xx"; "yy"; "xy"; "yz"; "zx"; "zz"];

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
        feature_idx = find(strain_vector == varargin{2});
    case 'FE_surf'
        FE_surf_plot = true;
        feature_idx = find(strain_vector == varargin{1});
    case 'FE_vol'
        FE_vol_plot = true;
        feature_idx = find(strain_vector == varargin{1});
    otherwise
        error('Unknown keyword: %s', plot_type);
end
%% ====== Calculation and plot ======
% --- Strain at a specific point in 3D ---
if point_plot
    point_strain = cal_point_strain_2D_RPT(IGA,Plate,point_coord);
    point_strain_norm = compute_strain_normalization(point_strain,norm_method,Plate);
    disp("> Normalized point strain field (" + strjoin(cellstr(strain_vector.'), ', ') + ") = (" + sprintf('%.4f, ', point_strain_norm) + ")")
end

% --- Strain distribution at a specific point in 2D ---
if point_dist_plot
    FE_strain_curve = {};
    FE_strain_curve.disp_idx = feature_idx;
    FE_strain_curve.npoint = 51;
    FE_strain_curve.z_coords = linspace(-h/2,h/2,FE_strain_curve.npoint);

    point_strain_dist = cal_point_strain_distribution_2D_RPT(IGA,Plate,point_coord,FE_strain_curve.z_coords);
    point_strain_dist_norm = compute_strain_normalization(point_strain_dist,norm_method,Plate);
    
    FE_strain_curve.FE_grid = FE_strain_curve.z_coords;
    FE_strain_curve.FE_feat = point_strain_dist_norm(:,FE_strain_curve.disp_idx);
    FE_strain_curve.data = point_strain_dist_norm;
    plot_FE_curve(FE_strain_curve, 'all')
end

% --- Strain FE surf ---
if FE_surf_plot
    FE_strain_surf = {};
    FE_strain_surf.disp_idx = feature_idx; 
    FE_strain_surf.z_coords = [0];

    FE_strain = cal_strain_field_2D_RPT(IGA,Plate,FE_strain_surf.z_coords);
    FE_strain_norm = compute_strain_normalization(FE_strain,norm_method,Plate);

    FE_strain_surf.FE_grid = IGA.NURBS.FE_grid; 
    FE_strain_surf.FE_feat = FE_strain_norm(:,:,1,FE_strain_surf.disp_idx);
    FE_strain_surf.data = FE_strain_norm;
    plot_FE_surf(FE_strain_surf, 'all')
end

% --- Strain FE vol ---
if FE_vol_plot
    FE_strain_vol = {};
    FE_strain_vol.disp_idx = feature_idx; 
    FE_strain_vol.npoint = 7;
    FE_strain_vol.z_coords = linspace(-h/2,h/2,FE_strain_vol.npoint);

    FE_strain = cal_strain_field_2D_RPT(IGA,Plate,FE_strain_vol.z_coords);
    FE_strain_norm = compute_strain_normalization(FE_strain,norm_method,Plate);
    
    idx_x = round(linspace(1, size(IGA.NURBS.FE_grid,1), 2*FE_strain_vol.npoint+1));
    idx_y = round(linspace(1, size(IGA.NURBS.FE_grid,2), 2*FE_strain_vol.npoint+1));
    FE_strain_vol.FE_2Dgrid = IGA.NURBS.FE_grid(idx_x,idx_y,:);
    FE_strain_vol.FE_3Dgrid.X = repmat(FE_strain_vol.FE_2Dgrid(:,:,1), 1, 1, FE_strain_vol.npoint);
    FE_strain_vol.FE_3Dgrid.Y = repmat(FE_strain_vol.FE_2Dgrid(:,:,2), 1, 1, FE_strain_vol.npoint);
    FE_strain_vol.FE_3Dgrid.Z = repmat(reshape(linspace(-Plate.geo.h/2, Plate.geo.h/2, FE_strain_vol.npoint), 1, 1, []), length(idx_x), length(idx_y), 1);
    FE_strain_vol.FE_feat = FE_strain_norm(idx_x,idx_y,:,FE_strain_vol.disp_idx);
    FE_strain_vol.data = FE_strain_norm;
    plot_FE_vol(FE_strain_vol, 'all')
end

switch plot_type
    case 'point'
        output = point_strain_norm;
    case 'point_dist'
        output = FE_strain_curve;
    case 'FE_surf'
        output = FE_strain_surf;
    case 'FE_vol'
        output = FE_strain_vol;
end
end
