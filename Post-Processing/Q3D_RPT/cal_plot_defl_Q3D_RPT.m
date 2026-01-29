function [output] = cal_plot_defl_Q3D_RPT(IGA, Plate, plot_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and plot deflection (transverse displacement) of rectangular plate with 4-Dof Quasi-3D RPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'cen_defl' = Deflection at the central point
% 'NURBS_defl' = Deflection NURBS mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
norm_method = IGA.result.norm_method(1);

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h; 

%% Check plot options
cen_defl = false; NURBS_defl = false;
switch plot_type
    case 'cen_defl'
        cen_defl = true;
    case 'NURBS_defl'
        NURBS_defl = true;
    otherwise
        error('Unknown keyword: %s', plot_type);
end

%% ====== Calculation and plot ======
% --- Deflection at the central point ---
if cen_defl
    cen_point_disp = cal_point_displacement_Q3D_RPT(IGA,Plate,[L/2,W/2,0]);
    cen_defl_norm = compute_displacement_normalization(cen_point_disp(3),norm_method,Plate);
    disp("> Normalized central deflection = " + sprintf('%.4f, ', cen_defl_norm))
end

% --- Deflection NURBS mesh ---
if NURBS_defl
    NURBS_defl = cal_deflection_Q3D_RPT(IGA,Plate);
    NURBS_defl_norm = compute_displacement_normalization(NURBS_defl,norm_method,Plate);
    NURBS_defl_plot = IGA.NURBS; plot_size_factor = 1;
    NURBS_defl_plot.CP(:,:,3) = NURBS_defl_plot.CP(:,:,3) + plot_size_factor * NURBS_defl_norm;
    NURBS_defl_plot.data = NURBS_defl_norm;
    NURBS_defl_plot = gen_Iee_Ine_surf(NURBS_defl_plot); NURBS_defl_plot = gen_FE_approx_surf(NURBS_defl_plot);
    plot_NURBS_surf(NURBS_defl_plot, 'physical_surf', 'plot_physical_element', 'plot_control_net')
end

switch plot_type
    case 'cen_defl'
        output = cen_defl_norm;
    case 'NURBS_defl'
        output = NURBS_defl_plot;
end

end
