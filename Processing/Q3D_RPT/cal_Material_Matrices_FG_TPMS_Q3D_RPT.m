function [D] = cal_Material_Matrices_FG_TPMS_Q3D_RPT(Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate material matrices for FG-TPMS for 4-Dof Quasi-3D RPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
h = Plate.geo.h;
shear_func = Plate.theory.shear_func; n_iterated = Plate.theory.iterated_order; stretch_func = Plate.theory.stretch_func;

%% ===== Initial material matrices =====
for i = 3:-1:1
    for j = 3:-1:1
        D.Db{i,j} = zeros(3,3);
    end
end
D.Ds = zeros(2,2);
for i = 3:-1:1
    D.Dbst{i} = zeros(3,1);
end
D.Dst = zeros(1,1);

%% ===== Gauss integration =====
[Pg, Wg] = gen_Gauss_point(1,-1,20) ;

%% ===== Calculate material matrices =====
z_lower = -h/2; z_upper = h/2;

for iGauss = 1:size(Wg,1)  % Loop over the integration points
    % --- Gauss points & weights ---
    gpt = Pg(iGauss);
    gwt = Wg(iGauss);
    
    % --- Map the point to global ---
    gpt = (z_upper-z_lower)/2*gpt + (z_upper+z_lower)/2;
    gwt = (z_upper-z_lower)/2*gwt;
   
    % --- Consitutive matrix ---
    [Q] = cal_point_constitutive_FG_TPMS_3D(Plate, gpt);
    Qb = Q.Qb; Qs = Q.Qs; Qbst = Q.Qbst; Qst = Q.Qst;
    
    % --- Shear deformation function ---
    [fz, dfz, ~] = compute_iterated_shear_deformation_function(gpt,h,shear_func,n_iterated);
    
    % --- Stretching effect deformation function ---
    [gz, dgz] = compute_iterated_stretching_deformation_function(gpt,h,shear_func,stretch_func,n_iterated);

    % --- Bending stiffness matrix ---
    h_e = [1, gpt, fz];
    for i = 1:3
        for j = 1:3
            D.Db{i,j} = D.Db{i,j} + h_e(i)*h_e(j)*Qb*gwt;
        end
    end

    % --- Shear stiffness matrix ---
    D.Ds = D.Ds + (dfz + gz)^(2)*Qs*gwt;

    % --- Bending-stretching stiffness matrix ---
    for i = 1:3
        D.Dbst{i} = D.Dbst{i} + dgz*h_e(i)*Qbst*gwt;
    end

    % --- Stretching stiffness matrix ---
    D.Dst = D.Dst + dgz^(2)*Qst*gwt;
end
end
