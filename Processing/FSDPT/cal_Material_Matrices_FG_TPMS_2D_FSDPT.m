function [D] = cal_Material_Matrices_FG_TPMS_2D_FSDPT(Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate material matrices for FG-TPMS for 5-Dof 2D RM's FSDPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
h = Plate.geo.h;
k_scf = Plate.theory.k_scf;

%% ===== Initial material matrices =====
for i = 2:-1:1
    for j = 2:-1:1
        D.Db{i,j} = zeros(3,3);
    end
end
D.Ds = zeros(2,2);

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
    [Q] = cal_point_constitutive_FG_TPMS_2D(Plate, gpt);
    Qb = Q.Qb; Qs = Q.Qs;
    
    % --- Bending stiffness matrix ---
    h_e = [1, gpt];
    for i = 1:2
        for j = 1:2
            D.Db{i,j} = D.Db{i,j} + h_e(i)*h_e(j)*Qb*gwt;
        end
    end

    % --- Shear stiffness matrix ---
    D.Ds = D.Ds + k_scf*Qs*gwt;
end
end
