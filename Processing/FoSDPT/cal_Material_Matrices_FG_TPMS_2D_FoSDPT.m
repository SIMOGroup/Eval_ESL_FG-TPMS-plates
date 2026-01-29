function [D] = cal_Material_Matrices_FG_TPMS_2D_FoSDPT(Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate material matrices for FG-TPMS for 7-Dof 2D SS's FoSDPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
h = Plate.geo.h;

%% ===== Initial material matrices =====
for i = 4:-1:1
    for j = 4:-1:1
        D.Db{i,j} = zeros(3,3);
    end
end
for i = 2:-1:1
    for j = 2:-1:1
        D.Ds{i,j} = zeros(2,2);
    end
end

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
    h_e = [1, 4/3*(gpt^3)/h^2, gpt - 4/3*(gpt^3)/h^2, gpt^2 + 2*(gpt^4)/h^2];
    for i = 1:4
        for j = 1:4
            D.Db{i,j} = D.Db{i,j} + h_e(i)*h_e(j)*Qb*gwt;
        end
    end

    % --- Shear stiffness matrix ---
    h_s = [(1 - 4*(gpt^2)/h^2), 2*gpt*(1 - 4*(gpt^2)/h^2)];
    for i = 1:2
        for j = 1:2
            D.Ds{i,j} = D.Ds{i,j} + h_s(i)*h_s(j)*Qs*gwt;
        end
    end
end
end
