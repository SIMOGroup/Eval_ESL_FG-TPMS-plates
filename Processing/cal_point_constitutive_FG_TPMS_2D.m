function [Q] = cal_point_constitutive_FG_TPMS_2D(Plate, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate consitutive model of a specific material point for FG-TPMS for 2D plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
h = Plate.geo.h;
mat_type = Plate.mat.type;
porous_type = Plate.por_mat.type; % RD = Plate.por_mat.RD;
por_dis = Plate.por_dis.type; RD_avg = Plate.por_dis.RD_avg; RD_max = Plate.por_dis.RD_max; RD_0 = Plate.por_dis.RD_0;

%% ===== Initial material matrices =====
Q.Qb = zeros(3,3);
Q.Qs = zeros(2,2);

%% ===== Define material properties =====
% --- TPMS parameters ---
switch porous_type
    case 1 % Primitive
        k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
        k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
        k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
    case 2 % Gyroid
        k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
        k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
        k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
    case 3 % IWP
        k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
        k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
        k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
end
C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
d_1 = 0.3 - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = 0.3 - a_2*(1)^2 - b_2(1);

%% Calculate n_{A}, n_{B}, psi_{C} of porous ditribution 1, 2, 3
switch por_dis
    case 1
        n_A = (RD_max - RD_avg) / (RD_avg - RD_max*(1-RD_0));
    case 2
        try
            fun = @(n) integral(@(u) 2/pi * RD_max * (1-RD_0 + RD_0*(1 - u).^n) ./ sqrt(1 - u.^2), 0, 1) - RD_avg;
            n_B = fzero(fun, 1);
        catch
            disp('Error in finding n_B');
            pause
        end
    case 3
        psi = (RD_avg - RD_max*(1-RD_0)) / (RD_max*RD_0);
end

%% ===== Calculate material matrices =====
% --- Material properties ---
[E_s, nu_s, rho_s, ~] = compute_basic_material(mat_type);
G_s = E_s / (2*(1+nu_s));

% --- Porosity distribution ---
switch por_dis
    case 1
        psi_z = (z/h + 1/2)^n_A;
    case 2
        psi_z = (1 - cos(pi*z/h))^n_B;
    case 3
        psi_z = psi;
end
RD = RD_max * (1-RD_0 + RD_0*psi_z);

% --- Porous material properties ---
if RD == 1
    E = E_s; nu = nu_s;
    G = E /(2*(1+nu));
else
    e = (RD <= k_e) * (C_1e*RD^(n_1e)) + ...
        (RD >  k_e) * (C_2e*RD.^(n_2e) + C_3e);
    g = (RD <= k_g) * (C_1g*RD^(n_1g)) + ...
        (RD >  k_g) * (C_2g*RD.^(n_2g) + C_3g);
    nu =(RD <= k_nu) * (a_1.*exp(b_1*RD) + d_1) + ...
        (RD >  k_nu) * (a_2*RD^(2) + b_2*RD + d_2);
    E = E_s*e;
    G = G_s*g;
end
rho = rho_s * RD;
        
% --- Consitutive matrix ---
Q_11 = E / (1-nu^(2)); Q_12 = (E*nu) / (1-nu^(2)); C_44 = G;
Q.Qb = [Q_11, Q_12,    0; ...
        Q_12, Q_11,    0; ...
           0,    0, C_44];
Q.Qs = [C_44,    0; ...
           0, C_44];

end
