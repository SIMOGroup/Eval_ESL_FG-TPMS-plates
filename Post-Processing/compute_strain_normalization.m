function [strain_norm] = compute_strain_normalization(strain,norm_method,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute normalization for strain vectors %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
[E_m, nu_m, ~, ~] = compute_basic_material(Plate.mat.type);
D_m = E_m*h^(3)/(12*(1-nu_m^2));

if exist('Plate.prob.q_uniform', 'var')
    q_uniform = Plate.prob.q_uniform;
else
    q_uniform = 1;
end

%% === Normalization ===
switch norm_method
    case 0
        strain_norm = strain;
    case 1
        norm = E_m/q_uniform * h^2/L^2;
        strain_norm = strain*norm;
    otherwise
        disp('Error')
        pause
end

end
