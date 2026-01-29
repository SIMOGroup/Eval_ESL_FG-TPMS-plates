function [stress_norm] = compute_stress_normalization(stress,norm_method,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute normalization for stress vectors %%%
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
        stress_norm = stress;
    case 1
        norm = h/(L*q_uniform);
        stress_norm = stress*norm;
    case 2
        norm = (h^2)/(L^2*q_uniform);
        stress_norm = stress*norm;
    otherwise
        disp('Error')
        pause
end

end
