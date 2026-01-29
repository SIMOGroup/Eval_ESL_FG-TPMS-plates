function [disp_norm] = compute_displacement_normalization(displ,norm_method,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute normalization for displacement vectors %%%
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
        disp_norm = displ;
    case 1
        norm = 100*D_m/(q_uniform*L^4);
        disp_norm = displ*norm;
    case 2
        norm = 1/h;
        disp_norm = displ*norm;
    case 3
        norm = E_m*h^2 / (q_uniform*L^3);
        disp_norm = displ*norm;
    otherwise
        disp('Error')
        pause
end

end
