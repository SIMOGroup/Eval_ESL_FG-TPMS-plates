function [NURBS_defl] = cal_deflection_Q3D_RPT(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate deflection (transverse displacement) of all control points of rectangular plate with 4-Dof Quasi-3D RPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
CP = IGA.NURBS.CP;
ndof = IGA.params.ndof; 
U = IGA.result.U;

%% Used parameters from Plate
h = Plate.geo.h;
shear_func = Plate.theory.shear_func; n_iterated = Plate.theory.iterated_order; 
stretch_func = Plate.theory.stretch_func;

%% ===== Initial deformation NURBS =====
NURBS_defl = zeros(size(CP,1),size(CP,2));

%% ====== Plot deflection ======
% Stretching effect deformation function
[gz, ~] = compute_iterated_stretching_deformation_function(0,h,shear_func,stretch_func,n_iterated);

for j = 1:size(CP,2)
    for i = 1:size(CP,1) 
        idx = (j-1)*size(CP,1) + i;
        
        disp = U(ndof*(idx-1)+3) + gz*U(ndof*(idx-1)+4);
        
        NURBS_defl(i,j) = disp;
    end
end
end
