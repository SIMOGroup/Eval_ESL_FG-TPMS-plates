function [NURBS_defl] = cal_deflection_2D_RPT(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate deflection (transverse displacement) of all control points of rectangular plate with 4-Dof 2D RPT plate models %%%
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

%% ===== Initial deformation NURBS =====
NURBS_defl = zeros(size(CP,1),size(CP,2));

%% ====== Plot deflection ======
for j = 1:size(CP,2)
    for i = 1:size(CP,1) 
        idx = (j-1)*size(CP,1) + i;
        
        disp = U(ndof*(idx-1)+3) + U(ndof*(idx-1)+4);
        
        NURBS_defl(i,j) = disp;
    end
end
end
