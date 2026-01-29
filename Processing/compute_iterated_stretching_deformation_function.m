function [gz, dgz] = compute_iterated_stretching_deformation_function(z,h,shear_func,stretch_func,n_iterated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute stretching effect deformation function at point z from iterated shear deformation function %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% === Shear deformation function ===
[~, dfz, ddfz] = compute_iterated_shear_deformation_function(z,h,shear_func,n_iterated);

%% === Stretching deformation function ===
% [1]: f'z, [2]: 3/20*f'z, [3]: 1/8*f'z, [4]: 1/12*f'z, [5]: 1/20*f'z
switch stretch_func
    case 1
        gz = dfz;
        dgz = ddfz;
    case 2
        gz = 3/20*dfz;
        dgz = 3/20*ddfz;
    case 3
        gz = 1/8*dfz;
        dgz = 1/8*ddfz;
    case 4
        gz = 1/12*dfz;
        dgz = 1/12*ddfz;
    case 5
        gz = 1/20*dfz;
        dgz = 1/20*ddfz;
    otherwise
        disp('Error')
        pause
end

end
