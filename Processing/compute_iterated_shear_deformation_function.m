function [fz, dfz, ddfz] = compute_iterated_shear_deformation_function(z, h, shear_func, n_iterated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute iterated shear deformation function at point z/h %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% ["Poly"]: Polynomial, ["Trig"]: Trigonometric, ["Hype"]: Hyperbolic, ["Expo"]: Exponential, ["Comb"]: Combined or Hybrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
n_iterated = int32(n_iterated);  % Integer number
gz = NaN(1, n_iterated+1); dgz = gz; ddgz = gz;

%% === Iterated shear deformation function ===
gz(1) = z; dgz(1) = 1; ddgz(1) = 0;
for i_iter = 2:n_iterated+1
    [gz(i_iter), dgz(i_iter), ddgz(i_iter)] = compute_shear_deformation_function(gz(i_iter-1), h, shear_func);
end

fz = gz(end);
dfz = prod(dgz);
ddfz = sum(ddgz(2:end) .* [1, cumprod(dgz(2:end-1).^2)] .*  flip([1, cumprod(flip(dgz(3:end)))]));

% ddfz = ddgz(1);
% for k = 2:n_iterated+1
%     mul_1(k) = prod(dgz(2:k-1).^2);
%     mul_2(k) = prod(dgz(k+1:end));
%     ddfz = ddfz + ddgz(k) * mul_1(k) * mul_2(k);
% end

end
