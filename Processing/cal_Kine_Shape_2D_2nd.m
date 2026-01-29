function [R, dRdxi, d2Rdxi2, dRdx, d2Rdx2, detJ] = cal_Kine_Shape_2D_2nd(NURBS,ni,nj,u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate 2D kinetic shape functions from 2nd gradients of all non-zero NURBS basis functions %%%
%%% with respect to parameter variables (xi, eta) and physical variables (x, y) %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA.NURBS
p = NURBS.p; mcp = NURBS.mcp; uKnot = NURBS.uKnot; 
q = NURBS.q; ncp = NURBS.ncp; vKnot = NURBS.vKnot;
B_net = NURBS.B_net; 
nsd = NURBS.nsd; nshl = NURBS.nshl;

%% ===== 1D basis functions and their derivatives =====
M = eval_ders_basis_func(ni,p,u,uKnot);
N = eval_ders_basis_func(nj,q,v,vKnot);

%% ===== NURBS functions and their derivatives w.r.t. parametric domain (R, dR/dxi, d2R/dxi2) =====
% === Numerator, A = N*M*w ===
A = zeros(nshl,1); dAdxi = zeros(nshl,2); d2Adxi2 = zeros(nshl,3);
for j = 0:q
    for i = 0:p
        i_shl = j*(p+1) + i + 1;

        % --- Basis functions --- 
        A(i_shl,1) = M(1,p+1-i) * N(1,q+1-j) * B_net(ni-i,nj-j,nsd+1);
        
        % --- First derivatives, [dA/dxi, dA/deta] = [dN/dxi*M*w, N*dM/deta*w] ---
        dAdxi(i_shl,1) = M(2,p+1-i) * N(1,q+1-j) * B_net(ni-i,nj-j,nsd+1);  % dA/dxi
        dAdxi(i_shl,2) = M(1,p+1-i) * N(2,q+1-j) * B_net(ni-i,nj-j,nsd+1);  % dA/deta

        % --- Second derivatives, [d2A/dxi2, d2A/deta2, d2A/dxideta] = [d2N/dxi2*M*w, N*d2M/deta2*w, dN/dxi*dM/deta*w] ---
        d2Adxi2(i_shl,1) = M(3,p+1-i) * N(1,q+1-j) * B_net(ni-i,nj-j,nsd+1);  % d2A/dxi2
        d2Adxi2(i_shl,2) = M(1,p+1-i) * N(3,q+1-j) * B_net(ni-i,nj-j,nsd+1);  % d2A/deta2
        d2Adxi2(i_shl,3) = M(2,p+1-i) * N(2,q+1-j) * B_net(ni-i,nj-j,nsd+1);  % d2A/dxideta
    end
end

% === Denominator, B = sum(N*M*w) ===
% --- Basis functions ---
B = sum(A);

% --- First derivatives, [dB/dxi, dB/deta] = [sum(dN/dxi*M*w), sum(N*dM/deta*w)] ---
dBdxi = sum(dAdxi);

% --- Second derivatives, [d2B/dxi2, d2B/deta2, d2B/dxideta] = [sum(d2N/dxi2*M*w), sum(N*d2M/deta2*w), sum(dN/dxi*dM/deta*w)] ---
d2Bdxi2 = sum(d2Adxi2);

% === Rational functions, R = N*M*w / sum(N*M*w)===
% --- Basis functions, R = A/B --- 
R = A/B;

% --- First derivatives, [dR/dxi, dR/deta] ---
dRdxi = zeros(nshl,nsd);
dRdxi(:,1) = dAdxi(:,1)/B - dBdxi(:,1)*R/B;
dRdxi(:,2) = dAdxi(:,2)/B - dBdxi(:,2)*R/B;

% --- Second derivatives, [d2R/dxi2, d2R/deta2, d2R/dxideta] ---
d2Rdxi2 = zeros(nshl,nsd+1);
d2Rdxi2(:,1) = d2Adxi2(:,1)/B - 2*dAdxi(:,1)*dBdxi(:,1)/(B^2) - d2Bdxi2(:,1)*R/B + 2*(dBdxi(:,1)^2)*R/(B^2);
d2Rdxi2(:,2) = d2Adxi2(:,2)/B - 2*dAdxi(:,2)*dBdxi(:,2)/(B^2) - d2Bdxi2(:,2)*R/B + 2*(dBdxi(:,2)^2)*R/(B^2);
d2Rdxi2(:,3) = d2Adxi2(:,3)/B - (dAdxi(:,1)*dBdxi(:,2) + dAdxi(:,2)*dBdxi(:,1))/(B^2) - d2Bdxi2(:,3)*R/B + 2*(dBdxi(:,1)*dBdxi(:,2))*R/(B^2);

%% ===== Derivative of physical domain w.r.t. parametric domain (dxdxi, d2xdxi2) =====
dxdxi = zeros(nsd,nsd); d2xdxi2 = zeros(nsd,nsd+1);
for j = 0:q
    for i = 0:p
        i_shl = j*(p+1) + i + 1;

        % --- First derivatives, [dx/dxi, dx/deta; dy/dxi, dy/deta] ---
        dxdxi = dxdxi + (dRdxi(i_shl,:)' * reshape(B_net(ni-i,nj-j,1:2),1,[]))';

        % --- Second derivatives, [d2x/dxi2, d2x/deta2, d2x/dxideta; d2y/dxi2, d2y/deta2, d2y/dxideta] ---
        d2xdxi2 = d2xdxi2 + (d2Rdxi2(i_shl,:)' * reshape(B_net(ni-i,nj-j,1:2),1,[]))';
    end
end

%% ===== Derivative of parametric domain w.r.t. physical domain (dxidx, d2xidx2) =====
% --- First derivatives, [dxi/dx, dxi/dy; deta/dx, deta/dy] ---
dxidx = dxdxi \ eye(2); detJ = det(dxdxi);

% --- Second derivatives, [d2xi/dx2, d2xi/dy2, d2xi/dxdy; d2eta/dx2, d2eta/dy2, d2eta/dxdy] ---
T = [         dxidx(1,1)^2,          dxidx(2,1)^2,                       2*dxidx(1,1)*dxidx(2,1);
              dxidx(1,2)^2,          dxidx(2,2)^2,                       2*dxidx(1,2)*dxidx(2,2);
     dxidx(1,1)*dxidx(1,2), dxidx(2,1)*dxidx(2,2), dxidx(1,1)*dxidx(2,2) + dxidx(1,2)*dxidx(2,1)];  % Quadratic transformation matrix
d2xidx2 = - dxidx*d2xdxi2*T';

%% ===== Derivative of NURBS functions w.r.t. physical domain (R, dR/dx, d2R/dx2) =====
% --- First derivatives, [dR/dx, dR/dy] ---
dRdx = dRdxi * dxidx;

% --- Second derivatives, [d2R/dx2, d2R/dy2, d2R/dxdy] ---
d2Rdx2 = d2Rdxi2*T' + dRdxi*d2xidx2;

end
