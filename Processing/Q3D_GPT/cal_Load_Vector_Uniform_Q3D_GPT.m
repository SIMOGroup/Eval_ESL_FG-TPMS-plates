function F = cal_Load_Vector_Uniform_Q3D_GPT(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate load vector with 6-Dof Quasi-3D GPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Inn = IGA.NURBS.Inn; Ien = IGA.NURBS.Ien; nel = IGA.NURBS.nel;
sdof = IGA.params.sdof; ndof = IGA.params.ndof; nGauss = IGA.params.nGauss; 

%% Used parameters from Plate
q_uniform = Plate.prob.q_uniform;
h = Plate.geo.h;
shear_func = Plate.theory.shear_func; n_iterated = Plate.theory.iterated_order; stretch_func = Plate.theory.stretch_func;

%% ===== Initial load vector =====
I_idx = zeros(ndof*numel(Ien),1); V_idx = I_idx;

%% ===== Initial Gauss integration =====
[Pg, Wg] = gen_Gauss_point(1,-1,nGauss);
nel_nza = 0; % Number of non-zero-area elements
tol = 1e-8;

%% ===== Stiffness matrices =====
for iel = 1:nel  % Loop over the elements
    % --- Element parameters ---
    sctr = Ien(iel,:);  % Control points indexes
    nn = length(sctr);  % Number of control points in the element
    nedof = ndof*nn;    % Number of dofs in the element
    for idof = 1:ndof  % Dofs of control points
        sctrF(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
    end
    ni = Inn(Ien(iel,1),1);  % Index of the element in parametric domain
    nj = Inn(Ien(iel,1),2);
    
    Fe = zeros(nedof, 1);  % Element load vector

    % --- Gauss integration ---
    if abs(uKnot(ni)-uKnot(ni+1)) > tol && abs(vKnot(nj)-vKnot(nj+1)) > tol  % Check if the current element has nonzero area in the parametric domain
        nel_nza = nel_nza + 1;
        detJ2_xi = (uKnot(ni+1) - uKnot(ni))/2;  % Mapping parametric domain into natural domain of [[-1, 1]; [-1, 1]]
        detJ2_eta = (vKnot(nj+1) - vKnot(nj))/2;

        for iGauss = 1: nGauss  % Loop over the integration points
            for jGauss = 1: nGauss
                % Gauss points & weights
                gpt_xi = Pg(iGauss); gwt_xi = Wg(iGauss);
                gpt_eta = Pg(jGauss); gwt_eta = Wg(jGauss);
                
                % Map the point to parametric domain
                gpt_xi = (uKnot(ni+1)-uKnot(ni))/2*gpt_xi + (uKnot(ni+1)+uKnot(ni))/2; gwt_xi = gwt_xi*detJ2_xi;
                gpt_eta = (vKnot(nj+1)-vKnot(nj))/2*gpt_eta + (vKnot(nj+1)+vKnot(nj))/2; gwt_eta = gwt_eta*detJ2_eta;
                gwt = gwt_xi*gwt_eta;
                
                % Kinematic matrices of NURBS shape function
                [N, ~, ~, ~, ~, detJ1] = cal_Kine_Shape_2D_2nd(IGA.NURBS,ni,nj,gpt_xi,gpt_eta);
                
                % Stretching effect deformation function
                [gz, ~] = compute_iterated_stretching_deformation_function(0,h,shear_func,stretch_func,n_iterated);
                
                % Force vector
                N0 = zeros(1,ndof*nn);
                N0(1,3:ndof:ndof*nn) = N';
                N0(1,6:ndof:ndof*nn) = gz*N';
                
                Fe = Fe + N0'*q_uniform*gwt*detJ1;
            end
        end

        % --- Assemble to global matrix ---
        idx = (iel-1)*(nedof) + (1:nedof);
        I_idx(idx) = sctrF; V_idx(idx) = Fe(:);

    end
end

% --- Assemble to global matrix ---
F = sparse(I_idx, ones(size(I_idx)), V_idx, sdof, 1);

end
