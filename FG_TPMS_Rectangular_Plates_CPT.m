%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title: Isogeometric analysis for Structural analyses of FG-TPMS plates with 3-Dof 2D classical plate theory (CPT) %%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% ! Please reference to paper: ............................................
% ! This work can be used, modified, and shared under the MIT License
% ! This work can be found in https://github.com/SIMOGroup/Eval-ESL-FG-TPMS-plates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================== Initialization =============================
tic
addpath(genpath('./'));
warning('off', 'all');

% clc
clear all
% close all
format long

%% ============================ Plate geometry ============================
% === Physical geometric properties ===
Plate.geo.L = 1;
Plate.geo.W = Plate.geo.L/1;
Plate.geo.h = Plate.geo.L/(10);
% Plate.geo.h = Plate.geo.W/20;

% === Define NURBS functions ===
IGA.NURBS.deg = 3; % Degree of basis functions
IGA.NURBS.ref = 11; % Number of mesh refinement

%% ============================ Material ==================================
% === Base material ===
% [1]: Steel, [2]: Alumium, [3]: Titanium, [4]: Copper, [5]: Brass
Plate.mat.type = 1;

% === Porous material ===
% [1]: Primitive, [2]: Gyroid, [3]: IWP
Plate.por_mat.type = 1;
% Plate.por_mat.RD = 0.35;

% === Transverse porosity distribution ===
% [1]: A (asymmetric), [2]: B (symmetric), [3]: C (uniform)
Plate.por_dis.type = 1;
Plate.por_dis.RD_avg = 0.35;
Plate.por_dis.RD_max = 1;
Plate.por_dis.RD_0 = 0.75;  % RD_min = RD_max*(1- RD_0);

%% ========================== Problem type ================================
% [1]: Static bending
Plate.prob.type = 1;
switch Plate.prob.type
    case 1  % Static bending
        Plate.prob.q_uniform = 1;
        IGA.result.norm_method = [1, 1, 1];  % [disp, strain, stress]
end

%% ========================= Boundary condition ===========================
% [Left, Lower, Right, Upper]
% [1]: Clamped (C), [2]: Hinged (freely movable) supported (S), [3]: Free (F), [4]: Soft (quasi-movable) supported (S2), [5]: Hard (immovable) supported (S3)
% Plate.bc.bc_case = [1,1,1,1];
Plate.bc.bc_case = [4,4,4,4];

%% ======================== Material matrices =============================
[Plate.mat_mat.D] = cal_Material_Matrices_FG_TPMS_2D_CPT(Plate);

%% =============================== IGA mesh ===============================
% === Generate NURBS mesh and properties ===    
IGA.NURBS = gen_rec_surf(Plate, IGA.NURBS);
IGA.NURBS = gen_Ien_Inn_surf(IGA.NURBS);
IGA.NURBS = gen_Iee_Ine_surf(IGA.NURBS);
IGA.NURBS = gen_FE_approx_surf(IGA.NURBS);

% Options: all, physical_shape, physical_surf, physical_element, control_net, plot_physical_element, note_physical_element, note_physical_edge, plot_control_net, note_control_net
% plot_NURBS_surf(IGA.NURBS, 'physical_surf')
% plot_NURBS_surf(IGA.NURBS, 'physical_shape', 'control_net')
% plot_NURBS_surf(IGA.NURBS, 'plot_physical_domain', 'control_net')

% === NURBS properties ===
IGA.NURBS.nsd   = 2;                                                             % Number of spatial dimension
IGA.NURBS.nnode = IGA.NURBS.mcp * IGA.NURBS.ncp;                                 % Number of control point
IGA.NURBS.nshl  = (IGA.NURBS.p + 1) * (IGA.NURBS.q + 1);                         % Number of local shape functions (= degree + 1 per element due to k refinement)
IGA.NURBS.nel   = (IGA.NURBS.mcp - IGA.NURBS.p) * (IGA.NURBS.ncp - IGA.NURBS.q); % Number of element

% === IGA properties ===
IGA.params.ndof   = 3;                                                           % Number of dofs of a control point
IGA.params.sdof   = IGA.NURBS.nnode * IGA.params.ndof;                           % Total number of dofs of the structure
IGA.params.nGauss = IGA.NURBS.p + 1;                                             % Number of gauss point in integration

%% ========================= IGA for linear geometry ======================
% === Building global matrices === 
IGA.result.K = cal_Stiffness_Matrices_2D_CPT(IGA,Plate);                    % Stiffness
switch Plate.prob.type
    case 1  % Static bending
        IGA.result.F = cal_Load_Vector_Uniform_2D_CPT(IGA,Plate);           % Extenal load
end

% === Imposing boundary conditions ===
[IGA.params.bcdof, IGA.params.bcval] = cal_bcdof_2D_CPT(IGA,Plate);
IGA.params.fdof = setdiff((1:IGA.params.sdof)', IGA.params.bcdof');  % Free dofs

%% === Solving weak form ===
switch Plate.prob.type
    case 1  % Static bending
        disp("--- Static bending responses ---")
    
        % --- Calculations ---
        bcdof = IGA.params.bcdof; bcval = IGA.params.bcval;
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        IGA.result.U = zeros(sdof, 1); 
        IGA.result.U(bcdof') = bcval';
        IGA.result.F(fdof) = IGA.result.F(fdof) - IGA.result.K(fdof, bcdof')*bcval';
        IGA.result.U(fdof) = IGA.result.K(fdof, fdof) \ IGA.result.F(fdof);

        clear bcdof bcval sdof fdof cen_disp L W plot_size_factor
end
toc

%% ========================== Post-Processing =============================
% === Deflection ===
% Options: cen_defl, NURBS_defl
cal_plot_defl_2D_CPT(IGA,Plate,'cen_defl');
% cal_plot_defl_2D_CPT(IGA,Plate,'NURBS_defl');

% === Displacement field ===
% Options: point, point_dist, FE_surf, FE_vol || ['u', 'v', 'w']
% cal_plot_disp_2D_CPT(IGA,Plate,'point',[1/2,1/4,-1/2]);
% cal_plot_disp_2D_CPT(IGA,Plate,'point_dist',[0/2,1/2],'u');
% cal_plot_disp_2D_CPT(IGA,Plate,'FE_surf','u');
% cal_plot_disp_2D_CPT(IGA,Plate,'FE_vol','u');

% === Strain field ===
% Options: point, point_dist, FE_surf, FE_vol || ['xx', 'yy', 'xy', 'yz', 'zx', 'zz']
% cal_plot_strain_2D_CPT(IGA,Plate,'point',[1/2,1/4,-1/2]);
% cal_plot_strain_2D_CPT(IGA,Plate,'point_dist',[0/2,1/2],'zx');
% cal_plot_strain_2D_CPT(IGA,Plate,'FE_surf','yy');
% cal_plot_strain_2D_CPT(IGA,Plate,'FE_vol','xx');

% === Stress field ===
% Options: point, point_dist, FE_surf, FE_vol || ['xx', 'yy', 'xy', 'yz', 'zx', 'zz']
% cal_plot_stress_2D_CPT(IGA,Plate,'point',[1/2,1/2,1/2]);
% cal_plot_stress_2D_CPT(IGA,Plate,'point',[0/2,1/2,0]);
% cal_plot_stress_2D_CPT(IGA,Plate,'point_dist',[1/2,1/2],'xx');
% cal_plot_stress_2D_CPT(IGA,Plate,'point_dist',[0/2,1/2],'zx');
% cal_plot_stress_2D_CPT(IGA,Plate,'FE_surf','yy');
% cal_plot_stress_2D_CPT(IGA,Plate,'FE_vol','xx');

