% ---------------------------------------------------------------
% 3-D Poisson  Δu + Q = 0  on Ω=(0,1)^3
%   Patch test: u_exact = 1 + x + 2y - 3z   (linear)
%   => Δu = 0,  Q = 0
% ---------------------------------------------------------------
%format long g 
% load your mesh
 %load('C:\Users\shiji\Documents\MATLAB\POISSON\mesh3data3.mat')
 % 应包含 node3, elem3
%load('C:\Users\shiji\Documents\MATLAB\POISSON\p1.mat')
option.FaceAlpha = 1;
showmesh(node3, elem3, option);

% ----------------------- 刚度矩阵 & 体力项 -----------------------
K  = GK(node3, elem3);

% 贴片测试：Q = 0
Qfun = @(x,y,z) 0.*x;           % 或者直接 0 也行，但保持函数接口
F  = GF(node3, elem3, Qfun);

% ------------------------ Neumann 面荷载 ------------------------
% 仅在 x=0 面施加一致的 Neumann：∂u/∂n = ∇u·n
% u = 1 + x + 2y - 3z  => ∇u = (1, 2, -3)
% x=0 面外法向 n = (-1,0,0)  => ∂u/∂n = (1,2,-3)·(-1,0,0) = -1
tol = 1e-9;
x0_nodes = find(abs(node3(:,1)) < tol);
gfun = @(y,z) -1 + 0.*y;        % 常数 -1，写成匿名函数以适配 GB 签名
F = F + GB(node3, elem3, x0_nodes, gfun);

% ---------------------- Dirichlet 边界处理 ----------------------
u_exact = @(p) 1 + p(:,1) + 2*p(:,2) - 3*p(:,3);

x1_nodes = find(abs(node3(:,1) - 1) < tol);
y0_nodes = find(abs(node3(:,2) - 0) < tol);
y1_nodes = find(abs(node3(:,2) - 1) < tol);
z0_nodes = find(abs(node3(:,3) - 0) < tol);
z1_nodes = find(abs(node3(:,3) - 1) < tol);

% 除了 x=0（已设 Neumann）外，其余五个面都设 Dirichlet
bdNode = unique([x1_nodes; y0_nodes; y1_nodes; z0_nodes; z1_nodes]);
uD     = u_exact(node3(bdNode,:));

allNode  = (1:size(node3,1))';
freeNode = setdiff(allNode, bdNode);

K_ff = K(freeNode, freeNode);
F_f  = F(freeNode) - K(freeNode, bdNode)*uD;

% -------------------------- 求解 -------------------------------
uh            = zeros(size(node3,1),1);
uh(freeNode)  = K_ff \ F_f;
uh(bdNode)    = uD;

% -------------------------- 误差 -------------------------------
u_val  = u_exact(node3);
e      = uh - u_val;
m_vec  = GF(node3, elem3, 1);           % lumped 质量向量

L2_err  = sqrt( sum( (e.^2) .* m_vec ) );
L2_norm = sqrt( sum( (u_val.^2) .* m_vec ) );
rel_L2  = L2_err / L2_norm;

fprintf('Patch test | Absolute L2 error : %.3e\n', L2_err);
fprintf('Patch test | Relative L2 error : %.3e\n', rel_L2);

% -------------------------- 显示 -------------------------------
showsolution3D(node3, elem3, uh);
title('Patch Test: u_h for linear exact solution');
