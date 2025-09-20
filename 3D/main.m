% ---------------------------------------------------------------
% 3-D Poisson  Δu + Q = 0  on Ω=(0,1)^3
%   u_exact = sin(2*x*y)*cos(z)
% ---------------------------------------------------------------
load('Q1.mat')
   option.FaceAlpha = 1;
     showmesh(node3,elem3,option); 

% ----------------------- 装配矩阵和体力项 -----------------------
K  = GK(node3, elem3);

Qfun = @(x, y, z) (4*(x.^2 + y.^2) + 1) .* sin(2*x.*y) .* cos(z);
F  = GF(node3, elem3, Qfun);

% ------------------------ Neumann 面荷载 ------------------------
tol = 1e-9;
x0_nodes = find(abs(node3(:,1)) < tol);               % x = 0  面

% ∂u/∂n = -u_x  (外法向 (-1,0,0) )
gfun = @(y,z) -2*y.*cos(z);   % 只接受 2 个参数                  % g(y,z) = -2y cos z

% GB  = 你的“面荷载”装配函数 ：把  g φ_i  积分到 F 中
F = F + GB(node3, elem3, x0_nodes, gfun);           

% ---------------------- Dirichlet 处理 --------------------------
u_exact = @(p) sin(2*p(:,1).*p(:,2)).*cos(p(:,3));

x1_nodes = find(abs(node3(:,1) - 1) < tol);
y0_nodes = find(abs(node3(:,2) - 0) < tol);
y1_nodes = find(abs(node3(:,2) - 1) < tol);
z0_nodes = find(abs(node3(:,3) - 0) < tol);
z1_nodes = find(abs(node3(:,3) - 1) < tol);

bdNode   = unique([x1_nodes; y0_nodes; y1_nodes; z0_nodes; z1_nodes]);
uD       = u_exact(node3(bdNode,:));

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
m_vec  = GF(node3, elem3, 1);     % lumped 质量向量：m_i = ∫ φ_i

L2_err  = sqrt( sum( (e.^2) .* m_vec ) );          % ‖e‖_L2
L2_norm = sqrt( sum( (u_val.^2) .* m_vec ) );      % ‖u‖_L2
rel_L2  = L2_err / L2_norm;

fprintf('Absolute  L2 error : %.3e\n', L2_err);
fprintf('Relative  L2 error : %.3e\n', rel_L2);

% -------------------------- 显示 -------------------------------
showsolution3D(node3, elem3, uh);
