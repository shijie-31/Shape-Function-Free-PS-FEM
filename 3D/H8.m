function poisson3D_hex(Nx)
% poisson3D_hex  Poisson solver on (0,1)^3 using trilinear HEX8 elements
%
%   poisson3D_hex         % default Nx = 32 subdivisions
%   poisson3D_hex(Nx)     % user-defined subdivisions per edge
%
%   Outputs L2 and H1 errors against analytical solution:
%           u(x,y,z) = sin(2*x*y)*cos(z)

% -------------------------------------------------------------------------
if nargin==0, Nx =64; end

% --- generate structured hex mesh ----------------------------------------
[coords, elems] = generate_hex_mesh(Nx);     % nodes & 8-node connectivity

% --- exact solution + gradient + RHS -------------------------------------
uex   = @(x,y,z) sin(2*x.*y).*cos(z);
gradU = @(x,y,z)[ 2*y.*cos(z).*cos(2*x.*y) ; ...
                  2*x.*cos(z).*cos(2*x.*y) ; ...
                 -sin(2*x.*y).*sin(z)       ];         % 3×N
f_rhs = @(x,y,z)(4*y.^2 + 4*x.^2 + 1).*sin(2*x.*y).*cos(z);

% --- global assembly ------------------------------------------------------
nn = size(coords,1);      ne = size(elems,1);
K  = sparse(nn,nn);
F  = zeros(nn,1);

% 2×2×2 Gauss rule on reference cube [-1,1]^3 ------------------------------
g  = 1/sqrt(3);                 % Gauss abscissa
xi = [-g  g];
[Xi,Eta,Zta] = ndgrid(xi,xi,xi);
gp   = [Xi(:), Eta(:), Zta(:)]; % 8×3
wgt  = ones(8,1);               % each weight = 1

for e = 1:ne
    vid = elems(e,:);                   % 1×8 node ids
    v   = coords(vid,:);                % 8×3 coordinates
    
    Ke = zeros(8,8);  Fe = zeros(8,1);
    for q = 1:8
        xi  = gp(q,1);  eta = gp(q,2);  zta = gp(q,3);
        [N, dNdr] = shape_hex8(xi,eta,zta);   % N:1×8 , dNdr:3×8
        
        J   = v.' * dNdr.';              % 3×3 Jacobian
        detJ= det(J);
        B   = (J \ dNdr).';              % 8×3 (each row grad φ_i)
        
        Ke = Ke + (B*B.') * detJ * wgt(q);
        
        % physical coordinates of GP
        xq = N * v;                      % 1×3
        fq = f_rhs(xq(1),xq(2),xq(3));
        Fe = Fe + N.' * fq * detJ * wgt(q);
    end
    K(vid,vid) = K(vid,vid) + Ke;
    F(vid)     = F(vid)     + Fe;
end

% --- Dirichlet boundary ---------------------------------------------------
bnd = find( any(coords==0 | coords==1 , 2) );
uD  = uex(coords(bnd,1), coords(bnd,2), coords(bnd,3));

free        = setdiff(1:nn, bnd);
F           = F - K(:,bnd)*uD;
uh          = zeros(nn,1);
uh(bnd)     = uD;
uh(free)    = K(free,free) \ F(free);

% --- error computation (reuse 2×2×2 rule) ---------------------------------
[L2err, H1err] = compute_errors_hex(coords, elems, uh, uex, gradU, gp, wgt);

fprintf('Nx = %d   DOF = %d   Hex = %d\n', Nx, nn, ne);
fprintf('‖u  – uh‖_L2  = %.3e\n', L2err);
fprintf('‖∇u – ∇uh‖_L2 = %.3e  (H1-seminorm)\n', H1err);
end
% ==========================================================================
function [coords, elems] = generate_hex_mesh(N)
% structured mesh on (0,1)^3 : (N+1)^3 nodes, N^3 bricks

[x,y,z]      = ndgrid(linspace(0,1,N+1));
coords       = [x(:), y(:), z(:)];                 % (N+1)^3 × 3
nodeID       = reshape(1:size(coords,1), N+1, N+1, N+1);
elems        = zeros(N^3, 8);   cnt = 0;

for k = 1:N
  for j = 1:N
    for i = 1:N
      cnt = cnt + 1;
      n000 = nodeID(i  ,j  ,k  );
      n100 = nodeID(i+1,j  ,k  );
      n110 = nodeID(i+1,j+1,k  );
      n010 = nodeID(i  ,j+1,k  );
      n001 = nodeID(i  ,j  ,k+1);
      n101 = nodeID(i+1,j  ,k+1);
      n111 = nodeID(i+1,j+1,k+1);
      n011 = nodeID(i  ,j+1,k+1);
      elems(cnt,:) = [n000 n100 n110 n010 n001 n101 n111 n011];
    end
  end
end
end
% ==========================================================================
function [N, dNdr] = shape_hex8(r,s,t)
% shape functions & derivatives at (r,s,t) in [-1,1]^3
N = 0.125*[(1-r)*(1-s)*(1-t), (1+r)*(1-s)*(1-t), ...
           (1+r)*(1+s)*(1-t), (1-r)*(1+s)*(1-t), ...
           (1-r)*(1-s)*(1+t), (1+r)*(1-s)*(1+t), ...
           (1+r)*(1+s)*(1+t), (1-r)*(1+s)*(1+t)];              % 1×8

dNdr = 0.125*[
 -(1-s)*(1-t)  ,  (1-s)*(1-t) ,  (1+s)*(1-t) , -(1+s)*(1-t) , ...
 -(1-s)*(1+t)  ,  (1-s)*(1+t) ,  (1+s)*(1+t) , -(1+s)*(1+t) ;  % ∂/∂r
 -(1-r)*(1-t)  , -(1+r)*(1-t) ,  (1+r)*(1-t) ,  (1-r)*(1-t) , ...
 -(1-r)*(1+t)  , -(1+r)*(1+t) ,  (1+r)*(1+t) ,  (1-r)*(1+t) ;  % ∂/∂s
 -(1-r)*(1-s)  , -(1+r)*(1-s) , -(1+r)*(1+s) , -(1-r)*(1+s) , ...
  (1-r)*(1-s)  ,  (1+r)*(1-s) ,  (1+r)*(1+s) ,  (1-r)*(1+s) ]; % ∂/∂t
end
% ==========================================================================
function [L2, H1] = compute_errors_hex(X, T, uh, uex, gradU, gp, wgt)
L2 = 0;  H1 = 0;
for e = 1:size(T,1)
    vid = T(e,:);    v = X(vid,:);      uh_e = uh(vid);

    for q = 1:8
        xi  = gp(q,1);  eta = gp(q,2);  zta = gp(q,3);
        [N, dNdr] = shape_hex8(xi,eta,zta);
        J   = v.' * dNdr.';                   detJ = det(J);
        B   = (J \ dNdr).';                   % 8×3

        xq  = N * v;                          % GP phys coord
        uq  = uex(xq(1),xq(2),xq(3));
        uhq = N * uh_e;

        gu  = gradU(xq(1),xq(2),xq(3));
        guh = B.' * uh_e;                     % 3×1

        L2 = L2 + wgt(q) * detJ * (uq - uhq)^2;
        H1 = H1 + wgt(q) * detJ * sum((gu - guh).^2);
    end
end
L2 = sqrt(L2);   H1 = sqrt(H1);
end
