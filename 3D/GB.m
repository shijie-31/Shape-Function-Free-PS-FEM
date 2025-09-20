function B = GB(node3, elem3, bdNodes, gfun)
%---------------------------------------------------------------
% (c) 2025  Shijie Zhao  <shijie31123@163.com>
% PS-FEM Neumann 面荷载装配：多面体→三角面（边+面心）扇形剖分
%
% B = GB(node3, elem3, bdNodes, gfun)
% --------------------------------------------------------------
% 输入:
%   node3   : (Nnode × 3) 全局结点坐标
%   elem3   : {Nelement × 1} 元胞，每元胞存若干面，面=顶点索引行向量
%   bdNodes : 属于 Neumann 边界 Γₙ 的节点索引列向量 (可由 find 得到)
%   gfun    : Neumann 边界通量，可为
%             • 常数        -> 均匀通量
%             • @(x,y,z)    -> 返回标量  g(x,y,z)
% 输出:
%   B       : (Nnode × 1)  全局 Neumann 荷载向量
%---------------------------------------------------------------
ndof     = 1;
Nnode    = size(node3,1);
B        = zeros(Nnode,1);          % 全局右端
area_eps = 1e-14;                   % 面积容差

% 若 gfun 不是函数，把它包成常数函数
if ~isa(gfun,'function_handle')
    gconst = gfun;
    gfun   = @(x,y,z) gconst;
end

for iel = 1:numel(elem3)
    faces         = elem3{iel};                 % 该多面体所有面
    allNumbers    = [faces{:}];
    uniqueNumbers = unique(allNumbers,'stable');
    E             = numel(uniqueNumbers);       % 本单元节点数
    fe            = zeros(E,1);                 % 局部荷载

    % ------------ 遍历各面 -------------- 
    for iface = 1:numel(faces)
        nod = faces{iface}(:)';                 % 本面的顶点编号
        % —— 判断该面是否完全位于 Neumann 边界 Γₙ
        if ~all(ismember(nod, bdNodes)),  continue;  end

        numv  = numel(nod);
        eNode = node3(nod,:);                   % (numv×3)
        Cface = mean(eNode,1);                  % 面心

        % N3: 本面节点在局部编号中的位置
        N3 = zeros(1,E);
        [~,loc_face] = ismember(nod, uniqueNumbers);
        loc_face     = loc_face(loc_face>0);
        N3(loc_face) = 1/numv;                  % 本面均分

        % ----------- 沿边形成三角面 ----------- 
        for j = 1:numv
            i1 = j;
            i2 = mod(j,numv)+1;
            v1 = eNode(i1,:);
            v2 = eNode(i2,:);

            % N1/N2: 边端点指示
            j1 = nod(i1);  j2 = nod(i2);
            N1 = zeros(1,E);  N2 = N1;
            [~,idx1] = ismember(j1, uniqueNumbers);
            [~,idx2] = ismember(j2, uniqueNumbers);
            if idx1>0,  N1(idx1)=1;  end
            if idx2>0,  N2(idx2)=1;  end

            % 构造三角面  [v1, v2, Cface]
            X3   = [v1; v2; Cface];             % (3×3)
            area = 0.5 * norm(cross(X3(2,:)-X3(1,:), ...
                                    X3(3,:)-X3(1,:)));
            if area <= area_eps,  continue;  end

            % ---------- 1 点高斯积分 ----------
            Xc   = mean(X3,1);                  % 三角面重心
            gval = gfun(Xc(2),Xc(3));     % 通量值
            Nc   = (N1 + N2 + N3)/3;            % 形函数在重心

            fe = fe + (Nc.' * (gval * area));   % 累加到局部
        end
    end

    % ----------- 单元→全局装配 ------------- 
    gdof   = uniqueNumbers(:);
    B(gdof)= B(gdof) + fe;
end
end
