function K = GK(node3, elem3)
%   --------------------------------------------------------------------
%   (c) 2025 Shijie Zhao <shijie31123@163.com>
%   多面体光滑有限元（PS-FEM）装配：多面体→四面体（边+面心+单元心）扇形剖分
%   本版采用简化形函数：
%       N1: 边起点顶点指示函数
%       N2: 边终点顶点指示函数
%       N3: 面上结点“均分” (每个面结点 1/num)
%       N4: 单元内所有结点“均分” (每个结点 1/E)
%   --------------------------------------------------------------------

ndof  = 1;                              % 每节点自由度（标量场）
Nnode = size(node3, 1);
K     = sparse(Nnode, Nnode);           % 全局刚度
vol_eps = 1e-14;                        % 体积容差，避免数值退化

for iel = 1:numel(elem3)
    faces = elem3{iel};                 % 该多面体的所有面(元胞数组)

    % 该单元涉及的结点集合（保持首次出现顺序）
    allNumbers    = [faces{:}];
    uniqueNumbers = unique(allNumbers, 'stable');   % 行/列均可，下面统一(:)
    E  = numel(uniqueNumbers);
    C  = node3(uniqueNumbers, :);       % (E×3) 单元使用到的结点坐标
    CV = mean(C, 1);                    % 单元中心
    Ke = sparse(E, E);                  % 单元局部刚度

    % N4：单元内“全域均分”
    N4 = (1/E) * ones(1, E);

    % 遍历每个面
    for iface = 1:numel(faces)
        nod   = faces{iface}(:)';       % 该面的结点编号（行向量）
        numv  = numel(nod);
        eNode = node3(nod, :);          % (numv×3) 该面各顶点坐标
        Cface = mean(eNode, 1);         % 面心

        % N3：该面“均分”，仅在本面结点上为 1/numv
        N3 = zeros(1, E);
        [~, loc_face] = ismember(nod, uniqueNumbers);
        loc_face = loc_face(loc_face > 0);    % 过滤掉 0（未找到）
        if ~isempty(loc_face)
            N3(loc_face) = 1/numv;
        end

        % 沿面边逐一构造四面体：[e1, e2, Cface, CV]
        for j = 1:numv
            i1 = j;
            i2 = mod(j, numv) + 1;
            e1 = eNode(i1, :);
            e2 = eNode(i2, :);

            % N1/N2：边两端的“顶点指示”
            j1 = nod(i1);
            j2 = nod(i2);

            N1 = zeros(1, E);
            [~, idx1] = ismember(j1, uniqueNumbers);
            if idx1 > 0, N1(idx1) = 1; end

            N2 = zeros(1, E);
            [~, idx2] = ismember(j2, uniqueNumbers);
            if idx2 > 0, N2(idx2) = 1; end

            % 构造 T4 并获取几何量
            X4 = [e1; e2; Cface; CV];   % (4×3)
            [VT4, ST4, nMT4] = G_T4(X4);

            Veff = abs(VT4);
            if Veff <= vol_eps
                continue;               % 跳过退化四面体
            end

            ST4 = ST4(:);               % 期望为 4×1
            if numel(ST4) ~= 4 || ~isequal(size(nMT4), [4, 3])
                error('G_T4 输出尺寸不符合预期：ST4 应为 4×1，nMT4 应为 4×3。');
            end

            % ------------------ 当前四面体的 B 矩阵 (3×E) ------------------
            % 四个面分别为：(1,2,3), (1,2,4), (1,3,4), (2,3,4)
            % 面平均权重 1/3，乘以对应面面积并除以体积
            B = zeros(3, E);
            face_sum_1 = N1 + N2 + N3;  % (1,2,3)
            B = B + (ST4(1)/(3*Veff)) * (nMT4(1, :)' * face_sum_1);

            face_sum_2 = N1 + N2 + N4;  % (1,2,4)
            B = B + (ST4(2)/(3*Veff)) * (nMT4(2, :)' * face_sum_2);

            face_sum_3 = N1 + N3 + N4;  % (1,3,4)
            B = B + (ST4(3)/(3*Veff)) * (nMT4(3, :)' * face_sum_3);

            face_sum_4 = N2 + N3 + N4;  % (2,3,4)
            B = B + (ST4(4)/(3*Veff)) * (nMT4(4, :)' * face_sum_4);

            % ------------------ 刚度累计（Poisson/扩散型） ------------------
            Ke = Ke + (B.' * B) * Veff;
        end
    end

    % 单元→全局装配
    gdof = uniqueNumbers(:);            % (E×1)
    K(gdof, gdof) = K(gdof, gdof) + Ke;
end
end
