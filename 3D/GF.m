function Q = GQ(node3, elem3, qfun)
%---------------------------------------------------------------
%  (c) 2025  Shijie Zhao  <shijie31123@163.com>
%  PS-FEM 右端向量装配：多面体→四面体（边+面心+单元心）扇形剖分
%
%  Q = GQ(node3, elem3, qfun)
%  --------------------------------------------------------------
%  输入:
%     node3 : (Nnode × 3) 全局结点坐标
%     elem3 : {Nelement × 1} 元胞, 每元胞存该多面体的若干面,
%             每个面的顶点索引行向量
%     qfun  : 热源项, 可以是
%             • 常数        ->  均匀热源
%             • @(x,y,z)    ->  返回标量
%  输出:
%     Q     : (Nnode × 1)  全局热源向量
%---------------------------------------------------------------
ndof    = 1;                           % 每节点自由度
Nnode   = size(node3,1);
Q       = zeros(Nnode,1);              % 全局右端
vol_eps = 1e-14;                       % 体积容差

% 若 qfun 不是函数把它包成常数函数
if ~isa(qfun,'function_handle')
    qconst = qfun;
    qfun   = @(x,y,z) qconst;
end

for iel = 1:numel(elem3)
    faces         = elem3{iel};                 % 该多面体所有面
    allNumbers    = [faces{:}];
    uniqueNumbers = unique(allNumbers,'stable');% 单元涉及的结点
    E             = numel(uniqueNumbers);
    Ce            = node3(uniqueNumbers,:);     % (E×3)
    CV            = mean(Ce,1);                 % 单元心
    fe            = zeros(E,1);                 % 局部热源

    N4 = (1/E)*ones(1,E);                       % N4: 全域均分

    % ------------ 遍历各面 --------------
    for iface = 1:numel(faces)
        nod   = faces{iface}(:)';               % 本面的顶点编号
        numv  = numel(nod);
        eNode = node3(nod,:);                   % (numv×3)
        Cface = mean(eNode,1);                  % 面心

        % N3: 本面均分
        N3 = zeros(1,E);
        [~,loc_face] = ismember(nod,uniqueNumbers);
        loc_face = loc_face(loc_face>0);
        if ~isempty(loc_face)
            N3(loc_face) = 1/numv;
        end

        % ----------- 沿边形成四面体 -----------
        for j = 1:numv
            i1 = j;
            i2 = mod(j,numv)+1;
            e1 = eNode(i1,:);
            e2 = eNode(i2,:);

            % N1/N2: 边端点指示
            j1 = nod(i1);      j2 = nod(i2);
            N1 = zeros(1,E);   N2 = N1;
            [~,idx1] = ismember(j1,uniqueNumbers);
            [~,idx2] = ismember(j2,uniqueNumbers);
            if idx1>0, N1(idx1)=1; end
            if idx2>0, N2(idx2)=1; end

            % 构造四面体  [e1, e2, Cface, CV]
            X4 = [e1; e2; Cface; CV];           % (4×3)
            [VT4,~,~] = G_T4(X4);               % 体积
            Veff = abs(VT4);
            if Veff <= vol_eps,  continue;  end % 退化跳过

            % ---------- 1 点高斯积分 ----------
            Xc   = mean(X4,1);                  % 四面体重心
            qval = qfun(Xc(1),Xc(2),Xc(3));     % 热源值
            Nc   = (N1 + N2 + N3 + N4)/4;       % 形函数在重心

            fe = fe + (Nc.' * (qval * Veff));   % 累加到局部
        end
    end

    % ----------- 单元→全局装配 -------------
    gdof     = uniqueNumbers(:);
    Q(gdof)  = Q(gdof) + fe;
end
end
