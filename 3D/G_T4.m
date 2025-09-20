function [VT4, ST4, nMT4] = G_T4_single(nodes_tetra)
    % 功能：计算单个四面体的几何属性
    % 输入：
    % nodes_tetra - 4 × 3 矩阵，单个四面体的节点坐标。
    %               例如：
    %               nodes_tetra = [
    %                   x1 y1 z1; % 节点1
    %                   x2 y2 z2; % 节点2
    %                   x3 y3 z3; % 节点3
    %                   x4 y4 z4  % 节点4
    %               ];
    % 输出：
    % VT4   - 标量，四面体体积
    % ST4   - 1 × 4 向量，四个面的面积
    % nMT4  - 4 × 3 矩阵，四个面的单位外法向量

    % 检查输入维度
    if ~isequal(size(nodes_tetra), [4, 3])
        error('输入 nodes_tetra 必须是 4x3 的矩阵，表示单个四面体的四个节点坐标。');
    end

    % 提取节点坐标
    a = nodes_tetra(1, :); % 节点 1 (1x3)
    b = nodes_tetra(2, :); % 节点 2 (1x3)
    c = nodes_tetra(3, :); % 节点 3 (1x3)
    d = nodes_tetra(4, :); % 节点 4 (1x3)

    % 1. 计算体积 VT4
    ba = a - b;
    ca = a - c;
    da = a - d;
    cross_ca_da = cross(ca, da);
    VT4 = abs(dot(ba, cross_ca_da)) / 6;

    % 2. 计算面外法向量、面积
    % 四面体有 4 个面：(1,2,3), (1,2,4), (1,3,4), (2,3,4)
    face_nodes = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
    nT4_face = zeros(4, 3); % 临时存储面法向量，现在是 4x3 矩阵
    ST4 = zeros(1, 4);      % 1x4 向量，存储面积

    % 计算四面体质心
    centroids = mean(nodes_tetra, 1); 

    for f = 1:4
        p1_idx = face_nodes(f, 1);
        p2_idx = face_nodes(f, 2);
        p3_idx = face_nodes(f, 3);

        p1 = nodes_tetra(p1_idx, :);
        p2 = nodes_tetra(p2_idx, :); 
        p3 = nodes_tetra(p3_idx, :);
        
        % 计算法向量和面积
        v1 = p2 - p1;
        v2 = p3 - p1;
        normal = cross(v1, v2); % normal 是 1x3 向量
        
        ST4(f) = 0.5 * sqrt(sum(normal.^2));
        
        % 调整法向量方向，使其指向四面体外部
        to_centroid = centroids - p1;
        direction = dot(normal, to_centroid);
        
        if direction > 0 
            normal = -normal;
        end
        
        nT4_face(f, :) = normal; % 直接赋值给 nT4_face 的第 f 行
    end

    % 3. 计算单位外法向量 nMT4
    % nT4_face 现在是 4x3 矩阵
    % sum(nT4_face.^2, 2) 会对每行求平方和，得到 4x1 向量
    norms = sqrt(sum(nT4_face.^2, 2)); % norms 现在是 4x1 向量
    
    % nMT4 = nT4_face ./ norms;
    % MATLAB 会自动进行广播（implicit expansion），将 norms (4x1) 扩展为 4x3，
    % 然后逐元素除以 nT4_face (4x3)。
    nMT4 = nT4_face ./ norms; 
end
