function K = GK(node3, elem3)
%   --------------------------------------------------------------------
%   (c) 2025 Shijie Zhao <shijie31123@163.com>
%   --------------------------------------------------------------------
ndof = 1;                          % 每节点自由度
sdof = ndof * size(node3,1);        % 总自由度
K    = sparse(sdof, sdof);         % 全局刚度


for iel = 1 : numel(elem3)
    face   = elem3{iel};             
    nface    = numel(  face );

    %单元中心的信息
    allNumbers = [face{:}];
    uniqueNumbers = unique(allNumbers); 
    E=length(uniqueNumbers);
    C=  node3(uniqueNumbers,:);
    CV=mean (C);
    NV= ones(1, E)*   (1/E) ;    %NODE平均就可以

    %-------------

    for iface = 1 : nface
           nod=face{iface};
           num=length(nod);
           eNode = node3(nod,:); 
             NF = zeros(1, nface); % 首先创建一个 1 行 nface 列的零矩阵
              NF(nod) = 1/num;  
          
           
           Cface=mean(eNode);
          for it = 1 : num   %这里的num也会代表四面体的数量
              for j=1: 4
                e1 = eNode(j, :);
                e2_idx = mod(j, num) + 1;
                e2 = eNode(e2_idx, :);
              N1 = zeros(1, nface); % 首先创建一个 1 行 nface 列的零矩阵
              N1(j) = 1;  
              N2 = zeros(1, nface); % 首先创建一个 1 行 nface 列的零矩阵
              N2(e2_idx) = 1
              nodet4=[e1;e2;Cface;CV]    %四点坐标 E1 E2 CFACE CV

              [LT4, nT4, VT4, ST4, nMT4] = G_T4(nodet4)
              end
                 
              end
    end



end