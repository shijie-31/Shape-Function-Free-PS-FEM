function showsolution3D(node,elem,u)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
%
% Copyright (C) Terence Yu.

if iscell(elem)
    if iscell(elem{1}), elem = vertcat(elem{:}); end
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    patch('Faces', tpad, 'Vertices', node,'facevertexcdata',u, 'facecolor','interp');
end
grid off;


view(3); 
grid on; 
 view(150,30);

%colormap('parula');  % 适合科学绘图的内置色彩图
numColors = 1024;  % 色彩图的分辨率
mid = floor(numColors/2);

% 从蓝色到白色
%red = [0.7, 0, 0];
%white = [0.8, 0.9, 1];
%blue = [0, 0, 0.7];

red = [88, 104, 161]/255;
white = [255, 255, 255]/255;
blue = [237, 164, 130]/255;

blue_to_white = [linspace(blue(1), white(1), mid)', ...
                linspace(blue(2), white(2), mid)', ...
                linspace(blue(3), white(3), mid)'];

% 从白色到红色
white_to_red = [linspace(white(1), red(1), numColors - mid)', ...
               linspace(white(2), red(2), numColors - mid)', ...
               linspace(white(3), red(3), numColors - mid)'];

% 合并两个部分
coolwarm = [blue_to_white; white_to_red];

colormap(coolwarm );
%camlight('headlight');
%lighting gouraud;  % 使用 Gouraud 光照实现平滑阴影

axis equal



