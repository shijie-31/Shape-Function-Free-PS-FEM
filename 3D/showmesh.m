function showmesh(node,elem,options)
% Showmesh displays a mesh in 2-D and 3-D, ensuring each element has a consistent random color.

if nargin==2, options.FaceAlpha = 0.4; end
if nargin==3 && ~isfield(options,'FaceAlpha') 
    options.FaceAlpha = 0.4;
end

dim = size(node,2);
if ~iscell(elem)
    if size(elem,2)==3  % triangles
        h = patch('Faces', elem, 'Vertices', node);
    elseif size(elem,2)==4 % tetrahedrons
        h = tetramesh(elem,node,ones(size(elem,1),1));
    end
else
    if iscell(elem{1}), elem = vertcat(elem{:}); end  % Flatten the cell array if needed
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h = patch('Faces', tpad, 'Vertices', node);
end

if dim==3
    view(3); set(h,'FaceAlpha',options.FaceAlpha); % transparency
end

% Generate the colormap from blue to white to red
numColors = 1024;  % Number of colors
mid = floor(numColors / 2);  % Middle index

% Define the color stops: blue, white, and red
red = [88, 104, 161]/255;
white = [255, 255, 255]/255;
blue = [237, 164, 130]/255;


% Create a gradient from blue to white, then from white to red
colors = [linspace(blue(1), white(1), mid)', linspace(blue(2), white(2), mid)', linspace(blue(3), white(3), mid)'; ...
          linspace(white(1), red(1), mid)', linspace(white(2), red(2), mid)', linspace(white(3), red(3), mid)'];

% Assign the same color to each element
num_elems = size(elem, 1);  % Number of elements
color_indices = randi([1 numColors], num_elems, 1);  % Random color indices for each element
element_colors = colors(color_indices, :);  % Random colors for each element

% Handle 'elem' as a cell array if necessary
if iscell(elem)
    num_vertices = size(node, 1);  % Number of vertices
    vertex_colors = zeros(num_vertices, 3);  % Initialize vertex colors

    for i = 1:num_elems
        vertex_indices = elem{i};  % Extract vertex indices for this element
        vertex_colors(vertex_indices, :) = repmat(element_colors(i, :), length(vertex_indices), 1);  % Apply the same color to all vertices of the element
    end
else
    % If 'elem' is not a cell, proceed with regular array processing
    num_vertices = size(node, 1);  % Number of vertices
    vertex_colors = zeros(num_vertices, 3);  % Initialize vertex colors
    for i = 1:num_elems
        vertex_indices = elem(i, :);  % Indices of vertices for this element
        vertex_colors(vertex_indices, :) = repmat(element_colors(i, :), length(vertex_indices), 1);  % Apply the same color to all vertices of the element
    end
end

% Apply colors to the mesh
set(h, 'FaceVertexCData', vertex_colors, 'FaceColor', 'flat');  % Assign colors with no interpolation

% Optional: Apply default face color if not specified in options
if isfield(options,'facecolor') 
    facecolor = options.facecolor;
else
    facecolor = [204,226,219]/255;
end

set(h,'edgecolor','k');
axis equal; 
axis off;  % Hide the axis
