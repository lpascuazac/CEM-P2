%{
 copyfile(fullfile(matlabroot,'toolbox','shared','cmlink','git','auxiliary_files', ...
'mwgitattributes'),fullfile(pwd,'.gitattributes'))
%}
clc; clear; close all; 

%obj = readObj('examples\example10.obj')

obj = readObj('Shapes\Surf_4.obj');
%obj = readObj('Shapes\Surf_2.obj')
%obj = readObj('Shapes\Surf_4.obj')

%OBJ = read_wobj('examples\example1.obj')

V = obj.v(:, 1:3);
F = [obj.f.v]; 


%[x y] = meshgrid(-1:0.1:1); % Generate x and y data
%z = zeros(size(x, 1)); % Generate z data
%surf(x, y, z) % Plot the surface

%V = [1 1 1; 2 1 1; 2 2 0];
%F = [4, 3, 2, 4; 2, 1, 4, 1];
figure();
FC = [1 1 1];
patch('Faces',F,'Vertices',V, 'FaceColor', FC)
axis vis3d
axis padded
grid on;
grid minor;