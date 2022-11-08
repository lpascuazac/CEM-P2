clc; close all; clear; tStart = cputime;

% Input Paramenters
% --- Define the geometry
% --- Charge Q
% --- Geometry Arguments
% --- Geometry File

shape = 'Sphere';
Q = 1;
radius = 5;
arguments = [Q, radius];
obj = readObj('Shapes\Surf_5.obj');

% --- Define the observation points
r = [0, 0, 100]';
R = [zeros(100, 1), zeros(100, 1), linspace(10,100, 100)']';


% --- Get the .obj file parameters
meshVertices = obj.v(:, 1:3)';
faces = [obj.f.v]';

% --- Compute the vertices of each triangle
A = [meshVertices(:, faces(1,:))];
B = [meshVertices(:, faces(2,:))];
C = [meshVertices(:, faces(3,:))];

Icalc = zeros(size(R,2), 1);

% --- Compute the integral due to a uniform source distribution
for jj = 1:size(R,2)
    for ii=1:size(A, 2)
        Icalc(jj, :) = Icalc(jj, :) + computeIntegral(R(:, jj), A(:,ii), B(:,ii), C(:,ii));
    end
end

% --- Compute the potencial
epsilon_0 = 8.8541878128E-12;
ke = 1/(4*pi*epsilon_0);

theoricVoltage = computeTheoricV(R(3,:),shape,arguments); 
dq = Q/size(faces,2);
computePotencial = dq*ke*Icalc;

% --- Graphics

% ------ Shape
figure('Name','Simulation');
subplot(2, 2, 1)
hold on;
patch('Faces',faces','Vertices',meshVertices', 'FaceColor', [1 1 1])
view([45 35.264]);
axis image
grid on;
grid minor;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
zlabel('z-axis [m]')
title("Shape: "+string(size(faces,2)) + " triangles")

% ------ Observation points
subplot(2, 2, 2)
hold on;
patch('Faces',faces','Vertices',meshVertices', 'FaceColor', [1 1 1])
plot3([R(1,1), R(1,end)], [R(2,1), R(2,end)], [R(3,1), R(3,end)], '-or', 'MarkerFaceColor', 'r')
view([45 35.264]);
axis image
grid on;
grid minor;
title('Observation Points')
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')

% ------ Results Comparison
subplot(2,2,[3,4])
hold on
plot(R(3, :)', computePotencial, '.:b', 'LineWidth', 1.5)
plot(R(3, :)', theoricVoltage, '--r', 'LineWidth', 1.5)
legend("Result", "Theoric")
axis padded
grid on;
grid minor;
xlabel("r distance [m]")
ylabel("Voltage [V]"),
title("Results Comparison")

simulationTime = cputime - tStart;
disp("Simulation time: "+num2str(simulationTime)+" s")
