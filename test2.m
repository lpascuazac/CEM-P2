clc; close all; clear; tStart = cputime;

% Input Paramenters
% --- Define the geometry
% --- Charge Q
% --- Geometry Arguments
% --- Geometry File

shape = 'Sphere';
Q = 1;
sphere.charge = Q;
sphere.radius = 5;

body1 = PEC('Shapes\Disk.obj');

meshVertices = body1.element.v(:, 1:3)';
faces = [body1.element.f.v]';

% --- Define the observation points
r = [0, 0, 100]';
R = [zeros(100, 1), zeros(100, 1), linspace(10,100, 100)']';


Icalc = body1.computeIntegral(R);

% --- Compute the potencial
epsilon_0 = 8.8541878128E-12;
ke = 1/(4*pi*epsilon_0);

theoricVoltage = body1.computeTheoricalV(R(3,:),shape,sphere); 
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