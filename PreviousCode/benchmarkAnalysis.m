clc; close all; clear; tStart = cputime;


% --- Load the Bechmark
load("Benchmark\MoMbenchmarkData1.mat");

% --- Define the analysis axis
analysisAxis = r(:,3);

% --- Get the vertices of the geometry
A = triverts(1,:)';
B = triverts(2,:)';
C = triverts(3,:)';
v = [A, B, C];

% --- Compute the potencial integral 
Icalc = zeros(length(r),1);
for ii=1:length(r)
   Icalc(ii, :) = computeIntegral(r(ii,:)', A, B, C); 
end

% --- Plot Geometry and observation points

figure('Name','Observation points')
hold on;
patch(v(1,:), v(2,:), v(3,:), 'FaceColor', [0, 0, 0], 'linestyle', ':')
plot3(r(:, 1), r(:, 2), r(:, 3), '-r', 'MarkerFaceColor', 'r')
grid on, grid minor;
axis image;
view([-37 -43]);
camroll(-360);
title('Observation Points')
xlabel('x-axis [m]')
ylabel('y-axis [m]')
zlabel('z-axis [m]')

% --- Plot Results
figure('Name','Results')
subplot(2, 1, 1)
hold on
plot(analysisAxis, Ival, '--b', 'LineWidth', 1.2)
plot(analysisAxis, Icalc, ':r', 'LineWidth', 1.5)
legend("Benchmark", "Results"),
title("Results Comparison")
grid on, grid minor;
xlabel("Position [m]")
ylabel("Integral [dS/R]")
M = 100*(Ival-Icalc)./Ival;
[~, ind] = max(abs(M), [], 1, 'linear');

subplot(2, 1, 2)
hold on
plot(analysisAxis, M, 'b')
yline(M(ind), '--r', string(M(ind)),...
      'LabelHorizontalAlignment','center',...
      'LabelVerticalAlignment','middle');
legend("Error", "Max. error"),
title("Results Conformance")
xlabel("x position [m]")
ylabel("Porcentual error [%]")
grid on, grid minor;

%whos
simulationTime = cputime - tStart;
disp("Simulation time: "+num2str(simulationTime)+" s")

