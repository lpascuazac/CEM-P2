clc; close all; clear; tStart = cputime;

% Input Paramenters
% --- Define the geometry
% --- Charge Q
% --- Geometry Arguments
% --- Geometry File

body1.shape = 'Sphere'; 
body1.charge = 1;
body1.radius = 5;

body1 = PEC('Shapes\Surf_5.obj',body1);


% --- Define the observation points
r = [0, 0, 100]';
R = [zeros(100, 1), zeros(100, 1), linspace(10,100, 100)']';


% --- Compute the potencial

body1.plotResultPotencial(R)


simulationTime = cputime - tStart;
disp("Simulation time: "+num2str(simulationTime)+" s")

characteristics.body1 = body1; 
characteristics.body2 = body1;

simulation = Simulation('Probe Shapes',characteristics);
simulation.patchElements()

