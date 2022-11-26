clc; close all; clear; tStart = cputime;

% Input Paramenters
% --- Define the geometry
% --- Charge Q
% --- Geometry Arguments
% --- Geometry File

propertiesBody1.shape = 'Sphere'; 
propertiesBody1.charge = 1;
propertiesBody1.radius = 5;

body1 = PEC('Shapes\Surf_5.obj',propertiesBody1);


% --- Define the observation points
r = [0, 0, 100]';
R = [zeros(100, 1), zeros(100, 1), linspace(10,100, 100)']';


% --- Compute the potencial

body1.plotResultPotencial(R)


simulationTime = cputime - tStart;
disp("Simulation time: "+num2str(simulationTime)+" s")

%% Multi-geometries simulation

clc; close all; clear; tStart = cputime;
 

propertiesBody1.shape = 'Sphere';
propertiesBody1.boundaryCondition = 1;
propertiesBody1.radius = 4.2;

propertiesBody2.shape = 'Sphere';
propertiesBody2.boundaryCondition = -1;
propertiesBody2.radius = 4.2;

characteristics.body1 = propertiesBody1;
characteristics.body2 = propertiesBody2;

simulation = Simulation('Probe Shapes',characteristics);
simulation.patchElements()
capacitanceMatrix = simulation.computeCapacitanceMatrix();

%disp(capacitanceMatrix)

disp("Simulation time: "+num2str(cputime - tStart)+" s")

