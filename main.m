%% Multi-geometries simulation

clc; close all; clear; tStart = cputime;

inputParameters = importdata('InputParameters.txt');

for m = 1:length(inputParameters.data)
    characteristics.("body"+m).boundaryCondition = inputParameters.data(m); 
end 

simulation = Simulation('Simulation Shapes',characteristics);
simulation.patchElements()
capacitanceMatrix = simulation.computeCapacitanceMatrix();

disp("Simulation time: "+num2str(cputime - tStart)+" s")

%% Potencial over an observation point

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
