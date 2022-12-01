%% Multi-geometries simulation

clc; close all; clear; tStart = cputime;
simulation = Simulation('Simulation Shapes','InputParameters.txt');
simulation.patchElements();
outputData = simulation.computeCapacitanceMatrix();

disp("Simulation time: "+num2str(cputime - tStart)+" s")
