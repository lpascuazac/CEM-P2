classdef Simulation
    properties
        folder
        files
        geometries
        elements
        epsilon_0
        ke 
    end

    methods    
        function obj = Simulation(folder, characteristics)
            obj.files = dir(fullfile(folder,'*.obj'));
            obj.files = string({obj.files.name})';
            obj.folder = folder;
            obj.elements = length(obj.files);
            obj.epsilon_0  = 8.8541878128E-12;
            obj.ke = 1/(4*pi*obj.epsilon_0);
            for i = 1:length(obj.files)
                obj.geometries.("body"+i) = PEC(obj.folder+"\"+obj.files(i),characteristics.("body"+i));
            end
        end

        function [] = patchElements(obj)
            figure();
            hold on,
            FC = [1 1 1];

            for i = 1:obj.elements
                V = obj.geometries.("body"+i).geometry.v(:, 1:3);
                F = [obj.geometries.("body"+i).geometry.f.v, obj.geometries.("body"+i).geometry.f.v(:, 1)];
                patch('Faces',F,'Vertices',V, 'FaceColor', FC)
            end
            view([45 35.264]);
            axis equal
            grid on;
            grid minor;
        end
        
        function capMatrix = computeCapacitanceMatrix(obj)
            capMatrix = zeros(obj.elements, obj.elements);
            GND = [0, 0, 0]';

            % Maxwell Capacitance Matrix
%             for n = 1:size(capMatrix,2)
%                 for m = 1:size(capMatrix, 2)
%                     voltageBodym = obj.geometries.("body"+m).computeVoltage(GND);
%                     voltageBodyn = obj.geometries.("body"+n).computeVoltage(GND);
%                     if m == n
%                         capMatrix(m, n) = obj.geometries.("body"+m).characteristics.charge...
%                                           /(voltageBodym);
%                     else
%                         capMatrix(m, n) = obj.geometries.("body"+m).characteristics.charge...
%                                           /voltageBodyn;
%                         %disp(obj.geometries.("body"+n).characteristics.charge)
%                     end
%                 end
%             end

            % Mutual capacitance
            charges = zeros(1,obj.elements);
            for n = 1:size(capMatrix,2)
                for m = 1:size(capMatrix, 2)
                    charges(m) = obj.geometries.("body"+m).characteristics.charge;
                end
                voltageBodyn = obj.geometries.("body"+n).computeVoltage(GND);
                capMatrix(:,n) = (charges/voltageBodyn)';
            end


        end
    
    end

end
