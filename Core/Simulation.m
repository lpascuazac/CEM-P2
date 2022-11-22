classdef Simulation
    properties
        folder
        files
        elements
        geometries
    end

    methods    
        function obj = Simulation(folder, characteristics)
            obj.files = dir(fullfile(folder,'*.obj'));
            obj.files = string({obj.files.name})';
            obj.folder = folder;
            obj.geometries = length(obj.files);
            for i = 1:length(obj.files)
                obj.elements.("body"+i) = PEC(obj.folder+"\"+obj.files(i),characteristics.("body"+i));
            end
        end
        


        function [] = patchElements(obj)
            figure();
            hold on,
            FC = [1 1 1];

            for i = 1:obj.geometries
                V = obj.elements.("body"+i).element.v(:, 1:3);
                F = [obj.elements.("body"+i).element.f.v, obj.elements.("body"+i).element.f.v(:, 1)];
                patch('Faces',F,'Vertices',V, 'FaceColor', FC)
            end
            axis equal
            grid on;
            grid minor;
        end
        
        function capMatrix = computeCapacitanceMatrix(obj)
            
            capMatrix = 2;
        end
    
    end

end