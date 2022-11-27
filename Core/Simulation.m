classdef Simulation
    properties
        folder
        files
        bodies
        elements
        epsilon_0
        ke 
        totalTriangles
    end

    methods    
        function obj = Simulation(folder, characteristics)
            obj.files = dir(fullfile(folder,'*.obj'));
            obj.files = string({obj.files.name})';
            obj.folder = folder;
            obj.elements = length(obj.files);
            obj.epsilon_0  = 8.8541878128E-12;
            obj.ke = 1/(4*pi*obj.epsilon_0);
            obj.totalTriangles = 0;
            for i = 1:length(obj.files)
                obj.bodies.("body"+i) = PEC(obj.folder+"\"+obj.files(i),characteristics.("body"+i));
                obj.totalTriangles = obj.totalTriangles + obj.bodies.("body"+i).triangles;
            end

        end

        function [] = patchElements(obj)
            figure();
            hold on,
            FC = [1 1 1];

            for i = 1:obj.elements
                V = obj.bodies.("body"+i).geometry.v(:, 1:3);
                F = [obj.bodies.("body"+i).geometry.f.v, obj.bodies.("body"+i).geometry.f.v(:, 1)];
                patch('Faces',F,'Vertices',V, 'FaceColor', FC)
            end
            view([45 35.264]);
            axis equal
            grid on;
            grid minor;
        end
        
        function C = computeCapacitanceMatrix(obj)
            % Maxwell Capacitance Matrix
            areas = zeros(1, obj.totalTriangles);
            centers = zeros(3, obj.totalTriangles);
            vertex.A = zeros(3, obj.totalTriangles);
            vertex.B = zeros(3, obj.totalTriangles);
            vertex.C = zeros(3, obj.totalTriangles);

            position = 0;
            % Get the centers and areas of all the simulation triangles

            for M = 1:obj.elements
                areas(position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).areas*1e-6;

                centers(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).centers*1e-3;

                vertex.A(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.A*1e-3;

                vertex.B(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.B*1e-3;

                vertex.C(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.C*1e-3;

                position = obj.bodies.("body"+M).triangles;
            end
            
            % compute the MoM matrix
            MoM = zeros(obj.totalTriangles);

            for m = 1:obj.totalTriangles
                for n = 1:obj.totalTriangles
                    MoM(m,n) = obj.bodies.("body1")...
                        .computeOneIntegral(centers(:,m), ...
                        vertex.A(:,n), vertex.B(:,n), vertex.C(:,n), "")/(4*pi*obj.epsilon_0);
                end
            end
            
            % Compute charge density

            sigma = zeros(obj.totalTriangles, obj.elements);
            Q = zeros(obj.elements);
            position = 0;
            for M = 1:obj.elements
                voltages = zeros(obj.totalTriangles,1);
                voltages(position+1:obj.bodies.("body"+M).triangles+position) ...
                        = obj.bodies.("body"+M).characteristics.boundaryCondition;
                position = position + obj.bodies.("body"+M).triangles;
                sigma(:,M) = MoM\voltages;
            end

            % Add the contribution of each area
            position = 0;
            
%             ii = 0;
%             for M = 1:obj.elements
%                 areas_aux(M,:) = [areas(ii+1:obj.bodies.("body"+M).triangles+ii)];
%                 ii = obj.bodies.("body"+M).triangles;
%             end
%             
%             disp(areas_aux)
%             disp(size(areas_aux))

            for M = 1:obj.elements
                for N = 1:obj.totalTriangles
                    a = zeros(obj.elements,obj.totalTriangles);
                    a(1,position+1:obj.bodies.("body"+M).triangles+position) = ...
                        areas(position+1:obj.bodies.("body"+M).triangles+position);
                    a(2,position+1:obj.bodies.("body"+M).triangles+position) = ...
                        areas(position+1:obj.bodies.("body"+M).triangles+position);
                end
                position = position + obj.bodies.("body"+M).triangles;
                Q(:,M) = a(M,:)*sigma;
            end
            
            disp("Q")
            disp(Q)

            V = [1 0;0 1e-9];
            C = V\Q;
            disp("C")
            disp(C)

%             for n = 1:size(capMatrix,2)
%                 for m = 1:size(capMatrix, 2)
%                     chargeBodym = obj.bodies.("body"+m).characteristics.charge;
%                     voltageBodyn = obj.bodies.("body"+n).computeVoltage(GND);
%                         capMatrix(m, n) = chargeBodym/voltageBodyn;
%                         %disp(obj.geometries.("body"+n).characteristics.charge)
%                 end
%             end

            % Mutual capacitance
           
%             for n = 1:size(capMatrix,2)
%                 for m = 1:size(capMatrix, 2)
%                     charges(m) = obj.geometries.("body"+m).characteristics.charge;
%                 end
%                 voltageBodyn = obj.geometries.("body"+n).computeVoltage(GND);
%                 disp(voltageBodyn)
%                 capMatrix(:,n) = (charges/voltageBodyn)';
%             end


        end
    
    end

end
