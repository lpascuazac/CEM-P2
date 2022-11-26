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
        
        function MoM = computeCapacitanceMatrix(obj)
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
                    = obj.bodies.("body"+M).areas;

                centers(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).centers;

                vertex.A(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.A;

                vertex.B(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.B;

                vertex.C(:,position+1:obj.bodies.("body"+M).triangles+position)...
                    = obj.bodies.("body"+M).vertex.C;

                position = obj.bodies.("body"+M).triangles;
            end
            
            % compute the MoM matrix

            MoM = zeros(obj.totalTriangles);

            for m = 1:obj.totalTriangles
                for n = 1:obj.totalTriangles
                    MoM(m,n) = obj.bodies.("body1")...
                        .computeOneIntegral(centers(:,m), ...
                        vertex.A(:,n), vertex.B(:,n), vertex.C(:,n), "")/4*pi*obj.epsilon_0;
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
                sigma(:,M) = MoM\voltages;
            end

            % Add the contribution of each area
            
            for M = 1:obj.elements
                for N = 1:obj.totalTriangles
                    a = zeros(obj.elements,obj.totalTriangles);
                    a(:,N) = areas(:,N);
                end
                Q(:,M) = a(M,:)*sigma;
            end
            
            disp("Q")
            disp(Q)

            V = [1 2;2 -1];
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

        function Isum = computeOneIntegral(~, r, A, B, C)
            v1 = B-A;  
            v2 = C-B;
            
            % Compute normal vector
            n = cross(v1, v2)/norm(cross(v1, v2));
            
            % Define: 
            % - The positive and negative sections of the segment of each triangle
            % - The position vectors (r)
            % - The pojection of the position vectors (rho)
            % - The paralell vector to the segment C (l)
            % - The normal vector to l (u)
            % - The distance from the observation point to the plane of each triangle (d)
            % - The vectors P0 and R0
            % - The distance P, P+- and R, R+-
            
            l_mm = [A, B, C];
            l_pp = [B, C, A];
            
            r_p = l_pp;
            r_m = l_mm;
            
            rho_p = r_p - n*(n'*r_p);
            rho_m = r_m - n*(n'*r_m);
            
            rho = r - n*(n'*r);
            
            l = [(rho_p(:,1)-rho_m(:,1))/norm(rho_p(:,1)-rho_m(:,1)),...
                (rho_p(:,2)-rho_m(:,2))/norm(rho_p(:,2)-rho_m(:,2)),...
                (rho_p(:,3)-rho_m(:,3))/norm(rho_p(:,3)-rho_m(:,3))];
            
            u = [cross(l(:,1), n), cross(l(:,2), n), cross(l(:,3), n)];
            
            l_p = [(rho_p(:,1)-rho)'*l(:,1), (rho_p(:,2)-rho)'*l(:,2), (rho_p(:,3)-rho)'*l(:,3)];
            l_m = [(rho_m(:,1)-rho)'*l(:,1), (rho_m(:,2)-rho)'*l(:,2), (rho_m(:,3)-rho)'*l(:,3)];
            
            d = [-n'*(r(:,1)-l_mm(:,1))*n, -n'*(r-l_mm(:,2))*n, -n'*(r-l_mm(:,2))*n];
            
            P0 = [norm((rho_p(:,1)-rho)'*u(:,1)), norm((rho_p(:,2)-rho)'*u(:,2)), norm((rho_p(:,3)-rho)'*u(:,3))];
            
            P0_n = [((rho_p(:,1)-rho)-l_p(:,1)*l(:,1))/P0(:,1),...
                ((rho_p(:,2)-rho)-l_p(:,2)*l(:,2))/P0(:,2),...
                ((rho_p(:,3)-rho)-l_p(:,3)*l(:,3))/P0(:,3)];
            
            R0 = [sqrt(P0(:,1).^2 + norm(d(:,1))^2),...
                sqrt(P0(:,2).^2 + norm(d(:,2))^2),...
                sqrt(P0(:,3).^2 + norm(d(:,3))^2)];
            
            P_p = [norm(rho_p(:,1)-rho), norm(rho_p(:,2)-rho), norm(rho_p(:,3)-rho)];
            P_m = [norm(rho_m(:,1)-rho), norm(rho_m(:,2)-rho), norm(rho_m(:,3)-rho)];
            
            R_p = [sqrt(P_p(:,1)^2 + norm(d(:,1))^2),...
                sqrt(P_p(:,2)^2 + norm(d(:,2))^2),...
                sqrt(P_p(:,3)^2 + norm(d(:,3))^2)];
            
            R_m = [sqrt(P_m(:,1)^2 + norm(d(:,1))^2),...
                sqrt(P_m(:,2)^2 + norm(d(:,2))^2),...
                sqrt(P_m(:,3)^2 + norm(d(:,3))^2)];
            
            % --- Calculation of terms involved in the integral
            T1 = [P0(:,1)*log((R_p(:,1)+l_p(:,1))/(R_m(:,1)+l_m(:,1))),...
                P0(:,2)*log((R_p(:,2)+l_p(:,2))/(R_m(:,2)+l_m(:,2))),...
                P0(:,3)*log((R_p(:,3)+l_p(:,3))/(R_m(:,3)+l_m(:,3)))];
            
            T2 = [atan((P0(:,1)*l_p(:,1))/(R0(:,1)^2 + norm(d(:,1))*R_p(:,1))),...
                atan((P0(:,2)*l_p(:,2))/(R0(:,2)^2 + norm(d(:,2))*R_p(:,2))),...
                atan((P0(:,3)*l_p(:,3))/(R0(:,3)^2 + norm(d(:,3))*R_p(:,3)))];
            
            T3 = [atan((P0(:,1)*l_m(:,1))/(R0(:,1)^2 + norm(d(:,1))*R_m(:,1))),...
                atan((P0(:,2)*l_m(:,2))/(R0(:,2)^2 + norm(d(:,2))*R_m(:,2))),...
                atan((P0(:,3)*l_m(:,3))/(R0(:,3)^2 + norm(d(:,3))*R_m(:,3)))];
            
            I = [P0_n(:,1)'*u(:,1)*(T1(:,1)-norm(d(:,1))*(T2(:,1)-T3(:,1))),...
                P0_n(:,2)'*u(:,2)*(T1(:,2)-norm(d(:,2))*(T2(:,2)-T3(:,2))),...
                P0_n(:,3)'*u(:,3)*(T1(:,3)-norm(d(:,3))*(T2(:,3)-T3(:,3)))];
            
            % --- Finally the sum of the results of the apport of each triangle segment
            
            Isum = sum(I);
            
    end
    
    
    end

end
