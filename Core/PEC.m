classdef PEC
    properties
        file
        name
        geometry
        characteristics

        triangles
        areas
        centers
        vertex

    end

    methods
        
    function obj = PEC(file,characteristics)
            obj.file = file;
            obj.characteristics = characteristics;
            obj.geometry = readObj(obj.file);

            meshVertices = obj.geometry.v(:, 1:3)';
            faces = [obj.geometry.f.v]';
            
            obj.triangles = size(faces,2);

            obj.vertex.A = [meshVertices(:, faces(1,:))];
            obj.vertex.B = [meshVertices(:, faces(2,:))];
            obj.vertex.C = [meshVertices(:, faces(3,:))];
            
            obj.areas = zeros(obj.triangles,1);
            for n=1:size(obj.vertex.A, 2)
                v1 = obj.vertex.B(:,n)-obj.vertex.A(:,n);  
                v2 = obj.vertex.C(:,n)-obj.vertex.B(:,n);
                obj.areas(n) = norm(cross(v1,v2))/2;           
            end
            obj.centers = (obj.vertex.A + obj.vertex.B + obj.vertex.C)/3;
       
    end
    
    function voltage = computeVoltage(obj, observationPoint)
        Icalc = obj.computeIntegral(observationPoint, "unitCharge");
        dq = obj.characteristics.boundaryCondition/size(obj.geometry.f.v,1);
        epsilon_0  = 8.8541878128E-12;
        ke = 1/(4*pi*epsilon_0);
        voltage = dq*ke*Icalc;
    end

    function voltage = computeTheoricalVoltage(obj, observationPoint)
        % --- Here it is computed the theoric potencial from an analytic result
        % --- Some simple geometries have defined potencial equations, hence
        % --- they work as a validation result from the computacional process.        
        
        Q = obj.characteristics.boundaryCondition;
        epsilon_0  = 8.8541878128E-12;
        ke = 1/(4*pi*epsilon_0);

        switch(obj.characteristics.shape)
            case 'Sphere'
                R = obj.characteristics.radius;
                if min(observationPoint)>R
                    voltage = ke*Q./observationPoint;
                else
                    voltage = ke*Q/R;
                    warning("Some observation points are inside the sphere.")
                end
            case 'Disk'
                R = obj.characteristics.radius;
                sigma = Q/(2*pi*R^2);
                %vDisk = Q*obj.ke*log(abs(rDisk+sqrt(rDisk^2+R(3, :)'.^2))./(abs(R(3, :)')));
                voltage = 2*pi*ke*sigma*(sqrt(R^2+observationPoint.^2')-abs(observationPoint'));
            case 'infPlane'
                A = obj.characteristics.area;
                sigma = Q/A;
                voltage = -sigma*observationPoint/(2*epsilon_0);
            case 'Ring'
                R = obj.characteristics.radius;
                voltage = ke*Q./sqrt(R^2+observationPoint.^2');
            otherwise
                disp('The entered figure does not have an analytical formula for the potential defined.')
        end
    end
    
    function Icompute = computeIntegral(obj, observationPoint, format)
    % --- Here it is computed the integral equation from a uniform source
    % --- distribution, due to a triangular mesh grid
        Icompute = zeros(size(observationPoint,2), 1);

        % Define reference vectors
        meshVertices = obj.geometry.v(:, 1:3)';
        faces = [obj.geometry.f.v]';

        % --- Compute the vertices of each triangle
        A = [meshVertices(:, faces(1,:))];
        B = [meshVertices(:, faces(2,:))];
        C = [meshVertices(:, faces(3,:))];

        R = observationPoint;
        for jj = 1:size(observationPoint,2)
            for ii=1:size(A, 2)
                Icompute(jj, :) = Icompute(jj, :) + ...
                    obj.computeOneIntegral(R(:, jj), A(:,ii), B(:,ii), C(:,ii), format);
            end
        end
    end
        
    function Isum = computeOneIntegral(~, observationPoint, A, B, C, format)

        R = observationPoint;

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

        while 1
            rho = R - n*(n'*R);
            
            l = [(rho_p(:,1)-rho_m(:,1))/norm(rho_p(:,1)-rho_m(:,1)),...
                (rho_p(:,2)-rho_m(:,2))/norm(rho_p(:,2)-rho_m(:,2)),...
                (rho_p(:,3)-rho_m(:,3))/norm(rho_p(:,3)-rho_m(:,3))];
            
            u = [cross(l(:,1), n), cross(l(:,2), n), cross(l(:,3), n)];
            
            l_p = [(rho_p(:,1)-rho)'*l(:,1), (rho_p(:,2)-rho)'*l(:,2), (rho_p(:,3)-rho)'*l(:,3)];
            l_m = [(rho_m(:,1)-rho)'*l(:,1), (rho_m(:,2)-rho)'*l(:,2), (rho_m(:,3)-rho)'*l(:,3)];
            
            d = [-n'*(R(:,1)-l_mm(:,1))*n, -n'*(R-l_mm(:,2))*n, -n'*(R-l_mm(:,2))*n];
            
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
            switch (format)
                case "unitCharge"
                    Isum = 2*sum(I)/norm(cross(v1,v2));
                otherwise
                    Isum = sum(I);
                    
            end
    
            if isnan(Isum)
                R = R + 1e-9;
            else
                break;
            end
        end
    end
    
    function [] = plotResultPotencial(obj, observationPoint)
        % --- Graphics
        meshVertices = obj.geometry.v(:, 1:3)';
        faces = [obj.geometry.f.v];

        theoricVoltage = obj.computeTheoricalVoltage(observationPoint(3,:)); 

        computeVoltage = obj.computeVoltage(observationPoint(3,:));


        % ------ Shape
        figure('Name',obj.file);
        subplot(2, 2, 1)
        hold on;
        patch('Faces',faces,'Vertices',meshVertices', 'FaceColor', [1 1 1])
        view([45 35.264]);
        axis image
        grid on;
        grid minor;
        xlabel('x-axis [m]')
        ylabel('y-axis [m]')
        zlabel('z-axis [m]')
        title("Shape: "+string(size(faces',2)) + " triangles")
        
        % ------ Observation points
        subplot(2, 2, 2)
        hold on;
        patch('Faces',faces,'Vertices',meshVertices', 'FaceColor', [1 1 1])
        plot3([observationPoint(1,1), observationPoint(1,end)], ...
                [observationPoint(2,1), observationPoint(2,end)], ...
                [observationPoint(3,1), observationPoint(3,end)], '-or', 'MarkerFaceColor', 'r')
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
        plot(observationPoint(3, :)', computeVoltage, '.:b', 'LineWidth', 1.5)
        plot(observationPoint(3, :)', theoricVoltage, '--r', 'LineWidth', 1.5)
        legend("Result", "Theoric")
        axis padded
        grid on;
        grid minor;
        xlabel("r distance [m]")
        ylabel("Voltage [V]"),
        title("Results Comparison")
    end
    end

end
