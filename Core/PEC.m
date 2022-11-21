classdef PEC
    properties
        file
        name
        element
        normal
        p1,p2,p3
    end

    methods
        
    function obj = PEC(file)
            obj.file = file;
            
            %
            % obj = readObj(fname)
            %
            % This function parses wavefront object data
            % It reads the mesh vertices, texture coordinates, normal coordinates
            % and face definitions(grouped by number of vertices) in a .obj file 
            % 
            % INPUT: fname - wavefront object file full path
            %
            % OUTPUT: init.element.v - mesh vertices
            %       : init.element.vt - texture coordinates
            %       : init.element.vn - normal coordinates
            %       : init.element.f - face definition assuming faces are made of of 3 vertices
            %
            % Bernard Abayowa, Tec^Edge
            % 11/8/07
            
            % set up field types
            vertex = []; textureCoor = []; normCoor = []; f.v = []; f.vt = []; f.vn = [];
            
            fid = fopen(file);
            % parse .obj file 
            while 1    
                tline = fgetl(fid);
                if ~ischar(tline),   break,   end  % exit at end of file 
                 ln = sscanf(tline,'%s',1); % line type 
                 %disp(ln)
                switch ln
                    case 'v'   % mesh vertexs
                        vertex = [vertex; sscanf(tline(2:end),'%f')'];
                    case 'vt'  % texture coordinate
                        textureCoor = [textureCoor; sscanf(tline(3:end),'%f')'];
                    case 'vn'  % normal coordinate
                        normCoor = [normCoor; sscanf(tline(3:end),'%f')'];
                    case 'f'   % face definition
                        fv = []; fvt = []; fvn = [];
                        str = textscan(tline(2:end),'%s'); str = str{1};
                   
                       nf = length(findstr(str{1},'/')); % number of fields with this face vertices
            
            
                       [tok str] = strtok(str,'//');     % vertex only
                        for k = 1:length(tok) fv = [fv str2num(tok{k})]; end
                       
                        if (nf > 0) 
                        [tok str] = strtok(str,'//');   % add texture coordinates
                            for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
                        end
                        if (nf > 1) 
                        [tok str] = strtok(str,'//');   % add normal coordinates
                            for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
                        end
                         f.v = [f.v; fv]; f.vt = [f.vt; fvt]; f.vn = [f.vn; fvn];
                end
            end
            fclose(fid);
            
            % set up matlab object 
            obj.element.v = vertex; obj.element.vt = textureCoor; 
            obj.element.vn = normCoor; obj.element.f = f;
            
            obj.p1 = obj.element.v(obj.element.f.v(:,1),1:3);
            obj.p2 = obj.element.v(obj.element.f.v(:,2),1:3);
            obj.p3 = obj.element.v(obj.element.f.v(:,3),1:3);
            obj.normal = cross(obj.p2-obj.p1,obj.p3-obj.p1);
            obj.normal = obj.normal./norm(obj.normal);
        
    end
        
    function Icompute = computeIntegral(obj, observationPoint)
    % --- Here it is computed the integral equation from a uniform source
    % --- distribution, due to a triangular mesh grid
        Icompute = zeros(size(observationPoint,2), 1);

        % Define reference vectors
        meshVertices = obj.element.v(:, 1:3)';
        faces = [obj.element.f.v]';

        % --- Compute the vertices of each triangle
        A = [meshVertices(:, faces(1,:))];
        B = [meshVertices(:, faces(2,:))];
        C = [meshVertices(:, faces(3,:))];
        
        for jj = 1:size(observationPoint,2)
            for ii=1:size(A, 2)
                Icompute(jj, :) = Icompute(jj, :) + obj.computeOneIntegral(observationPoint(:, jj), A(:,ii), B(:,ii), C(:,ii));
            end
        end

    end
        
    function Isum = computeOneIntegral(~, observationPoint, A, B, C)
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
        
        rho = observationPoint - n*(n'*observationPoint);
        
        l = [(rho_p(:,1)-rho_m(:,1))/norm(rho_p(:,1)-rho_m(:,1)),...
            (rho_p(:,2)-rho_m(:,2))/norm(rho_p(:,2)-rho_m(:,2)),...
            (rho_p(:,3)-rho_m(:,3))/norm(rho_p(:,3)-rho_m(:,3))];
        
        u = [cross(l(:,1), n), cross(l(:,2), n), cross(l(:,3), n)];
        
        l_p = [(rho_p(:,1)-rho)'*l(:,1), (rho_p(:,2)-rho)'*l(:,2), (rho_p(:,3)-rho)'*l(:,3)];
        l_m = [(rho_m(:,1)-rho)'*l(:,1), (rho_m(:,2)-rho)'*l(:,2), (rho_m(:,3)-rho)'*l(:,3)];
        
        d = [-n'*(observationPoint(:,1)-l_mm(:,1))*n, -n'*(observationPoint-l_mm(:,2))*n, -n'*(observationPoint-l_mm(:,2))*n];
        
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
        
        Isum = 2*sum(I)/norm(cross(v1,v2));

end
    
    function V = computeTheoricalV(~,observationPoint,shape,characteristics)
        % --- Here it is computed the theoric potencial from an analytic result
        % --- Some simple geometries have defined potencial equations, hence
        % --- they work as a validation result from the computacional process.        
            epsilon_0 = 8.8541878128E-12;
            ke = 1/(4*pi*epsilon_0);
            Q = characteristics.charge;
            switch(shape)
                case 'Sphere'
                    R = characteristics.radius;
                    if min(observationPoint)>R
                        V = ke*Q./observationPoint;
                    else
                        V = ke*Q/R;
                        warning("Some observation points are inside the sphere.")
                    end
                case 'Disk'
                    R = characteristics.radius;
                    sigma = Q/(2*pi*R^2);
                    %vDisk = Q*ke*log(abs(rDisk+sqrt(rDisk^2+R(3, :)'.^2))./(abs(R(3, :)')));
                    V = 2*pi*ke*sigma*(sqrt(R^2+observationPoint.^2')-abs(observationPoint'));
                case 'infPlane'
                    A = characteristics.area;
                    sigma = Q/A;
                    V = -sigma*observationPoint/(2*epsilon_0);
                case 'Ring'
                    R = characteristics(2);
                    V = ke*Q./sqrt(R^2+observationPoint.^2');
                otherwise
                    disp('The entered figure does not have an analytical formula for the potential defined.')
            end
    end
    end
end
