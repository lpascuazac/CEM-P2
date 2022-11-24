function [V] = computeTheoricV(r,shape,arguments)
% --- Here it is computed the theoric potencial from an analytic result
% --- Some simple geometries have defined potencial equations, hence
% --- they work as a validation result from the computacional process.

    epsilon_0 = 8.8541878128E-12;
    ke = 1/(4*pi*epsilon_0);
    Q = arguments(1);
    switch(shape)
        case 'Sphere'
            R = arguments(2);
            if min(r)>R
                V = ke*Q./r;
            else
                V = ke*Q/R;
                warning("Some observation points are inside the sphere.")
            end
        case 'Disk'
            R = arguments(2);
            sigma = Q/(2*pi*R^2);
            %vDisk = Q*ke*log(abs(rDisk+sqrt(rDisk^2+R(3, :)'.^2))./(abs(R(3, :)')));
            V = 2*pi*ke*sigma*(sqrt(R^2+r.^2')-abs(r'));
        case 'infPlane'
            A = arguments(2);
            sigma = Q/A;
            V = -sigma*r/(2*epsilon_0);
        case 'Ring'
            R = arguments(2);
            V = ke*Q./sqrt(R^2+r.^2');
        otherwise
            disp('The entered figure does not have an analytical formula for the potential defined.')
    end
end

