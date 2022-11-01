
clc;

r = [0.5, 2, 1]';

A = [1 0 0]';
B = [0 1 0]';
C = [1 1 0]';

v1 = B-A;
v2 = C-B;

n = cross(v1, v2)/norm(cross(v1, v2));

l_mm = [A, B, C];
l_pp = [B, C, A];

r_p = l_pp;
r_m = l_mm;

rho_p = r_p - n*(n'*r_p);
rho_m = r_m - n*(n'*r_m);

rho = r - n*(n'*r);

%l = (rho_p-rho_m)/norm(rho_p-rho_m)
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

% Integral terms calculation
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
     P0_n(:,3)'*u(:,3)*(T1(:,3)-norm(d(:,3))*(T2(:,3)-T3(:,3)))]

I_sum = sum(I)