
%clc

A = [1 0 0]';
B = [0 1 0]';
C = [1 1 0]';

%pnt = [0.5, -0.5, 0.5]'; 

V = [A'; B'; C';];
F = [1 2 3;];

v1 = B-A;
v2 = C-B;

x0 = A(1);
y0 = A(2);
z0 = A(3);

A0 = [x0; y0; z0];

n = cross(v1, v2)/norm(cross(v1, v2));
N = cross(v1, v2);

Op = n'*A0*n;
%pnt_p = ((n'*A0 - pnt)'*n)*n
%pnt_p = (pnt - (dot(pnt,n)/norm(n)^2)*n);
%pnt_p2 = pnt - pnt*n.*n

[P,Q] = meshgrid(0:0.1:1); % Provide a gridwork (you choose the size)
X = x0+v1(1,1)*P+v2(1,1)*Q; % Compute the corresponding cartesian coordinates
Y = y0+v1(2,1)*P+v2(2,1)*Q; %   using the two vectors in w
Z = z0+v1(3,1)*P+v2(3,1)*Q;
%figure()z
%surf(X,Y,Z)


l_mm = C;
l_pp = A;

r = [0.5, 2, 0]';

r_p = l_pp;
r_m = l_mm;

rho_p = r_p - n*(n'*r_p);
rho_m = r_m - n*(n'*r_m);

rho = r - n*(n'*r);

l = (rho_p-rho_m)/norm(rho_p-rho_m);

u = cross(l, n);

l_p = (rho_p-rho)'*l;
l_m = (rho_m-rho)'*l;

d = -n'*(r-l_mm)*n;

P0 = norm((rho_p-rho)'*u);
P0_n = ((rho_p-rho)-l_p*l)/P0;

R0 = sqrt(P0^2 + norm(d)^2);

P_p = norm(rho_p-rho);
P_m = norm(rho_m-rho);

R_p = sqrt(P_p^2 + norm(d)^2); 
R_m = sqrt(P_m^2 + norm(d)^2); 

% Integral terms calculation
T1 = P0*log((R_p+l_p)/(R_m+l_m));
T2 = atan((P0*l_p)/(R0^2 + norm(d)*R_p));
T3 = atan((P0*l_m)/(R0^2 + norm(d)*R_m));

I = P0_n'*u*(T1-norm(d)*(T2-T3))

figure();
FC = [1 1 1];
hold on;
plot3(0, 0, 0, 'ok', 'MarkerFaceColor', 'k')
plot3(A(1), A(2), A(3), 'or')
%plot3(Aq(1), Aq(2), Aq(3), 'or', 'MarkerFaceColor', 'r')

plot3(B(1), B(2), B(3), 'ob')
%plot3(Bq(1), Bq(2), Bq(3), 'ob', 'MarkerFaceColor', 'b')

plot3(C(1), C(2), C(3), 'om')
%plot3(Cq(1), Cq(2), Cq(3), 'om', 'MarkerFaceColor', 'm')

%plot3(P1(1), P1(2), P1(3), 'ok')
%plot3(Pq(1), Pq(2), Pq(3), 'ok', 'MarkerFaceColor', 'k')

%plot3([A(1), n(1)+A(1)], [A(2), n(2)+A(2)], [A(3), n(3)+A(3)], 'r')
plot3(Op(1), Op(2), Op(3), 'ok', 'MarkerFaceColor', 'none')
plot3(r(1), r(2), r(3), 'r', 'Marker', 'hexagram', 'MarkerFaceColor', 'r')
%plot3(pnt_p(1)+Op, pnt_p(2)+Op, pnt_p(3)+Op, 'r', 'Marker', 'hexagram', 'MarkerFaceColor', 'none')

plot3([0, n(1)], [0, n(2)], [0, n(3)], ':r')
%plot3([pnt(1), n(1)+pnt(1)], [pnt(2), n(2)+pnt(2)], [pnt(3), n(3)+pnt(3)], ':r')

plot3([Op(1), rho_p(1)+Op(1)], [Op(2), rho_p(2)+Op(2)], [Op(3), rho_p(3)+Op(3)], ':b')
plot3([Op(1), rho_m(1)+Op(1)], [Op(2), rho_m(2)+Op(2)], [Op(3), rho_m(3)+Op(3)], '--b')

plot3([l_mm(1), l_pp(1)], [l_mm(2), l_pp(2)], [l_mm(3), l_pp(3)], '-.m', LineWidth=1.5)

plot3([0, r(1)], [0, r(2)], [0, r(3)], '-.r')
plot3([Op(1), rho(1)+Op(1)], [Op(2), rho(2)+Op(2)], [Op(3), rho(3)+Op(3)], '-.r')

%plot3([pnt_p(1), d(1)+pnt_p(1)], [pnt_p(2), d(2)+pnt_p(2)], [pnt_p(3), d(3)+pnt_p(3)], '-.k')
plot3([r(1), d(1)+r(1)], [r(2), d(2)+r(2)], [r(3), d(3)+r(3)], '-.k')
plot3(rho(1)+Op(1), rho(2)+Op(2), rho(3)+Op(3), 'r', 'Marker', 'hexagram', 'MarkerFaceColor', 'none')

%plot3([r(1), d(1)+r(1)], [r(2), d(2)+r(2)], [r(3), d(3)+r(3)], '-.k')
%plot3([rho(1), P0(1)+rho(1)], [rho(2), P0(2)+rho(2)], [rho(3), P0(3)+rho(3)], '--g')
%plot3([0, P0(1)], [0, P0(2)], [0, P0(3)], '--g')
%plot3([P0_n(1)+Op(1)], [P0_n(2)+Op(2)], [P0_n(3)+Op(3)], 'sb')
plot3([rho(1)+Op(1), P0_n(1)+rho(1)+Op(1)], [rho(2)+Op(2), P0_n(2)+rho(2)+Op(2)], [rho(3)+Op(3), P0_n(3)+rho(3)+Op(3)], ':b')

patch('Faces',F,'Vertices',V, 'FaceColor', FC)

surf(X,Y,Z, FaceAlpha=0.3, FaceColor="none", EdgeColor="k", LineStyle=":")
%shading flat
%legend("A", "A'", "B", "B'", "C", "C'")
legend("O", "A", "B", "C", "O'")

%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

axis equal
%axis vis3d
%axis padded
grid on;
grid minor;
xlabel('x')
ylabel('y')
zlabel('z')