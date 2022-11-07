

clc;

obj = readObj('Test_shapes\Surf_5.obj');

V = obj.v(:, 1:3)';
F = [obj.f.v]';

epsilon_0 = 8.8541878128E-12;
ke = 1/(4*pi*epsilon_0);

r = [0, 0, 1000]';

R = [zeros(100, 1), zeros(100, 1), linspace(5, 100, 100)']';

FC = [1 1 1];

figure();
subplot(1, 2, 1)
hold on;
patch('Faces',F','Vertices',V', 'FaceColor', FC)
view([45 35.264]);
axis image
grid on;
grid minor;
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
title(string(size(F,2)) + " triangles figure")

subplot(1, 2, 2)
hold on;

patch('Faces',F','Vertices',V', 'FaceColor', FC)
%plot3([0, r(1)], [0, r(2)], [0, r(3)], 'or', 'MarkerFaceColor', 'r')
%plot3(r(1), r(2), r(3), 'or', 'MarkerFaceColor', 'r')
plot3([R(1,1), R(1,end)], [R(2,1), R(2,end)], [R(3,1), R(3,end)], '-or', 'MarkerFaceColor', 'r')
%axis vis3d
%axis padded
view([45 35.264]);
%camroll(-360);
axis image
grid on;
grid minor;
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')

A = [V(:, F(1,:))];
B = [V(:, F(2,:))];
C = [V(:, F(3,:))];

Icalc = zeros(size(R,2), 1);

for jj = 1:size(R,2)
    for ii=1:size(A, 2)
        Icalc(jj, :) = Icalc(jj, :) + triang_int(R(:, jj), A(:,ii), B(:,ii), C(:,ii));
    end
end

Q = 1/size(F,2);

figure()
hold on
plot(R(3, :)', Q*Icalc, '.:b', 'LineWidth', 1.5)
plot(R(3, :)', 1./R(3, :)', '--r', 'LineWidth', 1.5)
legend("Calc.", "1/r ref.")
axis padded
grid on;
grid minor;
xlabel("r distance [m]")
ylabel("Integral [(dS/R)*(1/#Triangles)]")

%%
figure()
hold on
plot(R(3, :)', 1*Icalc*Q)
plot(R(3, :)', 1*1./R(3, :)')


figure()
hold on
plot(R(3, :)', Icalc./max(Icalc))
plot(R(3, :)', 1./(R(3, :)') .* (1/max(1./(R(3, :)'))) )


% for ii=1:size(A, 2)
%    Icalc(ii, :) = triang_int(r, A(:,ii), B(:,ii), C(:,ii)); 
% end

%%
epsilon_0 = 8.8541878128E-12;
ke = 1/(4*pi*epsilon_0);

I = sum(Icalc)

V = ke/norm(r)
% mk = ["o", "x", "+"];
% 
% for ii = 1:size(A, 2)
%     plot3(A(1,ii), A(2,ii), A(3,ii), mk(mod(ii, 3)+1)+"r", 'MarkerSize', 12)
%     plot3(B(1,ii), B(2,ii), B(3,ii), mk(mod(ii, 3)+1)+"b", 'MarkerSize', 12)
%     plot3(C(1,ii), C(2,ii), C(3,ii), mk(mod(ii, 3)+1)+"c", 'MarkerSize', 12)
% end


