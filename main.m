

tStart = cputime;


clc;


load("MoMbenchmark\MoMbenchmarkData2.mat");


A = triverts(1,:)';
B = triverts(2,:)';
C = triverts(3,:)';

Icalc = zeros(length(r),1);

for ii=1:length(r)
   Icalc(ii, :) = triang_int(r(ii,:)', A, B, C); 
end

x_axis = r(:,1);

figure()
subplot(2, 1, 1)
hold on
plot(x_axis, Ival, '--')
plot(x_axis, Icalc, ':')
legend("Benchmark", "Calculation")
%grid on, grid minor;

M = Ival-Icalc;
[~, ind] = max(abs(M), [], 1, 'linear');


subplot(2, 1, 2)
hold on
plot(x_axis, (M))
yline(M(ind), '--r', string(M(ind)),...
      'LabelHorizontalAlignment','center',...
      'LabelVerticalAlignment','middle');
legend
%grid on, grid minor;

%whos
tEnd = cputime - tStart