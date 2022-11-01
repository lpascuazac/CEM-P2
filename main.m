


clc;

clc

load("MoMbenchmark\MoMbenchmarkData.mat")

A = triverts(1,:)';
B = triverts(2,:)';
C = triverts(3,:)';

for ii=1:length(r)
   Icalc(ii, :) = triang_int(r(ii,:)', A, B, C); 
end

figure()
subplot(2, 1, 1)
hold on
plot(r(:,3), Ival, ':')
plot(r(:,3), Icalc, '-.')

subplot(2, 1, 2)
hold on
plot(r(:,3), Ival-Icalc)