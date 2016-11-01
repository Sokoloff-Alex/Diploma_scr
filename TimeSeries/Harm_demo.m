

clear all
close all
clc

T1 = 365.25;
T2=  351.00;
T3=  353.00;


A1 = 1;
A2 = 1;
A3 = 1;

time=0:1:365*(12.5*2);

f1 = A1*sin(2*pi/T1*(time)) + A1*cos(2*pi/T1*(time));
f2 = A2*sin(2*pi/T2*(time)) + A2*cos(2*pi/T2*(time));
f3 = A3*sin(2*pi/T3*(time)) + A3*cos(2*pi/T3*(time));

figure(1)
hold on
grid on
plot(time,f1,'.-b')
plot(time,f2,'.-r')
% plot(time,f3,'.-r')
plot(time,f2+f1,'.-g')
xlabel('days')
ylabel('amplitude')
legend('Annual, 365.25 days', 'Draconitic, GPS 351.0 days', 'combined')



