function plot_residuals(resENU,Epoch)
% plot residuals plot and histogram at left side
%

figure(1)
subplot(3,3,[1,2])
hold on
grid on
plot(Epoch, resENU(:,2)*1000,'.')
legend(['Norts, res = ', num2str(rms(resENU(:,2))*1000,'%5.2f'), ' mm'])
ylim([-5 5])
ylabel('[mm]')

subplot(3,3,3)
hold on
grid on
[nelements,centers] = hist(resENU(:,2)*1000,20);
barh(centers,nelements)
ylim([-5 5])

subplot(3,3,[4,5])
hold on
grid on
plot(Epoch, resENU(:,1)*1000, '.')
legend(['East, RMS = ', num2str(rms(resENU(:,1))*1000,'%5.2f'), ' mm'])
ylim([-5 5])
ylabel('[mm]')

subplot(3,3,6)
hold on
grid on
[nelements,centers] = hist(resENU(:,2)*1000,20);
barh(centers,nelements)
ylim([-5 5])

subplot(3,3,[7,8])
hold on
grid on
plot(Epoch, resENU(:,3)*1000,'.')
legend(['Up, RMS = ', num2str(rms(resENU(:,3))*1000,'%5.2f'), ' mm'])
ylim([-20 20])
ylabel('[mm]')
xlabel('Epoch')

subplot(3,3,9)
hold on
grid on
[nelements,centers] = hist(resENU(:,3)*1000,20);
barh(centers, nelements)
ylim([-20 20])

