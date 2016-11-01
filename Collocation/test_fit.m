% fit_test
% 6,44, 100 km

clc
close all

% x = [0      200     300     400     500     600     700:100:1000]';  
% y = [0.8    0       -0.2    -0.2    -0      0       -0.2*ones(size(700:100:1000))]';

x = [ 0         50          100     150     200      250:50:500]';
y = [0.0331   -0.0010   -0.0067    0.0011   -0.0067     zeros(size(250:50:500))]';

% x = [ 0       33.3333   66.6667  100.0000 ]';
% y = [0.0045   -0.0005   -0.0020    0.0009 ]';

b = +abs(min(y));
coeff_NN  = fitExp([y(1)  , -0.001], x, y);
coeff_NN2 = fitExp([y(1)+b, -0.001], x, y+b);

f_NN  = fit(x,y,  'exp1','StartPoint',[y(1),   -0.001]);
f_NN2 = fit(x,y+b,'exp1','StartPoint',[y(1)+b, -0.001]);
coeff_NN3 = [f_NN.a  f_NN.b];
coeff_NN4 = [f_NN2.a f_NN2.b];

[K0_1,  a_1] = fitCov([y(1),100], x, y)
[K0_2,  a_2] = fitCov([y(1),100], x, y+b)

d = 0:1:max(x);

clr = lines(8);
figure(4)
hold on
grid on
% title(['point : lat ', num2str(lat0), ', long  ' num2str(long0), ' # obs.: ', num2str(length(lat)), ' ; step = ', num2str(scale) , ' km', ' # Classes ', num2str(nClasses)])
plot(x, y, 'o-b');
plot(x, y+b, 'o-g');

plot(d, y(1).*exp(-0.005*d ), '--k');
plot(d, y(1).*exp(-0.05*d), '--k');

%plot(d, coeff_NN(1).*exp(coeff_NN(2)*d),'--b');
plot(d, coeff_NN2(1).*exp(coeff_NN2(2)*d),'-g');
%plot(d, coeff_NN(1).*exp(coeff_NN2(2)*d),'--r');
pl1 = plot(d, y(1).*exp(coeff_NN2(2)*d),'-m');

plot(d, coeff_NN3(1).*exp(coeff_NN3(2)*d),'--r');
plot(d, coeff_NN4(1).*exp(coeff_NN4(2)*d),'.--r');
plot(d, coeff_NN3(1).*exp(coeff_NN4(2)*d),'-r');

% plot(d, K0_1./(1+(d/a_1).^2), '.-k' )
% plot(d, K0_1./(1+(d/20).^2), '-k')
plot(d, y(1)./(1+(d/a_2).^2), '.-m')



legend([pl1 ],['C_N_N  a = ',num2str(coeff_NN(1)), ' b = ', num2str(coeff_NN2(2)) ] )
hold off
                     