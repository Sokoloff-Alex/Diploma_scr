% Linearization

x = [0  1   2       3    4       5       6]';
y = [1  0.2 0.001   0.1 00.1    0.05    0.001]';

%%
[ones(length(x),1), x] \ [log(y)]

%% LS
clc
a = (sum(log(y))*sum(x.^2) - sum(x)*sum(x.*log(y))) / (length(x)*sum(x.^2) - (sum(x.^2))^2)
b = (length(x)*sum(x.*log(y)) - sum(x)*sum(log(y))) / (length(x)*sum(x.^2) - (sum(x.^2))^2)

%% LSW
a1 = ( (x.^2)'*y  *      y'*log(y) - x'*y * (x.*y)'*log(y) )/( sum(y) * (x.^2)'*y - (x'*y)^2 )
b1 = (     sum(y) * (x.*y)'*log(y) - x'*y *      y'*log(y) )/( sum(y) * (x.^2)'*y - (x'*y)^2 )

c = fitExp([y(1), -0,01], x, y)

%%
x1 = 0:0.1:5
close all
figure(1)
hold on; grid on
plot(x,y,'--ob')
plot(x1, 1*exp(-0.33*x1) , '-.r' )
plot(x1, a*exp(b*x1), '-.m')
plot(x1, a1*exp(b1*x1), '-.k')
plot(x1, y(1)*exp(b1*x1), '-.g')
plot(x1, c(1)*exp(c(2)*x1), '-.m')
hold off


%%

