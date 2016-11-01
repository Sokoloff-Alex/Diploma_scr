%% cov
clc
a = 1;
b = -0.04;
x = [a, b];

Max_Dist = 500; % km
for i = 1:4 % len
    point = [mean(long(Sets{i})),  mean(lat(Sets{i}))];
    arc = distance(point(1,2), point(1,1),lat(Sets{i}), long(Sets{i})) * 111 ; % km

    sel = range(arc < Max_Dist);

    tdata = arc(sort(sel));
    ydata = dVe(sort(sel))*1000;

    fun = @(x)sseval(x,tdata,ydata);

    x0 = rand(2,1);
    bestx = fminsearch(fun,x0)

    A = bestx(1);
    lambda = bestx(2);
    yfit = A*exp(lambda*tdata);
    close all
    figure(7)
    hold on
    grid on
    plot(tdata,ydata,'*');
    plot(tdata,yfit,'.r');
    xlabel('tdata')
    ylabel('Response Data and Curve')
    title('Data and Best Fitting Exponential Curve')
    legend('Data','Fitted Curve')
    hold off
    pause(0.5)
end