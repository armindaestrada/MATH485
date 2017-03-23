clear all
close all
graph_num = 1;
switch graph_num
    case 0
        N=500000; %100; 
        tau = 300;% 5;
        Ua = 5*10^-8;% 2*10^-6;
        Ub=Ua;
        d=.1;
        z= 15;%9;
        y=5;
        alpha=.5;%1;
        ymin = 10^-3;
        ymax = 1;
        t_max = 10000;
    case 1
        N=100; 
        t_max = 100;
        tau = 5;
        Ua = 2*10^-6;
        Ub=Ua;
        d=.1;
        z= 9;
        y=3;
        alpha=1;
        ymin = 10^-5;
        ymax = 10^-1;
end

for t=1:t_max
    j=t-1;
    f(t) = (1-(Ua/(2*tau))-alpha*Ua/tau)^j;
end

disp('start')
for t = 1:t_max
    if 1 %graph_num == 1
        for k=1:t-z
            first = ((2^z)-1)*Ub;
            second = exp((-((2^z)-1)*Ub))*(2^z)*Ub*(t-(z+k));
            G(k) = abs(1-exp(first + second)); %Always negative, since it's a probability, we'll take the absolute value
            F(k) = exp((-(2^z-1)*Ua*(k-1)))*(2^z-1)*Ua; %Always Greater than zero
            FG(k) = F(k)*G(k)*(1+d); %Absolute value always less than one
            L(k) = (f(k)*Ua/2)+(1-f(k));
            LG(k) = L(k)*G(k)*(1+d); 
        end
        if t>z
            inside2 = (-N/tau)*sum(LG); 
        else
            inside2 = 0;
        end
    end
    P4(t) = 1-exp(inside2);
    P1(t) = 1-exp((-N*Ub*t*(1+d)*(.5+alpha))/(tau)); %Use for Case Two
    P2(t) = ((N/(2*tau))*2^(2*z))*Ub*Ua*((t-z)^2)*(1+d); 
    P3(t) = (N/tau)*(2^(z+y-1))*(z+1)*Ua*Ub*(t-(z/2+y))*(1+d);

end
tvec = 1:t_max;
semilogy(tvec,P1,'blue')
hold on
tvec2 = 1:length(P2); 
semilogy(tvec2,P2,'black')
semilogy(tvec,P3,'green')
semilogy(tvec,P4,'red')
legend('P1','P2','P3','P4')
legend('Location','southeast')
axis([10 t_max ymin ymax+ymax*.1])
xlabel('Time, t');
ylabel('Probability of Cancer Initiation')
title('Left')
disp('Done')