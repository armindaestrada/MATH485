close all
clear all      
t= 0;        
for graph_num = [0,1] %Zero for graph on right from paper, one for graph on left.
clear P1 and P2 and P3 and P4 and L and LG
switch graph_num
    case 0                  %Sets parameters to values for graph on right 
        N=500000;
        tau = 300;
        Ua = 5*10^-8;
        Ub=Ua;
        d=.1;
        z= 15;
        y=5;            %Gamma
        alpha=.5;
        ymin = 10^-3;
        ymax = 1;
        t_max = 10000;
    case 1               %Sets parameters to values for graph on left 
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
tP4 = [0:1:t_max];      %Initializes a vector from 0 to t_max

f = [1];                %Initializes first value of f.

for j = 2:length(tP4)+1 %Create values for rest of vector f.

   update = (1 - Ua/(2*tau)-alpha*Ua/tau);

   fj = f(j-1)*update;

   f = [f;fj];

end

L = [];
FG=[];

for i = 2:length(tP4)+1

   Li = 1 - f(i-1)*(1-Ua/2);

   L = [L; Li];
end

for t = 1:t_max
    for k=1:t-z
        first = ((2^z)-1)*Ub;
        second = -exp((-((2^z)-1)*Ub))*(2^z)*Ub*(t-(z+k));
        G(k) = (1-exp(first + second)); 
        if graph_num == 0                   %Saves time if this equation isn't needed.
            F(k) = exp((-(2^z-1)*Ua*(k-1)))*(2^z-1)*Ua;
            FG(k) = F(k)*G(k)*(1+d); 
        end
        LG(k) = L(k)*G(k)*(1+d);
    end
    if t>z
        inside = (-N/tau)*sum(FG);
        inside2 = (-N/tau)*sum(LG); 
    else
        inside = 0;
        inside2 = 0;
    end
    P1(t) = 1-exp((-N*Ub*t*(1+d)*(.5+alpha))/(tau)); 
    if graph_num
        P2(t) = ((N/(2*tau))*2^(2*z))*Ub*Ua*((t-z)^2)*(1+d);    %Accurate only for certain parameter values. Start x-axis at 10.
    else
        P2(t) = 1-exp(inside);              %Actual equation, takes longer    
    end
    P3(t) = (N/tau)*(2^(z+y-1))*(z+1)*Ua*Ub*(t-(z/2+y))*(1+d);
    P4(t) = 1-exp(inside2);
end
tvec = 1:t_max;
subplot(1,2,(-graph_num+2))
semilogy(tvec,P1,'blue')
hold on                                     %Plots all on a single graph.
semilogy(tvec,P2,'black')
semilogy(tvec,P3,'green')
semilogy(tvec,P4,'red')
legend('P1','P2','P3','P4')                 %Adds a legend
legend('Location','southeast')
axis([10 t_max ymin ymax+ymax*.1])          %Sets axis to one similar to that used in the paper.
xlabel('Time, \it{t}');
ylabel('Probability of Cancer Initiation')  %Labels axes 
end