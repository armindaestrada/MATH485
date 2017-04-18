clear all
close all
   for   y= [0,8,20,-8,7]    ;        %Gamma
   clear H and P3
    t= 0;   
graph_num = 0;
        N=500000;
        tau = 300;
        Ua = 5*10^-8;
        Ub=Ua;
        d=.1;
        z= 15;
        alpha=.5;
        ymin = 10^-3;
        ymax = 1;
        t_max = 10000;

for t = 1:t_max
for i = 1:z
    H(i) = (2^(i-1))*Ub*(2^(z+y-i)-1)*Ua*(t-(z+y-i));
end
P3(t) = 1-exp((-N/tau)*sum(H)*(1+d));

end


tvec = 1:t_max;

semilogy(tvec,P3,'green')

legend('Location','southeast')
          %Sets axis to one similar to that used in the paper.
xlabel('Time, \it{t}');
ylabel('Probability of Cancer Initiation')
title(y) 
saveas(gcf,strcat('math485_gamma_equals_',int2str(y),'.png'))
   end
