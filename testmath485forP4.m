close all
clear all       %Ideally a loop or function could be used to condense this(You'll see what I mean), but lazy. 
t= 0;
t_max = 10000;
N = 10000; % number of stem cells about 5*10^5
tau = 10;
P_SC = 1/tau; % probability of SC dividing unevenly
z = 14; %number of times progenitor cell divides before dying
Ub = 10^-8; %probability of assymetric division gaining jak2 mutation
Ua=Ub;
alpha = 2; %ratio of symmetric vs. asymmetric
beta = 5; %additonal number of divisions
Z_Jak2 = z +beta; %number of divisions progenitor cells with jak2 mutations will undergo before dying
d = 0.05; %death rate of cell 
y=3; %Gamma
%%%
ymin = 10^-5;
ymax=1;
if 1
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
    control=1;
end
if 0
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
    control=0;
end


%for t=1:t_max
%    j=t-1;
%    f(t) = (1-(Ua/(2*tau))-alpha*Ua/tau)^j;
%end
%sum(f)/length(f)
%length(f)
%disp('done')

for t = 1:t_max
    for k=1:t-z
        first = ((2^z)-1)*Ub;
        second = exp((-((2^z)-1)*Ub))*(2^z)*Ub*(t-(z+k));
        G(k) = abs(1-exp(first + second)); %Always negative, since it's a probability, we'll take the absolute value
        if G(k)>1
            G(k);
            k;
        end
        %L(k) = (f(k)*Ua/2)+(1-f(k));
        L(k) = Ua/2;
        LG(k) = L(k)*G(k)*(1+d); 
    end
    if t>z
        inside2 = (-N/tau)*sum(LG); 

    else
        inside2 = 0;
    end
    P4(t) = 1-exp(inside2);
end
tvec = 1:t_max;
semilogy(tvec,P4,'red')
axis([10 t_max ymin ymax+ymax*.1])