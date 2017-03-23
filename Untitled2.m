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
    control=0
end
%%%

%may want to update death rate after mutation 
lineage_cell = z*N*P_SC; %number of progenitor cell at each time step
Total_Num_Prog = 2*exp((z+1) - 1 )*N*P_SC; %total number of progenitor cells 
prog_div = N*P_SC*(2*exp(z) -1); % number of progenitor cells that will divide in two
prog_death = N*P_SC*(2*exp(z)); % number of progenitor cells that will die
for t=1:t_max
    j=t-1;
    % base equation: f(j+1) = f(j)*(1-(Ua/(2*tau))-alpha*Ua/tau);
    f(t) = (1-(Ua/(2*tau))-alpha*Ua/tau)^j;
end
for k=1:t_max-z
    F(k) = exp((-(2^z-1)*Ua*(k-1)))*(2^z-1)*Ua;
    G(k) = 1-exp((-((2^z)-1)*Ub)+exp((-((2^z)-1)*Ub))*2^z*Ub*(t_max-z-k));
    S(k) = F(k)*G(k)*(1+d);
    L(k) = (f(k)*Ua/2)+(1-f(k));
    J(k) = L(k)*G(k)*(1+d);
    %J(k) = abs(J(k));
end
sum((J));

for t=1:t_max
    P1(t) = 1-exp((-N*Ub*t*(1+d)*(.5+alpha))/(tau)); %Use for Case Two (Left)
    if control
        P2(t) = 1-exp((-N/tau))*(sum(S(1:t-z))); %Use for Case One 
        (sum(S(1:t-z)))
    else
        P2(t) = ((N/(2*tau))*2^(2*z))*Ub*Ua*((t-z)^2)*(1+d); 
    end
    %if P2(t)>1     %real threshold is %((2^(z)-1)*Ua)*(t-1)<<1;
    %    P2(t) = 1-exp((-N/tau))*(sum(S(1:t-z)));
    %end
    
    P3(t) = (N/tau)*(2^(z+y-1))*(z+1)*Ua*Ub*(t-(z/2+y))*(1+d);
    %P4(t) = 1-exp((-N/tau)*sum...
    P4(t) = 1-exp((-N/tau)*(sum(J(1:t-z))));
    %    LsubK*GsubK*(1+d)) %sum from k=1 to t-z
    %    LsubK = expected nuber of progenitor 
end
tvec = 1:t_max;
subplot(1,2,2)
semilogy(tvec,P1,'blue')
hold on
tvec2 = 1:length(P2); 
semilogy(tvec2,P2,'black','LineWidth',2)
semilogy(tvec,P3,'green')
man = 0;

tvec3 = z:length(P2)+z-1; %Need to shift P4 to the right to match data.
semilogy(tvec2,P4)
legend('P1','P2','P3','P4')
legend('Location','southeast')
axis([10 t_max ymin ymax+ymax*.1])
%ts = 0:length(S)-1;
%plot(ts,(-N/tau).*S)
xlabel('Time, t');
ylabel('Probability of Cancer Initiation')
%title(['Investigates the importance of'; 'the trajectories for the most '; 'accurate parameter values of  ';'the hematopoietic system      '])
title(['Investigates the importance of the          ';
       'trajectories for the most accurate parameter';
       'values of the hematopoietic system          '])

clear all
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
if 0
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
if 1
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
%%%

%may want to update death rate after mutation 
lineage_cell = z*N*P_SC; %number of progenitor cell at each time step
Total_Num_Prog = 2*exp((z+1) - 1 )*N*P_SC; %total number of progenitor cells 
prog_div = N*P_SC*(2*exp(z) -1); % number of progenitor cells that will divide in two
prog_death = N*P_SC*(2*exp(z)); % number of progenitor cells that will die
for t=1:t_max
    j=t-1;
    % base equation: f(j+1) = f(j)*(1-(Ua/(2*tau))-alpha*Ua/tau);
    f(t) = (1-(Ua/(2*tau))-alpha*Ua/tau)^j;
end
%sum(f) == sum(abs(f)); %True, so all of f is non-negative

for k=1:t_max-z
    F(k) = exp((-(2^z-1)*Ua*(k-1)))*(2^z-1)*Ua;
    G(k) = 1-exp((-((2^z)-1)*Ub)+exp((-((2^z)-1)*Ub))*2^z*Ub*(t_max-z-k));
    S(k) = F(k)*G(k)*(1+d);
    L(k) = (f(k)*Ua/2)+(1-f(k));
    J(k) = L(k)*G(k)*(1+d);
    J(k) = abs(J(k));
end
sum((J));

for t=1:t_max
    P1(t) = 1-exp((-N*Ub*t*(1+d)*(.5+alpha))/(tau)); %Use for Case Two
    if control
        P2(t) = 1-exp((-N/tau))*(sum(J(1:t-z))); %Use for Case One
    else
        P2(t) = ((N/(2*tau))*2^(2*z))*Ub*Ua*((t-z)^2)*(1+d); 
    end
    %if P2(t)>1     %real threshold is %((2^(z)-1)*Ua)*(t-1)<<1;
    %    P2(t) = 1-exp((-N/tau))*(sum(S(1:t-z)));
    %end
    
    P3(t) = (N/tau)*(2^(z+y-1))*(z+1)*Ua*Ub*(t-(z/2+y))*(1+d);
    %P4(t) = 1-exp((-N/tau)*sum...
    P4(t) = 1-exp((-N/tau)*(sum(J(1:t-z))));
    %    LsubK*GsubK*(1+d)) %sum from k=1 to t-z
    %    LsubK = expected nuber of progenitor 
end
tvec = 1:t_max;
subplot(1,2,1)
semilogy(tvec,P1,'blue')
hold on
tvec2 = 1:length(P2); 
semilogy(tvec2,P2,'black')
semilogy(tvec,P3,'green')
semilogy(tvec2,P4)     %Totes the wrong equation
legend('P1','P2','P3','P4')
legend('Location','southeast')
axis([10 t_max ymin ymax])
%ts = 0:length(S)-1;
%plot(ts,(-N/tau).*S)
xlabel('Time, t');
ylabel('Probability of Cancer Initiation')
title('Left')
