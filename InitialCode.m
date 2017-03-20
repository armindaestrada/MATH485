t= 0;
t_max = 10000;
N = 100; % number of stem cells about 5*10^5
tau = 10;
P_SC = 1/tau; % probability of SC dividing unevenly
z = 14; %number of times progenitor cell divides before dying
Ub = 10^-8; %probability of assymetric division gaining jak2 mutation
Ua = 2*10^-6; % symetric division without altering cell's phenotype
alpha = 2; %ratio of symmetric vs. asymmetric
gamma = 5; %additonal number of divisions
Z_Jak2 = z +gamma; %number of divisions progenitor cells with jak2 mutations will undergo before dying
d = 0.05; %death rate of cell 
%may want to update death rate after mutation 
lineage_cell = z*N*P_SC; %number of progenitor cell at each time step
Total_Num_Prog = 2*exp((z+1) - 1 )*N*P_SC; %total number of progenitor cells 
prog_div = N*P_SC*(2*exp(z) -1); % number of progenitor cells that will divide in two
prog_death = N*P_SC*(2*exp(z)); % number of progenitor cells that will die
fk = 1; %


%Prog_cells = [];
%if
% newinfo = ;
%prog_cells = [prog_cells)

%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Different Scenarios%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Scenario (i). A Stem Cell Acquires JAK2V617F.%%%%%
P1 = 0;  %probability of Scenerio 1
dP1 = (1-P1)*(N/tau)*(1+d)*(1/2+alpha)*Ub;
%N/tau: number of stem cell divides
%dN/tau: expected # of births stem cells undergo per time step to compensate for cell death
% (1/2 + alpha)*ub: probablilty that mutation arises
% 1-P1: probability Jak2 mutation has not occured. 




%%%%%Scenario (ii). A Progenitor Acquires a Mutation Conferring
%%%%%Self-Renewal Followed by JAK2V617F%%%
% p2 = 1- exp(-(N/tau)*sum(Fk*Gk*(1-d),1,t-z)); 
%  Fk = exp(-(2^(z - 1)*ua*(K-1)))*((2^z)-1)*ua;   %expected number of cells in which a mutation conferring self-renewal emerged at time k,
%  Gk = 1- exp(-((2^z -1)*ub + exp(-((2^z) -1)*ub)*(2^z)*ub*(t-z-k))); %probability that the JAK2V617F mutation emerges in the clone carrying the first mutation between times k and t.
 
 
 p2 = (N/2)*P_SC*2^(2*z)*Ua*Ub*((t-z)^2)*(1+d);
 
 


%%%%%Scenario (iii). A Progenitor Acquires JAK2V617F Followed by a Mutation
%%%%%Conferring Self-Renewal.%%%%

p3 = N*P_SC*(2^(z + gamma - 1))*(z+1)*Ua*Ub*((t - (z/2 + gamma))*(1+d));


%Hi = 2^(i - 1)*Ub*(

%%%%%Scenario (iv). A Stem Cell Acquires a Mutation Conferring Self-Renewal
%%%%%to a Progenitor, Followed by the JAK2V617F Mutation Arising in a
%%%%%Progenitor.%%%%%%%%

fk = (1 - Ua/2 * (P_SC) -(alpha)*Ua*(P_SC));
Lk = (fk*Ua/2 + (1-fk)); %expected number of progenitors at the most undifferentiated stage carrying the self-renewal mutation
%p4 = 1 - exp(-(N+P_SC)*sum(Lk*Gk, 1, t-z)*(1-d));



















