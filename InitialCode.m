t= 0;
t_max = 10000;
N = 100; % number of stem cells about 5*10^5
tau = 10;
P_SC = 1/tau; % probability of SC dividing unevenly
z = 14; %number of times progenitor cell divides before dying
Ub = 10^-8; %probability of assymetric division gaining jak2 mutation
alpha = 2; %ratio of symmetric vs. asymmetric
beta = 5; %additonal number of divisions
Z_Jak2 = z +beta; %number of divisions progenitor cells with jak2 mutations will undergo before dying
d = 0.05; %death rate of cell 
%may want to update death rate after mutation 
lineage_cell = z*N*P_SC; %number of progenitor cell at each time step
Total_Num_Prog = 2*exp((z+1) - 1 )*N*P_SC; %total number of progenitor cells 
prog_div = N*P_SC*(2*exp(z) -1); % number of progenitor cells that will divide in two
prog_death = N*P_SC*(2*exp(z)); % number of progenitor cells that will die

%Prog_cells = [];
%if
% newinfo = ;
%prog_cells = [prog_cells)

%end