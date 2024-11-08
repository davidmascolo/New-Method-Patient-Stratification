% Iterated Generalized Louvain Algorithm

% Import Data
data1 = BREASTGeneExpression2;
data2 = BREASTMirnaExpression2;
data3 = BREASTMethyExpression2;
id    = BREASTSurvival;

%data1 = filteredgene;
%data2 = filteredmirna;
%data3 = filteredmethy;
%data4 = filteredmutations;
%id    = filteredidsurvival;

% Table to array
gene = table2array(data1)';
mirna = table2array(data2)';
metil = table2array(data3)';
%mutations = table2array(data4)';

% Set parameters
k = 20;
alpha = 0.5;
t = 20;
num_patients = 105;
num_layers = 3;

% Standardization 
data_1 = Standard_Normalization(gene);  
data_2 = Standard_Normalization(mirna);   
data_3 = Standard_Normalization(metil);  
%data_4 = Standard_Normalization(mutations);  

% Compute distances
Dist1 = dist2(data_1,data_1);  
Dist2 = dist2(data_2,data_2);  
Dist3 = dist2(data_3,data_3); 
%Dist4 = dist2(data_4,data_4); 

% Similarity Networks
w1 = affinityMatrix(Dist1,k,alpha); 
w2 = affinityMatrix(Dist2,k,alpha); 
w3 = affinityMatrix(Dist3,k,alpha);
%w4 = affinityMatrix(Dist4,k,alpha);

% Fused Network
%A = {w1,w2,w3,w4};
A = {w1,w2,w3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the algorithm one time
gamma = 0.7;
[B,twom]=multicat(A,gamma,1);
[S,Q,n_it] = iterated_genlouvain(B);
S=reshape(S,num_patients,num_layers);
L = S(:,1)';

tabulate(L)
writematrix(L,"C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/BREAST/Res/L.txt");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the algorithm 1000 times with for cycle
A = {w1,w2,w3};

rng(1)                           % genlouvain  depends by starting nodes, so set a seed
gamma = 0.7;
omega = 1;
iter = 1000;                     % iterations number
L = zeros(iter,num_patients);    % matrix iter*patients, where each row is the classification of an interation
equal = zeros(iter,1);  

% Iterated Generalized Louvain
for i =1:iter
[B,twom]=multicat(A,gamma,omega);
[S,Q,n_it] = iterated_genlouvain(B,10000,0,0,1);
S=reshape(S,num_patients,num_layers);
equal(i) = isequal(S(:,1),S(:,2))*isequal(S(:,1),S(:,3))==1;
L(i,:) = S(:,1);
end

% Check
tabulate(equal)

% quanti clusters ci sono nella matrice ricavata L?
len = zeros(iter,1);
for i = 1:iter 
    len(i) = size(tabulate(L(i,:)),1);
end

length(find(len==2))
length(find(len==3))
length(find(len==4))

writematrix(L,"C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/BREAST/Res/L.txt");