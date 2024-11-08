% Iterated Generalized Louvain ALgorithm

% Import Data
data1 = KIDNEYGeneExpression2;
data2 = KIDNEYMirnaExpression2;
data3 = KIDNEYMethyExpression2;
id    = KIDNEYSurvival;

%data1 = filteredgene;
%data2 = filteredmirna;
%data3 = filteredmethy;
%data4 = filteredmutations;
%id    = filteredidsurvival;

gene  = table2array(data1)';
mirna = table2array(data2)';
metil = table2array(data3)';
%mutations = table2array(data4);

% Set parameters
k     = 20;
alpha = 0.5;
t     = 20;
num_patients = 122; 
num_layers   = 3;

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

% Fused Network with Similarity Network Fusion
%W = SNF({w1,w2,w3,w4}, k, t); 
W = SNF({w1,w2,w3}, k, t); 

% For 2 clusters
C = 2;
lab = SpectralClustering(W,C);
tabulate(lab)
writematrix(lab,"C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/GLIO/Res/labg00.txt");

% For 3 clusters
C = 3;
lab = SpectralClustering(W,C);
tabulate(lab)
writematrix(lab,"C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/GLIO/Res/labg01.txt");


% Generalized Louvain

% Compute euclidean distance
Dist1 = dist2(data_1,data_1);  
Dist2 = dist2(data_2,data_2);  
Dist3 = dist2(data_3,data_3);
%Dist4 = dist2(data_4,data_4);


% Pearson
%Dist1 = 1-corrcoef(data_1');  
%Dist2 = 1-corrcoef(data_2');  
%Dist3 = 1-corrcoef(data_3'); 

% Spearman
%Dist1 = 1-corr(data_1',data_1','Type','Spearman');  
%Dist2 = 1-corr(data_2',data_2','Type','Spearman');  
%Dist3 = 1-corr(data_3',data_3','Type','Spearman'); 
%Dist4 = 1-corr(data_4',data_4','Type','Spearman'); 


% Similarity Networks
w1 = affinityMatrix(Dist1,k,alpha); 
w2 = affinityMatrix(Dist2,k,alpha); 
w3 = affinityMatrix(Dist3,k,alpha);
%w4 = affinityMatrix(Dist4,k,alpha);

% Set Parameters
A = {w1,w2,w3};
rng(0)                            % genlouvain depends by starting nodes, so we set a seed
gamma = 0.82;
omega = 1;
iter = 1000;                      % iterations number
L = zeros(iter,num_patients);     % matrix iter*patients, where each row is the classification of an interation
equal = zeros(iter,1);            % vector of 0,1 where if the partitions' compositions are equal at each iteration

% Iterated algorithm generalized louvain
for i = 1:iter
[B,twom] = multicat(A,gamma,omega);
[S,Q] = genlouvain(B,10000,1,0,1);
equal(i) = isequal(S(1:num_patients),S(num_patients+1:num_patients*2))*isequal(S(1:num_patients),S(num_patients*2+1:num_patients*3))==1;
L(i,:) = S(1:num_patients);
end

tabulate(equal) % equal to 1 means that all the lables identified are equal in each layer
writematrix(L,"C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Res/L.txt");
% The partition matrix is used as input for the K-Means algorithm in R.














% --------------------------------------------------------------------------------------------------------------
% PRIMA Classificazione più frequente 
[uA,~,uIdx] = unique(L,'rows');
modeIdx = mode(uIdx);
modeRow = uA(modeIdx,:);            % the first output argument
whereIdx = find(uIdx==modeIdx);     % the second output argument
length(whereIdx)/iter*100;          % percentuale iterazioni PARTIZIONE 1


lab = L(whereIdx(1),:);
lab = lab';
tabulate(lab)
writematrix(lab,'C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/labG1.txt');


% SECONDA Classificazione più frequente 
L2 = L;
L2(whereIdx,:) = [];
[uA2,~,uIdx2] = unique(L2,'rows');
modeIdx2 = mode(uIdx2);
modeRow2 = uA2(modeIdx2,:)          % the first output argument
whereIdx2 = find(uIdx2==modeIdx2)   % the second output argument
length(whereIdx2)/iter*100          % percentuale iterazioni PARTIZIONE 2



lab2 = L2(whereIdx2(1),:);
lab2 = lab2';
tabulate(lab2)
writematrix(lab2,'C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/labG2.txt');


% TERZA Classificazione più frequente 
L3 = L2;
L3(whereIdx2,:) = [];
[uA3,~,uIdx3] = unique(L3,'rows');
modeIdx3 = mode(uIdx3);
modeRow3 = uA3(modeIdx3,:)          % the first output argument
whereIdx3 = find(uIdx3==modeIdx3)   % the second output argument
length(whereIdx3)/iter*100          % percentuale iterazioni PARTIZIONE 3

lab3 = L3(whereIdx3(1),:);
lab3 = lab3';
tabulate(lab3)
writematrix(lab3,'C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/labG3.txt');



% QUARTA Classificazione più frequente 
L4 = L3;
L4(whereIdx3,:) = [];
[uA4,~,uIdx4] = unique(L4,'rows');
modeIdx4 = mode(uIdx4);
modeRow4 = uA4(modeIdx4,:)          % the first output argument
whereIdx4 = find(uIdx4==modeIdx4)   % the second output argument
length(whereIdx4)/iter*100          % percentuale iterazioni PARTIZIONE 4


lab4 = L4(whereIdx4(1),:);
lab4 = lab4';
tabulate(lab4)
writematrix(lab4,'C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/labG4.txt');



% QUINTA Classificazione più frequente 
L5 = L4;
L5(whereIdx4,:) = [];
[uA5,~,uIdx5] = unique(L5,'rows');
modeIdx5 = mode(uIdx5);
modeRow5 = uA5(modeIdx5,:)          % the first output argument
whereIdx5 = find(uIdx5==modeIdx5)   % the second output argument
length(whereIdx5)/iter*100          % percentuale iterazioni PARTIZIONE 5



lab5 = L5(whereIdx5(1),:);
lab5 = lab5';
tabulate(lab5)
writematrix(lab5,'C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/labG5.txt');