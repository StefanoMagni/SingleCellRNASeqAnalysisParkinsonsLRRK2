

folder_files='./jonas/';
day_comp=14; % 0, 10, 14, or 42
day_comp_index=3;
name_groups={'A','B'}; % A is mutated and B control
% The columns are A and B, and the rows 0, 10, 14 or 42 days
numb_cells=[500,500;250,250;505,350;400,300];


%%%
files_to_compare={[name_groups{1,1} '_' num2str(day_comp) '.txt'],...
    [name_groups{1,2} '_' num2str(day_comp) '.txt']};



all_tables_comp={};
for i=1:length(files_to_compare)
    table=readtable([folder_files files_to_compare{1,i}],'Delimiter','tab','ReadRowNames',1,'ReadVariableNames',1);
    all_tables_comp{1,i}=table;    
end

% Get tables with only common genes
names_genes={};
for i=1:length(all_tables_comp)
    names_g=all_tables_comp{1,i}.Properties.RowNames;
    names_genes{1,i}=names_g;
end
[list_of_genes,ind_inters1,ind_inters2]=intersect(names_genes{1,1},names_genes{1,2});

aux=table2array(all_tables_comp{1,1});
all_arrays_comp_red{1,1}=aux(ind_inters1,:);
aux=table2array(all_tables_comp{1,2});
all_arrays_comp_red{1,2}=aux(ind_inters2,:);

% Ordered by the level of expression, and then select the number indicated
% in numb_cells
% We first emulate function get_top_cells_matrix(M,n), to select only the
% cells with higher expressions

A_expr=all_arrays_comp_red{1,1};
A_expr_max=get_top_cells_matrix(A_expr,numb_cells(day_comp_index,1));

B_expr=all_arrays_comp_red{1,2};
B_expr_max=get_top_cells_matrix(B_expr,numb_cells(day_comp_index,2));

% Also, reduce the tables by only considering those genes expressed in both
% groups
genes_expA=sum(A_expr_max,2);
genes_expB=sum(B_expr_max,2);

ind_genes_expA=find(genes_expA>0);
ind_genes_expB=find(genes_expB>0);
ind_genes_exp_tot=intersect(ind_genes_expA,ind_genes_expB);

A_expr=A_expr(ind_genes_exp_tot,:);
B_expr=B_expr(ind_genes_exp_tot,:);
A_expr_max=A_expr_max(ind_genes_exp_tot,:);
B_expr_max=B_expr_max(ind_genes_exp_tot,:);

Both_expr_max=cat(2,A_expr_max,B_expr_max);
labels=[zeros(1,size(A_expr_max,2)),ones(1,size(B_expr_max,2))];

list_of_genes=list_of_genes(ind_genes_exp_tot);

% Then, we perform a normalization to obtain for each cell the number
% of transcripts per 10.000 (in python P.append(10000*M[:,i]/sum(M[:,i]))),
% and also obtain the zscore, i.e., normalize to have mean 0 and variance 1
[A_expr_norm,A_expr_sca]=get_normalization(A_expr_max);
[B_expr_norm,B_expr_sca]=get_normalization(B_expr_max);
Both_expr_sca=cat(2,A_expr_sca,B_expr_sca);

% From here, we really start now with the experiments

%% Checking Suresh results
% We first perform PCA, keep some PCs, and then tSNE
% [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X). Rows of X 
% correspond to observations and columns to variables.

%% Check the level of sparsity of each gene
if 0
    threshold_exp=0.3;
    
    sparsityA=zeros(1,size(A_expr_max,1));
    sparsityB=zeros(1,size(B_expr_max,1));
    for i=1:size(A_expr_max,1)
        sparsityA(i)=sum(A_expr_max(i,:)>0)/size(A_expr_max,2);
        sparsityB(i)=sum(B_expr_max(i,:)>0)/size(B_expr_max,2);
    end
    
    ind_nonSpA=find(sparsityA>threshold_exp);
    ind_nonSpB=find(sparsityB>threshold_exp);
    list_of_genesA=list_of_genes(ind_nonSpA);
    list_of_genesB=list_of_genes(ind_nonSpB);

    ind_nonSp_genes=intersect(ind_nonSpA,ind_nonSpB);
    
    A_expr_max=A_expr_max(ind_nonSpA,:);
    B_expr_max=B_expr_max(ind_nonSpB,:);
    A_expr_sca=A_expr_sca(ind_nonSpA,:);
    B_expr_sca=B_expr_sca(ind_nonSpB,:);
    Both_expr_sca=Both_expr_sca(ind_nonSp_genes,:);
    list_of_genes=list_of_genes(ind_nonSp_genes);
end

%% Check the level of dispersion of each gene
% Dispersion = variance/mean expression
if 1
    threshold_disp=0.1;
    
    mean_aver_expA=mean(A_expr_max,2);
    mean_aver_expB=mean(B_expr_max,2);
    dispersionA=var(A_expr_max,[],2)./mean_aver_expA;
    dispersionB=var(B_expr_max,[],2)./mean_aver_expB;
    
    plot(mean_aver_expA,zscore(dispersionA),'ob');
    plot(mean_aver_expB,zscore(dispersionB),'or');    
end

%% Check the significance between groups

if 1
    
    p_ext=0.001;
    A_expr_max_log=log10(A_expr_max+1);
    B_expr_max_log=log10(B_expr_max+1);
    
    all_p_values=zeros(1,size(A_expr_max_log,1));
    for i=1:size(A_expr_max_log,1)
        
        % Equivalent to the Mann-Whitney test
        all_p_values(i)=ranksum(A_expr_max_log(i,:),B_expr_max_log(i,:));
        
    end
    ind_less=find(all_p_values<p_ext);
    
    list_of_genes=list_of_genes(ind_less);    
    Both_expr_sca=Both_expr_sca(ind_less,:);
    
    [~,score,~,~,~,~]=pca(Both_expr_sca');
    Both_expr_sca_trans=score';
    
    [Both_expr_sca_tsne]=tsne(Both_expr_sca_trans',[],num_dim_tsne);
    % [Both_expr_sca_tsne]=tsne(Both_expr_sca',[],num_dim_tsne);
    scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);

    model = train(labels',sparse(Both_expr_sca'),'-s 6 -c 0.001 -v 10');
    
end

%% Let's apply PCA
%%%% Mutant group
[~,score,eigval,~,~,~]=pca(A_expr_sca');
norm_eig_s=cumsum((eigval/norm(eigval)).^2);
ind_minA=min(find(norm_eig_s>0.999));
A_expr_sca_trans=score';

%%%% Control group
[~,score,eigval,~,~,~]=pca(B_expr_sca');
norm_eig_s=cumsum((eigval/norm(eigval)).^2);
ind_minB=min(find(norm_eig_s>0.999));
B_expr_sca_trans=score';

%%%% Both groups
[coeff_feat,score,eigval,tsquared,explained_pc,estimated_mean]=pca(Both_expr_sca');
norm_eig_s=cumsum((eigval/norm(eigval)).^2);
ind_minBoth=min(find(norm_eig_s>0.999));
Both_expr_sca_trans=score';

%% And now tsne
number_pcs=20;
num_dim_tsne=3;
gene_colorcode='SOX2'; %'SOX2', 'MAP2'

ind_gene=cellfun(@(s) (strcmp(s,gene_colorcode)), list_of_genes);
ind_gene=find(ind_gene==1);

colorsA=jet(max(A_expr_max(ind_gene,:))+1);
colorsB=jet(max(B_expr_max(ind_gene,:))+1);
colorsBoth=jet(max([A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)])+1);
colors_label=[1,0,0;0,0,1];
%

%%%% Mutant group
%[A_expr_sca_tsne]=tsne(A_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
[A_expr_sca_tsne]=tsne(A_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
%plot(A_expr_sca_tsne(:,1),A_expr_sca_tsne(:,2),'o')
%plot3(A_expr_sca_tsne(:,1),A_expr_sca_tsne(:,2),A_expr_sca_tsne(:,3),'o')
figure
scatter_colorCoded(A_expr_sca_tsne,colorsA,A_expr_max(ind_gene,:),num_dim_tsne)

%%%% Control group
%[B_expr_sca_tsne]=tsne(B_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
[B_expr_sca_tsne]=tsne(B_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
figure
scatter_colorCoded(B_expr_sca_tsne,colorsB,B_expr_max(ind_gene,:),num_dim_tsne)

%%%% Both group
%[Both_expr_sca_tsne]=tsne(Both_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
[Both_expr_sca_tsne]=tsne(Both_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
figure
scatter_colorCoded(Both_expr_sca_tsne,colorsBoth,...
    cat(2,A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)),num_dim_tsne);
% Color code by label
%scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);


%% Characterize which genes are more expressed in each group and cluster

model = train(labels',sparse(Both_expr_sca'),'-s 6 -c 1 -v 10');
model = svmtrain(labels',Both_expr_sca','-v 5');

%% Instead, let's separate the clusters by fitting a GMM

options = statset('Display','final');
gmfitted=fitgmdist(Both_expr_sca_tsne,4,'Replicates',50,'Options',options);

%% TEST



