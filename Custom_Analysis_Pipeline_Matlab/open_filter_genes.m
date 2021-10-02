
function [all_expr,labels,list_of_genes]=open_filter_genes(folder_files,...
    day_comp,day_comp_index,name_groups,numb_cells)


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

A_expr_max_log=log10(A_expr_max+1);
B_expr_max_log=log10(B_expr_max+1);
Both_expr_max_log=log10(Both_expr_max+1);

all_expr.max={A_expr_max,B_expr_max,Both_expr_max};
all_expr.norm={A_expr_norm,B_expr_norm};
all_expr.sca={A_expr_sca,B_expr_sca,Both_expr_sca};
all_expr.log={A_expr_max_log,B_expr_max_log,Both_expr_max_log};


