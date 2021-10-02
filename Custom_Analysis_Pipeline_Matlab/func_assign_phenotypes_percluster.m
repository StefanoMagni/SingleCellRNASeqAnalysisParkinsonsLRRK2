

function [prob_cluster_expressed,all_genes_to_check,index_group_genes,index_group_genes_cluster,...
    pvalues_all_clusters,pval_and_cum_cluster_expressed]=...
    func_assign_phenotypes_percluster(list_pheno,folder_pheno,list_of_genes,all_expr_log,...
    all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes,flag_comp_log,...
    flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
    flag_plot,th_prop_cells,num_dim_tsne,num_perm,th_expr,flag_per_gene)

indexes_clusters=all_assignments{1,1};
all_phenotypes_assignment=zeros(length(list_pheno),length(indexes_clusters));

all_genes_to_check={};
all_indexes_gene=[];
index_group_genes=[];
genes_per_group=[];
for i_list=1:length(list_pheno)
    
    all_genes_to_explore=write_struct_fromtxt([folder_pheno list_pheno{i_list}]);    
    % Only choose the genes present in the list_if_genes
    ind_val=[];
    for i_gene=1:size(all_genes_to_explore,1)
        all_list_gene=all_genes_to_explore{i_gene};
        ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            ind_val=cat(2,ind_val,i_gene);
            all_indexes_gene=cat(2,all_indexes_gene,ind_gene);
            index_group_genes=cat(2,index_group_genes,i_list);
        end
    end
    all_genes_to_explore=all_genes_to_explore(ind_val,1);
    all_genes_to_check=cat(1,all_genes_to_check,all_genes_to_explore);    
    genes_per_group=[genes_per_group,length(all_genes_to_explore)];
end

thres_pval_corr=thres_pval/length(all_genes_to_check);

data_group_log=all_expr_log{1,3};
data_group_log_red=data_group_log(all_indexes_gene,:);

data_group=all_expr{1,3};
data_group_red=data_group(all_indexes_gene,:);
prob_cluster_expressed=zeros(length(list_pheno),max(indexes_clusters));

genes_highexpr_per_clust=cell(1,max(indexes_clusters));
pvalues_all_clusters=zeros(size(data_group_log_red,1),max(indexes_clusters));

pval_and_cum_cluster_expressed=zeros(length(list_pheno),max(indexes_clusters),3);

for i_cluster=1:max(indexes_clusters)
    i_cluster
    ind_cell_cluster=find(indexes_clusters==i_cluster);
    ind_cell_rest=setdiff(1:length(indexes_clusters),ind_cell_cluster);
    
    p_val_gene=ones(size(data_group_log_red,1),1);
    if flag_per_gene
        for i_gene=1:length(all_genes_to_check)
            
            aux_pval=permtest(data_group_log_red(i_gene,ind_cell_cluster),data_group_log_red(i_gene,ind_cell_rest),...
                num_perm,'conservative');
            % Also check if the gene is more highly express in that cluster
            if mean(data_group_log_red(i_gene,ind_cell_cluster))>mean(data_group_log_red(i_gene,ind_cell_rest))
                p_val_gene(i_gene,1)=aux_pval;
            else
                p_val_gene(i_gene,1)=1;
            end
            
        end        
    end
    
    pvalues_all_clusters(:,i_cluster)=p_val_gene;
    ind_valp=find(p_val_gene<=thres_pval_corr);
    genes_highexpr_per_clust{1,i_cluster}=ind_valp;
    index_group_genes_cluster=index_group_genes(ind_valp);
    number_genes_pergroup_cl=histc(index_group_genes_cluster,(unique(index_group_genes_cluster)));
    prob_norm=number_genes_pergroup_cl./genes_per_group(unique(index_group_genes_cluster));
    prob_norm=prob_norm/sum(prob_norm);
    prob_cluster_expressed(unique(index_group_genes_cluster),i_cluster)=prob_norm;
    
    % And now, for the cumulative expression
    for i_list=1:length(list_pheno)
       
        ind_list_ph=find(index_group_genes==i_list);
        all_genes_to_check_aux=all_genes_to_check(ind_list_ph);
        
        if flag_comb_cumgenes==1 % 1 for sum, 2 for median, 3 for mean
            both_cum_expression=sum(data_group_red(ind_list_ph,:));
        elseif flag_comb_cumgenes==2
            both_cum_expression=median(data_group_red(ind_list_ph,:));
        elseif flag_comb_cumgenes==3
            both_cum_expression=mean(data_group_red(ind_list_ph,:));
        end
        if flag_log_cum
            both_cum_expression=log10(both_cum_expression+1);
        end
        aux_pval=permtest(both_cum_expression(ind_cell_cluster),both_cum_expression(ind_cell_rest),...
            num_perm,'conservative');   
        
        pval_and_cum_cluster_expressed(i_list,i_cluster,1)=aux_pval;
        pval_and_cum_cluster_expressed(i_list,i_cluster,2)=mean(both_cum_expression(ind_cell_cluster));
        pval_and_cum_cluster_expressed(i_list,i_cluster,3)=mean(both_cum_expression(ind_cell_rest));
        
    end
    
end
