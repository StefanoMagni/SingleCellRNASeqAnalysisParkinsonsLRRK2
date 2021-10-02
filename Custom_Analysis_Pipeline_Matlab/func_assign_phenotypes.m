
function all_phenotypes_assignment=func_assign_phenotypes(list_pheno,folder_pheno,list_of_genes,all_expr_log,...
    all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes,flag_comp_log,...
    flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
    flag_plot,th_prop_cells,num_dim_tsne,min_perc_exp,th_expr)

indexes_clusters=all_assignments{1,1};
all_phenotypes_assignment=zeros(length(list_pheno),length(indexes_clusters));

for i_list=1:length(list_pheno)
    
    all_genes_to_explore=write_struct_fromtxt([folder_pheno list_pheno{i_list}]);
    
    % Only choose the genes present in the list_if_genes
    ind_val=[];
    all_indexes_gene=[];
    for i_gene=1:size(all_genes_to_explore,1)
        all_list_gene=all_genes_to_explore{i_gene};
        ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            ind_val=cat(2,ind_val,i_gene);
            all_indexes_gene=cat(2,all_indexes_gene,ind_gene);
        end
    end
    all_genes_to_explore=all_genes_to_explore(ind_val,1);
    all_genes_to_explore_in=all_genes_to_explore;
    % We only choose those genes express in at least a minimum number of
    % cells
    aux_all_expr=all_expr{1,3};
    aux_all_expr=sum((aux_all_expr(all_indexes_gene,:)>=th_expr),2)/size(aux_all_expr,2);
    ind_val=find(aux_all_expr>min_perc_exp);
    all_genes_to_explore=all_genes_to_explore(ind_val,1);
    
    if ~isempty(all_genes_to_explore_in)
        all_pvalues_and_comp=cell(3,1);
        indexes_all_genes=zeros(1,size(all_genes_to_explore,1));
        table_pvalues=zeros(number_comb,size(all_genes_to_explore,1));
        
        if flag_comp_log
            data_group=all_expr_log{1,3};
        else
            data_group=all_expr_sca_cluster{1,3};
        end
        
        % Get the cumulative effect for the list of genes
        Both_cum_effect_all_genes=[];
        data_group_expr=all_expr{1,3};
        for i_gene=1:size(all_genes_to_explore_in,1)
            all_list_gene=all_genes_to_explore_in{i_gene};
            ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
            ind_gene=find(ind_gene==1);
            
            if ~isempty(ind_gene)
                Both_cum_effect_all_genes=cat(1,Both_cum_effect_all_genes,...
                    data_group_expr(ind_gene,:));
            end
        end
        
        if flag_comb_cumgenes==1
            Both_cum_effect_all_genes=sum(Both_cum_effect_all_genes);
        elseif flag_comb_cumgenes==2
            Both_cum_effect_all_genes=median(Both_cum_effect_all_genes);
        elseif flag_comb_cumgenes==3
            Both_cum_effect_all_genes=mean(Both_cum_effect_all_genes);
        end
        Both_cum_effect_all_genes_log=log10(Both_cum_effect_all_genes+1);
        if flag_log_cum
            Both_cum_effect_all_genes_st=Both_cum_effect_all_genes_log;
        else
            Both_cum_effect_all_genes_st=Both_cum_effect_all_genes;
        end
        
        %%%%%
        
        for i_gene=1:size(all_genes_to_explore,1)
            
            all_list_gene=all_genes_to_explore{i_gene};
            ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
            ind_gene=find(ind_gene==1);
            
            [table_pvalues_gene,table_comparison_and_ind]=evaluate_comb(ind_gene,...
                each_num_comp,indexes_clusters,data_group,number_comb,0,number_group_per_cluster);
            table_pvalues(:,i_gene)=table_pvalues_gene;
            
        end
        % For the cumulative expression
        [table_pvalues_cum_exp,table_comparison_and_ind]=evaluate_comb(1,...
            each_num_comp,indexes_clusters,Both_cum_effect_all_genes_st,...
            number_comb,1,number_group_per_cluster);
        
        table_pvalues_cum_exp(:,4)=table_pvalues_cum_exp(:,2)./table_pvalues_cum_exp(:,3);
        table_pvalues_cum_exp(:,5)=table_pvalues_cum_exp(:,3)./table_pvalues_cum_exp(:,2);
        table_pvalues_cum_exp(:,6)=(table_pvalues_cum_exp(:,4)>th_prop_cells)...
            +(table_pvalues_cum_exp(:,5)>th_prop_cells);
        
        all_pvalues_and_comp{1,1}=table_pvalues;
        all_pvalues_and_comp{2,1}=table_comparison_and_ind;
        all_pvalues_and_comp{3,1}=table_pvalues_cum_exp;
        
        cluster_more_expressed=zeros(1,length(indexes_clusters));
        % Now we see which clusters are mainly highly expressing
        for i_gene=1:size(table_pvalues,2)
            
            all_list_gene=all_genes_to_explore{i_gene};
            ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
            ind_gene=find(ind_gene==1);
            
            aux_table_pvalues=table_pvalues(:,i_gene);
            if sum(aux_table_pvalues<(thres_pval/size(table_pvalues,2)))
                
                [min_p,ind_min]=min(aux_table_pvalues);
                
                aux_indexes_assign=table_comparison_and_ind{ind_min,3};
                
                data_group_gene=data_group(ind_gene,:);
                if mean(data_group_gene(aux_indexes_assign==1))>mean(data_group_gene(aux_indexes_assign==2))
                    cluster_highexp=table_comparison_and_ind{ind_min,1};
                else
                    cluster_highexp=table_comparison_and_ind{ind_min,2};
                end
                for i_c=1:length(cluster_highexp)
                    cluster_more_expressed(indexes_clusters==cluster_highexp(i_c))=...
                        cluster_more_expressed(indexes_clusters==cluster_highexp(i_c))+1;
                end
                if flag_plot
                    p_val_plot=min_p;
                    data_group_expr=all_expr{1,3};
                    data_group_red=all_expr_sca_cluster{1,3};
                    data_group_expr_log=all_expr_log{1,3};
                    
                    figure;
                    %FigHandle = figure('Position', [100, 100, 1200, 900]);
                    
                    %colorsA=jet(max(data_group_expr(ind_gene,:))+1);
                    colorsA=jet(256);
                    %aux_exp=data_group_expr(ind_gene,:);
                    aux_exp=round(255*(data_group_expr(ind_gene,:))/(max(data_group_expr(ind_gene,:))));
                    scatter_colorCoded(data_group_red,colorsA,aux_exp,num_dim_tsne);
                    title(['List: ' list_pheno{i_list} ', gene ' all_list_gene ', Pval ' num2str(p_val_plot) ', Clusters ' num2str(cluster_highexp)])
                    colorbar
                    colormap('jet')
                    %                     if num_dim_tsne==2
                    %                         gscatter(data_group_red(:,1),data_group_red(:,2),indexes_comp);
                    %                         title(['Gene ' all_list_gene ' Pval ' p_val_plot ' Clusters ' num2str(cluster_highexp)])
                    %                     elseif num_dim_tsne==3
                    %                         gscatter3_custom(data_group_red,data_group_expr(ind_gene,:),h1);
                    %                         title(['Gene ' all_list_gene ' Pval ' p_val_plot ' Clusters ' num2str(cluster_highexp)])
                    %                     end
                end
            end
            % STILL THE PLOT FOR THE CUMULATIVE EXPRESSION
            
        end
        
        % The plot for cluster_more_expressed
        figure
        colorsA=jet(max(cluster_more_expressed)+1);
        scatter_colorCoded(data_group_red,colorsA,cluster_more_expressed,num_dim_tsne);
        max_val=max(cluster_more_expressed);
        ind_all_max=find(cluster_more_expressed==max_val);
        clust_max_exp=unique(indexes_clusters(ind_all_max));
        title(['List: ' list_pheno{i_list} ', clusters highly expressed more times: ' num2str(clust_max_exp')]);
        colorbar
        colormap('jet')
        legend(num2str(unique(cluster_more_expressed)))
        all_phenotypes_assignment(i_list,:)=cluster_more_expressed;
        
        % The plot for the cumulative expression
        figure
        [min_p,ind_min]=min(aux_table_pvalues);
        aux_indexes_assign=table_comparison_and_ind{ind_min,3};
        if mean(data_group_gene(aux_indexes_assign==1))>mean(data_group_gene(aux_indexes_assign==2))
            cluster_highexp=table_comparison_and_ind{ind_min,1};
        else
            cluster_highexp=table_comparison_and_ind{ind_min,2};
        end
        %colorsA=jet(max(Both_cum_effect_all_genes)+1);
        colorsA=jet(256);
        %aux_exp=data_group_expr(ind_gene,:);
        aux_exp=round(255*(Both_cum_effect_all_genes)/(max(Both_cum_effect_all_genes)));
        scatter_colorCoded(data_group_red,colorsA,aux_exp,num_dim_tsne);
        title(['List: ' list_pheno{i_list} ', cumulative expression, cluster more highly expressed ' num2str(cluster_highexp)])
        colorbar
        colormap('jet')
        
    else
        
        var_empty=0;
        
    end
    
end