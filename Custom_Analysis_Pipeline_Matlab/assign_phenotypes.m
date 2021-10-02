
% Assign phenotypes to clusters
folder_files='./jonas/';
day_comp=[0,10,14,42]; % [0,10,14,42]
%day_comp_index=1:3; % 1:4
day_comp_index=3;
day_comp=day_comp(day_comp_index);
name_groups={'A','B'}; % A is mutated and B control
% The columns are A and B, and the rows 0, 10, 14 or 42 days
numb_cells=[500,500;250,250;505,350;400,300];
number_pcs=40;

num_dim_tsne=3;

folder_res='./Images_results/Exp1/';
max_clusters=6;

file_genes='./Drop_Seq/list_CellCycle.txt'; % list_stem / list_CellCycle / list_total
all_genes_to_explore=write_struct_fromtxt(file_genes);
empty_cols=sum(not(cellfun(@isempty,all_genes_to_explore)));
all_genes_to_explore=all_genes_to_explore(:,find(empty_cols>0));
thres_pval=0.01;
max_num_genes=10;
th_prop_cells=10;

folder_pheno='./Drop_Seq/';
list_pheno={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt'};
flag_comb_cumgenes=1; % 1 for sum, 2 for median, 3 for mean

% Flags to do stuff
flag_plot=1;
flag_save=0;
if flag_save
    flag_plot=1;
end
flag_useTSNE=1;
flag_comp_log=1;
flag_min_pval=1;

%%%%%%%%%%%%

% Open and filter the genes
all_days_expr=cell(1,length(day_comp_index));
all_days_labels=cell(1,length(day_comp_index));
all_days_genes=cell(1,length(day_comp_index));
for i_f=1:length(day_comp_index)
    i_f
    [all_expr_norm,labels,list_of_genes]=open_filter_genes(folder_files,...
        day_comp(i_f),i_f,name_groups,numb_cells);
    
    %% Let's apply PCA
    %%%% Mutant group
    [~,score,eigval,~,~,~]=pca(all_expr_norm.sca{1,1}');
    norm_eig_s=cumsum((eigval/norm(eigval)).^2);
    ind_minA=min(find(norm_eig_s>0.999));
    A_expr_sca_trans=score';
    
    %%%% Control group
    [~,score,eigval,~,~,~]=pca(all_expr_norm.sca{1,2}');
    norm_eig_s=cumsum((eigval/norm(eigval)).^2);
    ind_minB=min(find(norm_eig_s>0.999));
    B_expr_sca_trans=score';
    
    %%%% Both groups
    [~,score,eigval,~,~,~]=pca(all_expr_norm.sca{1,3}');
    norm_eig_s=cumsum((eigval/norm(eigval)).^2);
    ind_minBoth=min(find(norm_eig_s>0.999));
    Both_expr_sca_trans=score';
    
    %% And now tsne
    colors_label=[1,0,0;0,0,1];
    
    %%%% Mutant group
    %[A_expr_sca_tsne]=tsne(A_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
    [A_expr_sca_tsne]=tsne(A_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
    %plot(A_expr_sca_tsne(:,1),A_expr_sca_tsne(:,2),'o')
    %plot3(A_expr_sca_tsne(:,1),A_expr_sca_tsne(:,2),A_expr_sca_tsne(:,3),'o')
    
    %%%% Control group
    %[B_expr_sca_tsne]=tsne(B_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
    [B_expr_sca_tsne]=tsne(B_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
    
    %%%% Both group
    %[Both_expr_sca_tsne]=tsne(Both_expr_sca_trans(1:number_pcs,:)',[],num_dim_tsne);
    [Both_expr_sca_tsne]=tsne(Both_expr_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
    
    %%%%%%%%
    if flag_useTSNE
        all_expr_sca_cluster{1,1}=A_expr_sca_tsne;
        all_expr_sca_cluster{1,2}=B_expr_sca_tsne;
        all_expr_sca_cluster{1,3}=Both_expr_sca_tsne;
    else
        all_expr_sca_cluster{1,1}=A_expr_sca_trans(1:number_pcs,:)';
        all_expr_sca_cluster{1,2}=B_expr_sca_trans(1:number_pcs,:)';
        all_expr_sca_cluster{1,3}=Both_expr_sca_trans(1:number_pcs,:)';
    end
    
    all_expr_log{1,1}=log10(all_expr_norm.max{1,1}+1);
    all_expr_log{1,2}=log10(all_expr_norm.max{1,2}+1);
    all_expr_log{1,3}=log10(all_expr_norm.max{1,3}+1);
    all_expr{1,1}=all_expr_norm.max{1,1};
    all_expr{1,2}=all_expr_norm.max{1,2};
    all_expr{1,3}=all_expr_norm.max{1,3};
    
    %% Now, fit a GMM to the distribution obtained with tSNE
    
    each_num_comp=zeros(2,1);
    all_structures_gm=cell(max_clusters,1);
    i_g=3;
    
    data_group=all_expr_sca_cluster{i_g};
    AIC=zeros(1,max_clusters);
    for k=1:max_clusters
        all_structures_gm{k,1}= fitgmdist(data_group,k,'Replicates',50);
        AIC(k)= all_structures_gm{k,1}.AIC;
    end
    [minAIC,numComponents]=min(AIC);
    each_num_comp(:,1)=[minAIC;numComponents];
    
    
    %% Get the assignment to each cluster
    
    all_assignments=cell(1,1);
    
    data_group=all_expr_sca_cluster{i_g};
    obj=all_structures_gm{each_num_comp(2,1),1};
    idx=cluster(obj,data_group);
    all_assignments{1,1}=idx;
    
    
    %% Plot the assignments
    
    if flag_plot
        
        data_group=all_expr_sca_cluster{i_g};
        indexes_clusters=all_assignments{1,1};
        obj=all_structures_gm{each_num_comp(2,1),1};
        h1=figure;
        if num_dim_tsne==2
            h1 = gscatter(data_group(:,1),data_group(:,2),indexes_clusters);
        elseif num_dim_tsne==3
            gscatter3_custom(data_group,indexes_clusters,h1);
        end
        leg_cluster={};
        for i_leg=1:numComponents
            leg_cluster{1,i_leg}=['Cluster ' num2str(i_leg)];
        end
        legend(leg_cluster)
        
        plot_gaussians(data_group,obj,indexes_clusters,h1,num_dim_tsne);
        if(i_g==1); group_val='mutant'; elseif(i_g==2); group_val='control'; elseif(i_g==3); group_val='both'; end
        title(['Cluster with Gaussian to 99% for ' group_val])
        
        if flag_save
            print([folder_res 'Clusters_Day' num2str(day_comp) '_' group_val '_k' num2str(each_num_comp(1,1)) '_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
        end
        
    end
    
    %% Obtain some statistics of how many of each group are there in each cluster
    
    percent_group_per_cluster=zeros(2,numComponents);
    number_group_per_cluster=zeros(2,numComponents);
    for i=1:max(indexes_clusters)
        aux_ind=find(indexes_clusters==i);
        percent=sum(labels(aux_ind)==0)/length(aux_ind);
        number_group_per_cluster(:,i)=[sum(labels(aux_ind)==0);sum(labels(aux_ind)==1)];
        percent_group_per_cluster(:,i)=[percent;1-percent]*100;
    end
    
    
    number_comb=0;
    for i=1:floor(each_num_comp(2,1)/2)
        if i==each_num_comp(2,1)/2
            number_comb=number_comb+(nchoosek(each_num_comp(2,1),i)/2);
        else
            number_comb=number_comb+nchoosek(each_num_comp(2,1),i);
        end
    end
    
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
            end
        end
        all_genes_to_explore=all_genes_to_explore(ind_val,1);
        
        all_pvalues_and_comp=cell(3,1);
        indexes_all_genes=zeros(1,size(all_genes_to_explore,1));
        table_pvalues=zeros(number_comb,size(all_genes_to_explore,1));
        
        if flag_comp_log
            data_group=all_expr_log{1,3};
        else
            data_group=all_expr_sca_cluster{1,3};
        end
        indexes_clusters=all_assignments{1,1};
        
        % Get the cumulative effect for the list of genes
        Both_cum_effect_all_genes=[];
        data_group_expr=all_expr{1,3};
        for i_gene=1:size(all_genes_to_explore,1)
            all_list_gene=all_genes_to_explore{i_gene};
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
            each_num_comp,indexes_clusters,Both_cum_effect_all_genes_log,...
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
                    title(['Gene ' all_list_gene ',Pval ' num2str(p_val_plot) ',Clusters ' num2str(cluster_highexp)])                    
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
        title(['Clusters highly expressed more times: ' num2str(clust_max_exp')]);
        colorbar
        colormap('jet')
        legend
        
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
        title(['For the cumulative expression, cluster more highly expressed' num2str(cluster_highexp)])
        colorbar
        colormap('jet')
        
        
    end
    
end