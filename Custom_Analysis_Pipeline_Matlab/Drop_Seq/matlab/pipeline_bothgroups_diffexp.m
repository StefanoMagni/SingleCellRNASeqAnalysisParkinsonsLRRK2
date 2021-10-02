
function pipeline_bothgroups_diffexp(folder_files,day_comp,day_comp_index,...
    name_groups,numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,...
    flag_save,flag_useTSNE,flag_comp_log,folder_res,max_clusters,all_genes_to_explore,...
    thres_pval,th_prop_cells,flag_min_pval)

number_plot=10;

if 1
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
    [~,score,eigval,~,~,~]=pca(Both_expr_sca');
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
    
    for i_g=1:length(list_genes_plot)
        
        gene_colorcode=list_genes_plot{i_g};
        
        ind_gene=cellfun(@(s) (strcmp(s,gene_colorcode)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        colorsA=jet(max(A_expr_max(ind_gene,:))+1);
        colorsB=jet(max(B_expr_max(ind_gene,:))+1);
        colorsBoth=jet(max([A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)])+1);
        
        if flag_plot
            figure(1)
            scatter_colorCoded(Both_expr_sca_tsne,colorsBoth,...
                cat(2,A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)),num_dim_tsne);
            colorbar
            colormap('jet')
            figure(2)
            scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
        end
        % Color code by label
        %scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
        if flag_save
            figure(1)
            print([folder_res 'Both_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            figure(2)
            print([folder_res 'Both_Day' num2str(day_comp) '_Coloured_by_Group_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            
        end
        
        if num_dim_tsne==2
            close all
        end
        
    end
    
    if flag_useTSNE
        all_expr_sca_cluster{1,1}=A_expr_sca_tsne;
        all_expr_sca_cluster{1,2}=B_expr_sca_tsne;
        all_expr_sca_cluster{1,3}=Both_expr_sca_tsne;
    else
        all_expr_sca_cluster{1,1}=Both_expr_sca_trans(1:number_pcs,:)';
        all_expr_sca_cluster{1,2}=Both_expr_sca_trans(1:number_pcs,:)';
        all_expr_sca_cluster{1,3}=Both_expr_sca_trans(1:number_pcs,:)';
    end
    
    all_expr_log{1,1}=log10(A_expr_max+1);
    all_expr_log{1,2}=log10(B_expr_max+1);
    all_expr_log{1,3}=log10(Both_expr_max+1);
    all_expr{1,1}=A_expr_max;
    all_expr{1,2}=B_expr_max;
    all_expr{1,3}=Both_expr_max;
    
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
        indexes=all_assignments{1,1};
        obj=all_structures_gm{each_num_comp(2,1),1};
        h1=figure;
        if num_dim_tsne==2
            h1 = gscatter(data_group(:,1),data_group(:,2),indexes);
        elseif num_dim_tsne==3
            gscatter3_custom(data_group,indexes,h1);
        end
        
        plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne);
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
    for i=1:max(indexes)
        aux_ind=find(indexes==i);
        percent=sum(labels(aux_ind)==0)/length(aux_ind);
        number_group_per_cluster(:,i)=[sum(labels(aux_ind)==0);sum(labels(aux_ind)==1)];
        percent_group_per_cluster(:,i)=[percent;1-percent]*100;
    end
    
else
    load('all_data_test_both.mat')
end

%% Now, using these cluster, explore all possibilities of comparison between clusters
% For this, we will use non-parametric testing, and we will proceed gene by
% gene. For those most significative genes below a threshold, we will then 

all_pvalues_and_comp=cell(3,1);
indexes_all_genes=zeros(1,size(all_genes_to_explore,1));
    
number_comb=0;
for i=1:floor(each_num_comp(2,1)/2)
    if i==each_num_comp(2,1)/2
        number_comb=number_comb+(nchoosek(each_num_comp(2,1),i)/2);
    else
        number_comb=number_comb+nchoosek(each_num_comp(2,1),i);
    end
end
%     all_perms=perms(1:each_num_comp(2,i_g));
%     all_perms=all_perms(end:-1:1,:);

table_pvalues=zeros(number_comb,size(all_genes_to_explore,1));
table_pvalues_cum_exp=zeros(number_comb,3);
table_comparison_and_ind=cell(number_comb,4);

if flag_comp_log
    data_group=all_expr_log{1,i_g};
else
    data_group=all_expr_sca_cluster{1,i_g};
end
indexes_assign=all_assignments{1,1};

% Get the cumulative effect for the list of genes
Both_cum_effect_all_genes=zeros(1,length(labels));
for i_gene=1:size(all_genes_to_explore,1)
    all_list_gene=all_genes_to_explore{i_gene};
    ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
    ind_gene=find(ind_gene==1);
    
    if ~isempty(ind_gene)
        Both_cum_effect_all_genes=Both_cum_effect_all_genes+...
            10.^data_group(ind_gene,:);
    end
end
Both_cum_effect_all_genes=log10(Both_cum_effect_all_genes);

for i_gene=1:size(all_genes_to_explore,1)
    
    all_list_gene=all_genes_to_explore{i_gene};
    ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
    ind_gene=find(ind_gene==1);
    
    if ~isempty(ind_gene)
        indexes_all_genes(i_gene)=ind_gene;
        
        count_save=1;
        for i_comb1=1:floor(each_num_comp(2,1)/2)
            
            comp_to_do=combntns(1:each_num_comp(2,1),i_comb1);
            if i_comb1==floor(each_num_comp(2,1)/2)
                comp_to_do=comp_to_do(1:(size(comp_to_do,1)/2),:);
            end
            comp_to_do=cat(2,comp_to_do,zeros(size(comp_to_do,1),each_num_comp(2,1)-size(comp_to_do,2)));
            for i_aux=1:size(comp_to_do,1)
                comp_to_do(i_aux,(i_comb1+1):each_num_comp(2,1))=setdiff(1:each_num_comp(2,1),...
                    comp_to_do(i_aux,1:i_comb1));
                
                table_comparison_and_ind{count_save,1}=comp_to_do(i_aux,1:i_comb1);
                table_comparison_and_ind{count_save,2}=comp_to_do(i_aux,(i_comb1+1):end);
                
                max_ind=max(indexes_assign);
                aux_indexes_assign=changem(indexes_assign,kron(max_ind+1,ones(1,i_comb1)),comp_to_do(i_aux,1:i_comb1));
                aux_indexes_assign=changem(aux_indexes_assign,kron(max_ind+2,ones(1,each_num_comp(2,1)-i_comb1)),...
                    comp_to_do(i_aux,(i_comb1+1):end));
                aux_indexes_assign(aux_indexes_assign==(max_ind+1))=1;
                aux_indexes_assign(aux_indexes_assign==(max_ind+2))=2;
                
                table_comparison_and_ind{count_save,3}=aux_indexes_assign;
                
                vec_comp1=data_group(ind_gene,(aux_indexes_assign==1));
                vec_comp2=data_group(ind_gene,(aux_indexes_assign==2));
                
                Both_cum_effect_all_genes=Both_cum_effect_all_genes+...
                    data_group(ind_gene,:);
                
                table_pvalues(count_save,i_gene)=ranksum(vec_comp1,vec_comp2);
                if i_gene==1
                    table_pvalues_cum_exp(count_save,1)=ranksum(Both_cum_effect_all_genes(1,(aux_indexes_assign==1)),...
                        Both_cum_effect_all_genes(1,(aux_indexes_assign==2)));
                    aux_c1=table_comparison_and_ind{count_save,1};
                    aux_c2=table_comparison_and_ind{count_save,2};
                    table_pvalues_cum_exp(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2))';
                    %table_pvalues_cum_exp(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2)./sum(number_group_per_cluster,2))';
                    %table_pvalues_cum_exp(count_save,3)=sum(number_group_per_cluster(:,aux_c2),2)./sum(number_group_per_cluster,2);                    
                end
                
                count_save=count_save+1;
            end
        end
        
    end
    
end

table_pvalues_cum_exp(:,4)=table_pvalues_cum_exp(:,2)./table_pvalues_cum_exp(:,3);
table_pvalues_cum_exp(:,5)=table_pvalues_cum_exp(:,3)./table_pvalues_cum_exp(:,2);
table_pvalues_cum_exp(:,6)=(table_pvalues_cum_exp(:,4)>th_prop_cells)...
    +(table_pvalues_cum_exp(:,5)>th_prop_cells);

all_pvalues_and_comp{1,1}=table_pvalues;
all_pvalues_and_comp{2,1}=table_comparison_and_ind;
all_pvalues_and_comp{3,1}=table_pvalues_cum_exp;

%% Now, we just choose those cases with p_values below some threshold and plot them

if flag_plot
    
    if flag_min_pval
        [p_val_asc,ind_asc]=sort(table_pvalues_cum_exp(:,1),'ascend');
        row_to_plot=ind_asc(1:min(size(table_pvalues_cum_exp,1),number_plot));
    else
        ind_val=find(table_pvalues_cum_exp(:,6));
        row_to_plot=ind_val(table_pvalues_cum_exp(ind_val,1)<thres_pval);        
    end
    table_comparison_and_ind=all_pvalues_and_comp{2,1};
        
%     [p_val_asc,ind_asc]=sort(table_pvalues(:),'ascend');
%     ind_thres=max(find((p_val_asc<thres_pval)==1));
%     number_plot=min(ind_thres,max_num_genes);
%     ind_to_plot=ind_asc(1:number_plot);
%     [row_to_plot,col_to_plot]=ind2sub(size(table_pvalues),ind_to_plot);        
    
    for i_p=1:length(row_to_plot)
       
        p_val_plot=table_pvalues_cum_exp(row_to_plot(i_p),1);        
        indexes_comp=table_comparison_and_ind{row_to_plot(i_p),3};
        data_group_expr=all_expr{1,i_g};
        data_group=all_expr_sca_cluster{1,i_g};
        data_group_expr_log=all_expr_log{1,i_g};
        
        %figure;
        FigHandle = figure('Position', [100, 100, 1200, 900]);
        h1=subplot(2,2,1);
        if num_dim_tsne==2
            gscatter(data_group(:,1),data_group(:,2),indexes_comp);
        elseif num_dim_tsne==3
            gscatter3_custom(data_group,indexes_comp,h1); 
            title(['Cluster with ' num2str(table_pvalues_cum_exp(row_to_plot(i_p),2)) ' Mutant cells and '...
                num2str(table_pvalues_cum_exp(row_to_plot(i_p),3)) ' Control cells. Percentages of '...
                num2str(table_pvalues_cum_exp(row_to_plot(i_p),2)/sum(labels==0)) ' Mutant and '...
                num2str(table_pvalues_cum_exp(row_to_plot(i_p),3)/sum(labels==1)) ' Control'])
        end
        subplot(2,2,2)
        Both_cum_effect_all_genes_norm=round((Both_cum_effect_all_genes/max(Both_cum_effect_all_genes))*255);
        colorsA=jet(max(Both_cum_effect_all_genes_norm)+1);
        scatter_colorCoded(data_group,colorsA,Both_cum_effect_all_genes_norm,num_dim_tsne);
        plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne);
        title(['P value of ' num2str(p_val_plot)])
        colorbar;
        colormap('jet');
        subplot(2,2,3)
        scatter_colorCoded(data_group,colors_label,labels,num_dim_tsne);
        subplot(2,2,4)
        histogram(Both_cum_effect_all_genes(indexes_comp==1),20,'Normalization','probability');
        hold on;
        histogram(Both_cum_effect_all_genes(indexes_comp==2),20,'Normalization','probability');
        if flag_save
            print([folder_res 'Comp_Day' num2str(day_comp) '_' ...
                '_Comp' row_to_plot(i_p) 'Pval' num2str(p_val_plot) '.png'], '-dpng', '-r300')
        end

        
    end
    close all

end

%% Now, using these cluster, explore all possibilities of comparison between clusters
% For this, we will use non-parametric testing, and we will proceed gene by
% gene. For those most significative genes below a threshold, we will then


    function plot_gaussians(X0,gmfit,index,h1,num_dim_tsne)
        
        %fcontour(@(x1,x2)pdf(gmfit,[x1 x2]),[xlim(1) ylim(1) xlim(2) ylim(2)])
        hold on;
        if num_dim_tsne==2
            maxminXY=get(gca,{'XLim','YLim'});
            xlim=maxminXY{1};
            ylim=maxminXY{2};
            ezcontour(@(x1,x2)pdf(gmfit,[x1 x2]),xlim,ylim)
            plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        elseif num_dim_tsne==3
            for i_k=1:size(gmfit.mu,1)
                hold on;
                alpha_val=0.2;
                plot_ellipsoid(gmfit.Sigma(:,:,i_k),gmfit.mu(i_k,:),alpha_val);
            end
            plot3(gmfit.mu(:,1),gmfit.mu(:,2),gmfit.mu(:,3),'kx','LineWidth',2,'MarkerSize',10)
        end
        
    end

    function h1=gscatter3_custom(data,indexes,h1)
        
        [~,~,id] = unique(indexes);
        colors = 'rgbmckyrgbmckyrgbmcky';
        markers = 'osdosdosdosdosdosd';        
        for idx = 1:max(id)
            data_aux=data(indexes == idx,:);
            plot3(data_aux(:,1), data_aux(:,2), data_aux(:,3), [colors(idx) markers(idx)]);
            hold on;
        end
        grid
        
    end


end
