

function p_vals_perc_cells=func_perc_cells_cluster(all_genes_to_explore,list_of_genes,all_expr_log,...
    all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes,flag_comp_log,...
    flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
    flag_plot,number_plot,flag_min_pval,flag_save,labels,th_prop_cells,num_dim_tsne,obj,...
    num_perm,name_group_genes,az,el,folder_res,day_comp,all_lims,flag_nothing)

flag_comp_indCluster=1; % To 1 if we only want to compare one cluster vs the rest

if flag_comp_indCluster
    table_pvalues=ones(each_num_comp(2,1),size(all_genes_to_explore,1));
    table_pvalues_cum_exp=ones(each_num_comp(2,1),3);
    table_comparison_and_ind=cell(each_num_comp(2,1),4);
else
    table_pvalues=ones(number_comb,size(all_genes_to_explore,1));
    table_pvalues_cum_exp=ones(number_comb,3);
    table_comparison_and_ind=cell(number_comb,4);
end
i_g=3;
colors_label=[1,0,0;0,0,1]; % Red mutant, blue control

if flag_comp_log
    data_group=all_expr_log{1,i_g};
else
    data_group=all_expr_sca_cluster{1,i_g};
end
indexes_assign=all_assignments{1,1};

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
if flag_log_cum
    Both_cum_effect_all_genes_st=Both_cum_effect_all_genes_log;
else
    Both_cum_effect_all_genes_st=Both_cum_effect_all_genes;
end

% Calculate p values for all genes, and for the cumulative expression
if flag_comp_indCluster
    flag_cum_run=0;
    for i_gene=1:size(all_genes_to_explore,1)
        
        all_list_gene=all_genes_to_explore{i_gene};
        ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            indexes_all_genes(i_gene)=ind_gene;
            
            count_save=1;                
            comp_to_do=combntns(1:each_num_comp(2,1),1);
            comp_to_do=cat(2,comp_to_do,zeros(size(comp_to_do,1),each_num_comp(2,1)-1));
            
            for i_aux=1:size(comp_to_do,1)
                comp_to_do(i_aux,2:each_num_comp(2,1))=setdiff(1:each_num_comp(2,1),...
                    comp_to_do(i_aux,1));
                
                table_comparison_and_ind{count_save,1}=comp_to_do(i_aux,1);
                table_comparison_and_ind{count_save,2}=comp_to_do(i_aux,2:end);
                
                max_ind=max(indexes_assign);
                aux_indexes_assign=changem(indexes_assign,kron(max_ind+1,ones(1,1)),comp_to_do(i_aux,1));
                aux_indexes_assign=changem(aux_indexes_assign,kron(max_ind+2,ones(1,each_num_comp(2,1)-1)),...
                    comp_to_do(i_aux,2:end));
                aux_indexes_assign(aux_indexes_assign==(max_ind+1))=1;
                aux_indexes_assign(aux_indexes_assign==(max_ind+2))=2;
                
                table_comparison_and_ind{count_save,3}=aux_indexes_assign;
                
                vec_comp1=data_group(ind_gene,(aux_indexes_assign==1));
                vec_comp2=data_group(ind_gene,(aux_indexes_assign==2));
                
                %                 Both_cum_effect_all_genes=Both_cum_effect_all_genes+...
                %                     data_group(ind_gene,:);
                %
                table_pvalues(count_save,i_gene)=ranksum(vec_comp1,vec_comp2);
                if flag_cum_run==0                    
                    if mean(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)))>...
                            mean(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)))
                         table_pvalues_cum_exp(count_save,1)=ranksum(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)),...
                             Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)));
%                        table_pvalues_cum_exp(count_save,1)=permtest(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)),...
%                            Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)),20000,'conservative');                        
                    else
                        table_pvalues_cum_exp(count_save,1)=1;
                    end
                    %aux_p=permtest(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)),...
                    %    Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)),20000,'conservative');
                    aux_c1=table_comparison_and_ind{count_save,1};
                    aux_c2=table_comparison_and_ind{count_save,2};
                    table_pvalues_cum_exp(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2))';
                    %table_pvalues_cum_exp(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2)./sum(number_group_per_cluster,2))';
                    %table_pvalues_cum_exp(count_save,3)=sum(number_group_per_cluster(:,aux_c2),2)./sum(number_group_per_cluster,2);
                end
                
                count_save=count_save+1;
            end
            flag_cum_run=1;
            
        end
        
    end
    
else
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
                    
                    %                 Both_cum_effect_all_genes=Both_cum_effect_all_genes+...
                    %                     data_group(ind_gene,:);
                    %
                    table_pvalues(count_save,i_gene)=ranksum(vec_comp1,vec_comp2);
                    if i_gene==1
                        table_pvalues_cum_exp(count_save,1)=ranksum(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)),...
                            Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)));

                        %aux_p=permtest(Both_cum_effect_all_genes_st(1,(aux_indexes_assign==1)),...
                        %    Both_cum_effect_all_genes_st(1,(aux_indexes_assign==2)),20000,'conservative');
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
    
    table_pvalues_cum_exp_aux=table_pvalues_cum_exp;
    if flag_min_pval
        [p_val_asc,ind_asc]=sort(table_pvalues_cum_exp_aux(:,1),'ascend');
        row_to_plot=ind_asc(1:min(size(table_pvalues_cum_exp_aux,1),number_plot));
    else
        ind_val=find(table_pvalues_cum_exp_aux(:,6));
        row_to_plot=ind_val(table_pvalues_cum_exp_aux(ind_val,1)<thres_pval);
    end
    table_comparison_and_ind=all_pvalues_and_comp{2,1};
    
    %     [p_val_asc,ind_asc]=sort(table_pvalues(:),'ascend');
    %     ind_thres=max(find((p_val_asc<thres_pval)==1));
    %     number_plot=min(ind_thres,max_num_genes);
    %     ind_to_plot=ind_asc(1:number_plot);
    %     [row_to_plot,col_to_plot]=ind2sub(size(table_pvalues),ind_to_plot);
    
    for i_p=1:length(row_to_plot)
        
        p_val_plot=table_pvalues_cum_exp_aux(row_to_plot(i_p),1);
        indexes_comp=table_comparison_and_ind{row_to_plot(i_p),3};
        data_group_expr=all_expr{1,i_g};
        data_group=all_expr_sca_cluster{1,i_g};
        data_group_expr_log=all_expr_log{1,i_g};
        
        if mean(Both_cum_effect_all_genes(indexes_comp==2))>mean(Both_cum_effect_all_genes(indexes_comp==1))
            aux_indexes_comp=indexes_comp;
            indexes_comp=zeros(1,length(aux_indexes_comp));
            indexes_comp(aux_indexes_comp==1)=2;
            indexes_comp(aux_indexes_comp==2)=1;
        end
        
        %figure;
        FigHandle = figure('Position', [100, 100, 1200, 900],'PaperOrientation','landscape');
        h1=subplot(2,2,1);
        if num_dim_tsne==2
            gscatter(data_group(:,1),data_group(:,2),indexes_comp);
        elseif num_dim_tsne==3
            %gscatter3_custom(data_group,indexes_comp,h1);
            colors_exp=[1 0 0;0.7 0.7 0.7];
            scatter_colorCoded(data_group,colors_exp,indexes_comp-1,num_dim_tsne);
            ind_high_exp=find(indexes_comp==1);
            number_mut_cells=sum(labels(ind_high_exp)==0);
            number_control_cells=sum(labels(ind_high_exp)==1);
            title([num2str(number_mut_cells) '/' num2str(number_mut_cells/sum(labels==0)) ' Mutants, and '...
                num2str(number_control_cells) '/' num2str(number_control_cells/sum(labels==1)) ' Controls. Mutants = '...
                num2str(sum(labels==0)) ', Controls = ' num2str(sum(labels==1))],'FontSize',10)
        end
        view(az,el)
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        xlim(all_lims(1,:))
        ylim(all_lims(2,:))
        zlim(all_lims(3,:))
        box on
        subplot(2,2,2)
        Both_cum_effect_all_genes_norm=round((Both_cum_effect_all_genes/max(Both_cum_effect_all_genes))*255);
        if ~(sum(isnan(Both_cum_effect_all_genes_norm))>0)
            
            colorsA=jet(max(Both_cum_effect_all_genes_norm)+1);
            scatter_colorCoded(data_group,colorsA,Both_cum_effect_all_genes_norm,num_dim_tsne);
            %plot_gaussians(data_group,obj,indexes_assign,h1,num_dim_tsne);
            title(['P value of ' num2str(p_val_plot)])
            %colorbar;
            %colormap('jet');
            view(az,el)
            xticklabels(''); yticklabels('');zticklabels('');
            xticks([]);yticks([]);zticks([]);
            axis square
            xlim(all_lims(1,:))
            ylim(all_lims(2,:))
            zlim(all_lims(3,:))
            box on
            subplot(2,2,3)
            scatter_colorCoded(data_group,colors_label,labels,num_dim_tsne);
            view(az,el)
            xticklabels(''); yticklabels('');zticklabels('');
            xticks([]);yticks([]);zticks([]);
            axis square
            xlim(all_lims(1,:))
            ylim(all_lims(2,:))
            zlim(all_lims(3,:))
            box on
            %legend('Mutant','Control')
            subplot(2,2,4)
            histogram(Both_cum_effect_all_genes(indexes_comp==1),20,'Normalization','probability','FaceColor','g');
            hold on;
            histogram(Both_cum_effect_all_genes(indexes_comp==2),20,'Normalization','probability','FaceColor','m');
            title(['Distributions for cumulative effect of all genes in ' name_group_genes])
            %legend('Highly expressed','Rest of clusters')
        else
            stop=1;
        end
        if flag_save
            print([folder_res 'Comp_for' name_group_genes '_Day' num2str(day_comp) '_' ...
                '_Comp' num2str(row_to_plot(i_p)) 'Pval' num2str(p_val_plot) '.pdf'], '-dpdf','-fillpage')
        end
        
        % A figure conveying all the information
        if ~(sum(isnan(Both_cum_effect_all_genes_norm))>0)
        FigHandle = figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
        Both_cum_round=round((Both_cum_effect_all_genes/max(Both_cum_effect_all_genes))*10);
        sizes_markers=(1:(max(Both_cum_round)+1))+1;
        %colors_exp=[0,0.8,0.2;0.8,0.2,0.2];
        color_hot=autumn(max(Both_cum_round)+1);
        color_hot=1-color_hot.^2;
        color_hot=color_hot(end:-1:1,3:-1:1);
        color_cool=gray(max(Both_cum_round)+4);
        color_cool=color_cool((end-3):-1:1,:);
        colors_exp=[color_hot;color_cool];
        
        if size(Both_cum_round,1)~=size(indexes_comp,1)
            indexes_comp=indexes_comp';
        end
        
        markers_g={'s','o'};
        ind_g0=find(labels==0);
        ind_g1=find(labels==1);
        %         scatter3(data_group(ind_g0,1),data_group(ind_g0,2),data_group(ind_g0,3),...
        %             sizes_markers(Both_cum_round(ind_g0)+1)*20,colors_exp(indexes_comp(ind_g0),:),markers_g{1},'filled');
        %         hold on
        %         scatter3(data_group(ind_g1,1),data_group(ind_g1,2),data_group(ind_g1,3),...
        %             sizes_markers(Both_cum_round(ind_g1)+1)*20,colors_exp(indexes_comp(ind_g1),:)/2,markers_g{2});
        scatter3(data_group(ind_g0,1),data_group(ind_g0,2),data_group(ind_g0,3),...
            sizes_markers(Both_cum_round(ind_g0)+1)*20,colors_exp(size(color_hot,1)*(indexes_comp(ind_g0)-1)+(Both_cum_round(ind_g0)+1),:),...
            markers_g{1},'filled');
        hold on
        scatter3(data_group(ind_g1,1),data_group(ind_g1,2),data_group(ind_g1,3),...
            sizes_markers(Both_cum_round(ind_g1)+1)*20,colors_exp(size(color_hot,1)*(indexes_comp(ind_g1)-1)+(Both_cum_round(ind_g1)+1),:),...
            markers_g{2});
        view(az,el)
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        xlim(all_lims(1,:))
        ylim(all_lims(2,:))
        zlim(all_lims(3,:))
        box on
        title([num2str(number_mut_cells) '/' num2str(number_mut_cells/sum(labels==0)) ' Mutants, and '...
            num2str(number_control_cells) '/' num2str(number_control_cells/sum(labels==1)) ' Controls. Mutants = '...
            num2str(sum(labels==0)) ', Controls = ' num2str(sum(labels==1))],'FontSize',10)
        if flag_nothing
            grid off
            axis off
            box off
            title('')
            legend off
        end
        if flag_save
            print([folder_res 'CompAll_for' name_group_genes '_Day' num2str(day_comp) '_' ...
                'AllGroups_Smutant_Ocontrol.pdf'], '-dpdf','-fillpage')
        end
        end 
        %
        if i_p==1
            indexes_comp_end=indexes_comp;
            number_mut_cells=table_pvalues_cum_exp(row_to_plot(i_p),2);
            number_control_cells=table_pvalues_cum_exp(row_to_plot(i_p),3);
            ind_comp_end=row_to_plot(i_p);
        end
        
        
    end
    %close all
    
    % Now, we randomly permutate the labels for obtaining a p value for the
    % concentration of specific cells in that clusters. It works as
    % follows:
    % 1. From the previous experiment, I have a value for the number of
    % mutant and controls cells in the clusters
    % 2. Then, I run M permutations of the labels, checking now how many
    % cells lie on that clusters, and saving those
    % 3. Then, I have a distribution, and I cannot know how many times I
    % get more cells on that cluster that for the case computed in 1. This
    % cases that lie over 1. are a total number of N
    % 4. Then the p value is N/M
    table_num_cells_perm=zeros(num_perm,2); % Mutant label 0 and control label 1
    aux_indexes_assign=table_comparison_and_ind{ind_comp_end,3};
    if mean(Both_cum_effect_all_genes(aux_indexes_assign==1))>mean(Both_cum_effect_all_genes(aux_indexes_assign==2))
        cluster_highexp=table_comparison_and_ind{ind_comp_end,1};
    else
        cluster_highexp=table_comparison_and_ind{ind_comp_end,2};
    end
    ind_cells_clust=[];
    for i_clust=1:length(cluster_highexp)
        ind_cells_clust=[ind_cells_clust,find(indexes_assign==cluster_highexp(i_clust))'];
    end
    ind_cells_clust=sort(ind_cells_clust);
    
    for i_perm=1:num_perm
        labels_perm=labels(randperm(length(labels)));
        aux_labels_chosen=labels_perm(ind_cells_clust);
        table_num_cells_perm(i_perm,:)=[sum(aux_labels_chosen==0),sum(aux_labels_chosen==1)];
    end
    number_mut_cells=sum(labels(ind_cells_clust)==0);
    number_control_cells=sum(labels(ind_cells_clust)==1);
    p_val_mut=sum(table_num_cells_perm(:,1)>number_mut_cells)/num_perm;
    p_val_control=sum(table_num_cells_perm(:,2)>number_control_cells)/num_perm;
    
    p_vals_perc_cells{1,1}=table_num_cells_perm;
    p_vals_perc_cells{1,2}=[number_mut_cells,number_control_cells];
    p_vals_perc_cells{1,3}=[p_val_mut,p_val_control];
    p_vals_perc_cells{1,4}=cluster_highexp;
    
    if flag_plot
        h1 = figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
        histogram(table_num_cells_perm(:,1),'Normalization','probability')
        hold on
        histogram(table_num_cells_perm(:,2),'Normalization','probability')
        legend('Distribution for mutant','Distribution for control')
        plot([number_mut_cells,number_mut_cells],[0,0.6],'b','LineWidth',6);
        plot([number_control_cells,number_control_cells],[0,0.6],'r','LineWidth',6);
        title(['Pval for Mutant: ' num2str(p_val_mut) '/ Pval for Control: ' num2str(p_val_control)])
      
        if flag_save
            print([folder_res 'Distrib_and_Pvalues_for' name_group_genes '_Day' num2str(day_comp)...
                '.pdf'], '-dpdf','-fillpage')
        end
    end
    
    %     indexes_comp_end
    %     number_mut_cells
    %     number_control_cells
    
    
end