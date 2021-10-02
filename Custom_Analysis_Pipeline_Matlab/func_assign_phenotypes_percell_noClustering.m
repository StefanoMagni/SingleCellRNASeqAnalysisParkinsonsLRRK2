


function [both_cum_expression_all_phenos,all_genes_to_check,index_group_genes,ind_max_pheno]=...
    func_assign_phenotypes_percell_noClustering(list_pheno,folder_pheno,list_of_genes,all_expr,all_expr_log,...
    flag_comb_cumgenes,flag_plot,num_dim_tsne,data_group,...
    colors,title_eachfig,az,el,flag_save,name_save,labels,alpha_val,colormap_val,...
    flag_log_cum,day_comp,flag_nothing)

flag_legend=0;

all_genes_to_check={};
all_indexes_gene=[];
index_group_genes=[];
genes_per_group=zeros(1,length(list_pheno));
genes_checked={};
for i_list=1:length(list_pheno)
    
    all_genes_to_explore=write_struct_fromtxt([folder_pheno list_pheno{i_list}]);
    % Only choose the genes present in the list_if_genes
    ind_val=[];
    count_g=1;
    for i_gene=1:size(all_genes_to_explore,1)
        all_list_gene=all_genes_to_explore{i_gene};
        
        if 1
            ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
            ind_gene=find(ind_gene==1);
        else
            ind_gene=cellfun(@(s) (strfind(s,all_list_gene)),list_of_genes,'UniformOutput',0);
            ind_gene=find(cellfun(@(s) ~(isempty(s)),ind_gene));
        end;
        
        if ~isempty(ind_gene)
            genes_checked{i_list,count_g}=list_of_genes{ind_gene};
            count_g=count_g+1;
            ind_val=cat(2,ind_val,i_gene);
            all_indexes_gene=cat(2,all_indexes_gene,ind_gene);
            index_group_genes=cat(2,index_group_genes,i_list);
            genes_per_group(i_list)=genes_per_group(i_list)+1;
        end
    end
    all_genes_to_explore=all_genes_to_explore(ind_val,1);
    all_genes_to_check=cat(1,all_genes_to_check,all_genes_to_explore);
end

data_group_log=all_expr_log{1,1};
data_group_log_red=data_group_log(all_indexes_gene,:);

data_group_expr=all_expr{1,1};
data_group_expr=data_group_expr(all_indexes_gene,:);

both_cum_expression_all_phenos=zeros(length(list_pheno),size(data_group_log_red,2));

for i_list=1:length(list_pheno)
    
    ind_list_ph=find(index_group_genes==i_list);
    
    if flag_comb_cumgenes==1 % 1 for sum, 2 for median, 3 for mean
        both_cum_expression=sum(data_group_expr(ind_list_ph,:));
    elseif flag_comb_cumgenes==2
        both_cum_expression=median(data_group_expr(ind_list_ph,:));
    elseif flag_comb_cumgenes==3
        both_cum_expression=mean(data_group_expr(ind_list_ph,:));
    end
    if flag_log_cum
        both_cum_expression_all_phenos(i_list,:)=log10(both_cum_expression+1);
    else
        both_cum_expression_all_phenos(i_list,:)=both_cum_expression;
    end
    
    
end

[max_val,ind_max_pheno]=max(both_cum_expression_all_phenos);

perc_pheno_total_group=zeros(3,length(list_pheno));
ind_g0=find(labels==0);
ind_g1=find(labels==1);
for i_ph=1:length(list_pheno)
    perc_pheno_total_group(1,i_ph)=sum(ind_max_pheno==i_ph)/length(ind_max_pheno);
    perc_pheno_total_group(2,i_ph)=sum(ind_max_pheno(ind_g0)==i_ph)/sum(ind_max_pheno==i_ph);
    perc_pheno_total_group(3,i_ph)=sum(ind_max_pheno(ind_g1)==i_ph)/sum(ind_max_pheno==i_ph);
end

fprintf(['Proportions of cells for pheno, total and per group, Day ' num2str(day_comp) '\n'])
fprintf('%f ',perc_pheno_total_group(1,:))
fprintf('\n');
fprintf('%f ',perc_pheno_total_group(2,:))
fprintf('\n');
fprintf('%f ',perc_pheno_total_group(3,:))
fprintf('\n');

%sum(ind_max_pheno==3);
%N = histc(ind_max_pheno,0.5:8.5)
%list_pheno{3}

if 0
   
    both_exp_all_phenos_A=both_cum_expression_all_phenos(:,labels==0);
    both_exp_all_phenos_A=array2table(both_exp_all_phenos_A,'RowNames',list_pheno);
    writetable(both_exp_all_phenos_A,['A_' num2str(day_comp) '_CumExp_Pheno.txt'],'WriteRowNames',1)
    writetable(both_exp_all_phenos_A,['A_' num2str(day_comp) '_CumExp_Pheno.csv'],'WriteRowNames',1)
    
    both_exp_all_phenos_B=both_cum_expression_all_phenos(:,labels==1);
    both_exp_all_phenos_B=array2table(both_exp_all_phenos_B,'RowNames',list_pheno);
    writetable(both_exp_all_phenos_B,['B_' num2str(day_comp) '_CumExp_Pheno.txt'],'WriteRowNames',1)
    writetable(both_exp_all_phenos_B,['B_' num2str(day_comp) '_CumExp_Pheno.csv'],'WriteRowNames',1)

    ind_max_pheno_A=ind_max_pheno(:,labels==0);
    ind_max_pheno_A=array2table(ind_max_pheno_A);
    writetable(ind_max_pheno_A,['A_' num2str(day_comp) '_Highest_Pheno.txt'],'WriteRowNames',1)
    writetable(ind_max_pheno_A,['A_' num2str(day_comp) '_Highest_Pheno.csv'],'WriteRowNames',1)
    
    ind_max_pheno_B=ind_max_pheno(:,labels==1);
    ind_max_pheno_B=array2table(ind_max_pheno_B);
    writetable(ind_max_pheno_B,['B_' num2str(day_comp) '_Highest_Pheno.txt'],'WriteRowNames',1)
    writetable(ind_max_pheno_B,['B_' num2str(day_comp) '_Highest_Pheno.csv'],'WriteRowNames',1)
    
    data_group_exprs_A=data_group_expr(:,labels==0);
    data_group_exprs_A=array2table(data_group_exprs_A,'RowNames',list_of_genes);
    writetable(data_group_exprs_A,['A_' num2str(day_comp) '_Red_GenExprMat.txt'],'WriteRowNames',1)
    %writetable(data_group_exprs_A,['A_' num2str(day_comp) '_CumExp_Pheno.csv'],'WriteRowNames',1)
    
    data_group_expr_B=data_group_expr(:,labels==1);
    data_group_expr_B=array2table(data_group_expr_B,'RowNames',list_of_genes);
    writetable(data_group_expr_B,['B_' num2str(day_comp) '_Red_GenExprMat.txt'],'WriteRowNames',1)
    %writetable(data_group_expr_B,['B_' num2str(day_comp) '_CumExp_Pheno.csv'],'WriteRowNames',1)

end

if flag_plot
    
    % Change the color for a better visualization
    colors([5,8],:)=colors([8,5],:);
    
    % Now, we also directly plot the results and save them
    % Generate the figure
    h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
    % To obtain a first plot of single points to create the legend
    for i_leg=1:size(colors,1)
        aux_point=data_group(1,:);
        scatter3(aux_point(1,1),aux_point(1,2),aux_point(1,3),1,colors(i_leg,:),'filled')
        hold on
        if i_leg==1
            leg_cluster={'No particular phenotype'};
        else
            aux_list=list_pheno{i_leg-1};
            aux_list(strfind(aux_list,'_'))=' ';
            aux_list(strfind(aux_list,'.txt')+(0:3))=' ';
            leg_cluster{1,i_leg}=aux_list;
        end
    end
    view(az,el)
    % Plot the figure for the cumulative effect
    % HERE PERHAPS WE CAN CHOOSE THE ONE WITH P VALUE THAT SATISFIES,
    % AND WITH THE HIGHEST CUMULATIVE EXPRESSION
    if flag_legend
        leg_ob=legend(leg_cluster);
    end
    for i_ph=1:length(list_pheno)
        ind_pheno=find(ind_max_pheno==i_ph);

        group_belong=i_ph+1;

        aux_point=data_group(ind_pheno,:);
        if isempty(labels)
            scatter3(aux_point(:,1),aux_point(:,2),aux_point(:,3),30,colors(group_belong,:),'filled')
        else
            ind_cl0=ind_pheno((labels(ind_pheno)==0));
            ind_cl1=ind_pheno((labels(ind_pheno)==1));
            plot3(data_group(ind_cl0,1),data_group(ind_cl0,2),...
                data_group(ind_cl0,3),'s','MarkerEdgeColor',colors(group_belong,:),...
                'MarkerFaceColor',[1,1,1],'MarkerSize',6)
            hold on
            plot3(data_group(ind_cl1,1),data_group(ind_cl1,2),...
                data_group(ind_cl1,3),'o','MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',colors(group_belong,:),'MarkerSize',6)
        end
    end
    xticklabels(''); yticklabels('');zticklabels('');
    xticks([]);yticks([]);zticks([]);
    axis square
%     xlim(all_lims(1,:))
%     ylim(all_lims(2,:))
%     zlim(all_lims(3,:))
    box on
    title(title_eachfig)
    if flag_legend
        leg_ob.String=leg_ob.String(1:(length(list_pheno)+1));
    end
    if flag_nothing
        grid off
        axis off
        box off
        title('')
        legend off
    end
    if flag_save
        print(name_save, '-dpdf','-fillpage')
    end
end
