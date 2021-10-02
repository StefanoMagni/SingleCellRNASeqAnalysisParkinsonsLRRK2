


function [both_cum_expression_all_corr,data_group_log_red,all_genes_to_check,index_group_genes]=...
    func_correlation_lists(list_correlation,folder_lists,list_of_genes,all_expr_log,...
    all_assignments,flag_comb_cumgenes,flag_plot,folder_save,day_comp,...
    flag_save)

indexes_clusters=all_assignments{1,1};

all_genes_to_check={};
all_indexes_gene=[];
index_group_genes=[];
genes_per_group=[];
for i_list=1:length(list_correlation)
    
    all_genes_to_explore=write_struct_fromtxt([folder_lists list_correlation{i_list}]);
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

data_group_log=all_expr_log{1,3};
data_group_log_red=data_group_log(all_indexes_gene,:);

both_cum_expression_all_corr=zeros(length(list_correlation),size(data_group_log_red,2));

for i_list=1:length(list_correlation)
    
    ind_list_ph=find(index_group_genes==i_list);
    
    if flag_comb_cumgenes==1 % 1 for sum, 2 for median, 3 for mean
        both_cum_expression=sum(data_group_log_red(ind_list_ph,:));
    elseif flag_comb_cumgenes==2
        both_cum_expression=median(data_group_log_red(ind_list_ph,:));
    elseif flag_comb_cumgenes==3
        both_cum_expression=mean(data_group_log_red(ind_list_ph,:));
    end
    both_cum_expression_all_corr(i_list,:)=both_cum_expression;
    
end

aux_list={};
for i_l=1:length(list_correlation)
   
    aux_name=list_correlation{1,i_l};
    aux_name(strfind(aux_name,'_'))=' ';
    aux_name(strfind(aux_name,'.txt')+(0:3))=' ';
    aux_name(strfind(aux_name,'list ')+(0:4))='';    
    aux_list{1,i_l}=aux_name;
    
end

if flag_plot
% Now, we obtain two matrices of correlation, for the cumulative and for
% the individual genes
cov_genes = corr(data_group_log_red');
MYredgreencmap = redgreencmap.^10;

% clustergram(data_group_log_red, 'RowLabels', all_genes_to_check, 'DisplayRange', 1,...
%     'ClusterValue', 1, 'ColorMap', MYredgreencmap)
% 
% clustergram(cov_genes, 'RowLabels', all_genes_to_check,...
%     'ClusterValue', 1, 'ColorMap', MYredgreencmap)

cov_genes_gr = corr(both_cum_expression_all_corr');

% clustergram(cov_genes_gr, 'RowLabels', aux_list,...
%     'DisplayRange', 1, 'ColorMap', MYredgreencmap)
% colors_aux=hsv(length(list_correlation));
% objCol.Labels=aux_list;
% objCol.Colors=colors_aux;
% clustergram(cov_genes_gr, 'RowLabels', aux_list, 'ColumnLabels', aux_list,...
%     'DisplayRange', 1, 'RowLabelsColor', objCol, 'ColorMap', MYredgreencmap)

%HeatMap(data_group_log_red, 'RowLabels', all_genes_to_check, 'ColorMap', MYredgreencmap)
%HeatMap(both_cum_expression_all_corr, 'RowLabels', aux_list, 'ColorMap', MYredgreencmap)
h1=HeatMap(cov_genes, 'RowLabels', all_genes_to_check, 'ColumnLabels', all_genes_to_check, 'ColorMap', MYredgreencmap)
plot(h1)
%title('Correlation coefficient for all the genes')
if flag_save
    print([folder_save 'Corr_allGenes_' num2str(day_comp) '.pdf'],'-dpdf')
end

h1=HeatMap(cov_genes_gr, 'RowLabels', aux_list, 'ColumnLabels', aux_list, 'ColorMap', MYredgreencmap)
plot(h1)
%title('Correlation coefficient for the cumulative expressions')
if flag_save
    print([folder_save 'Corr_CumExpress_' num2str(day_comp) '.pdf'],'-dpdf')
end

end
