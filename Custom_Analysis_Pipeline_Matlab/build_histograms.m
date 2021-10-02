
folder_files='./jonas/';
day_comp=[0,10,14,42]; % [0,10,14,42]
day_comp_index=1:4; % 1:4
name_groups={'A','B'}; % A is mutated and B control
% The columns are A and B, and the rows 0, 10, 14 or 42 days
numb_cells=[500,500;250,250;505,350;400,300];

list_genes_pvals={'./Drop_Seq/list_CellCycle.txt','./Drop_Seq/list_stem.txt'};

folder_res='./Images_results/Exp1/';
max_clusters=6;

threshold_p=0.01;
min_tests_rej=2;

num_perm=10000;


%%%%%%
if 0
    
all_days_expr=cell(1,length(day_comp_index));
all_days_labels=cell(1,length(day_comp_index));
all_days_genes=cell(1,length(day_comp_index));
for i_f=1:length(day_comp_index)
    i_f
    [all_expr,labels,list_of_genes]=open_filter_genes(folder_files,...
        day_comp(i_f),i_f,name_groups,numb_cells);
    all_days_expr{1,i_f}=all_expr;
    all_days_labels{1,i_f}=labels;
    all_days_genes{1,i_f}=list_of_genes;
end

% Table of all p-values. The comparisons are:
% - For both groups, between days
% - For each day, between groups
% - For each group, against the previous days

struct_genes={};
num_genes_test=0;
for i_f=1:length(list_genes_pvals)
    aux_struct_genes={};
    file_genes=list_genes_pvals{i_f};
    all_genes_to_explore=write_struct_fromtxt(file_genes);
    empty_cols=sum(not(cellfun(@isempty,all_genes_to_explore)));
    all_genes_to_explore=all_genes_to_explore(:,find(empty_cols>0));
    [~,file_n,~]=fileparts(file_genes);    
    aux_struct_genes(1:length(all_genes_to_explore),2)={file_n};
    aux_struct_genes(1:length(all_genes_to_explore),1)=all_genes_to_explore;
    struct_genes=cat(1,struct_genes,aux_struct_genes);
end
num_genes_test=size(struct_genes,1);

%%%%
all_p_values_days=ones(num_genes_test,length(day_comp)-1);
for i_g=1:(length(day_comp)-1)
    i_g
    aux_all_expr1=all_days_expr{1,i_g}.log{1,3};
    aux_all_expr2=all_days_expr{1,i_g+1}.log{1,3};
    genes1=all_days_genes{1,i_g};
    genes2=all_days_genes{1,i_g+1};
    
    [list_of_genes,ind_inters1,ind_inters2]=intersect(genes1,genes2);

    aux_all_expr1=aux_all_expr1(ind_inters1,:);
    aux_all_expr2=aux_all_expr2(ind_inters2,:);
    tic
    for i_gene=1:num_genes_test
        
        gene_p=struct_genes{i_gene,1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            
            pval_gene=permtest(aux_all_expr1(ind_gene,:),aux_all_expr2(ind_gene,:),num_perm,'conservative');
            %pval_gene=ranksum(aux_all_expr1(ind_gene,:),aux_all_expr2(ind_gene,:));
            all_p_values_days(i_gene,i_g)=pval_gene;
            
        end
        
    end
    toc
end

%%%
all_p_values_groups=ones(num_genes_test,length(day_comp));
for i_g=1:length(day_comp)
    i_g
    aux_all_expr=all_days_expr{1,i_g}.log{1,3};
    list_of_genes=all_days_genes{1,i_g};    
    labels=all_days_labels{1,i_g};
    ind_g1=find(labels==0);
    ind_g2=find(labels==1);
    
    for i_gene=1:num_genes_test
        
        gene_p=struct_genes{i_gene,1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            
            pval_gene=permtest(aux_all_expr(ind_gene,ind_g1),aux_all_expr(ind_gene,ind_g2),num_perm,'conservative');
            all_p_values_groups(i_gene,i_g)=pval_gene;
            
        end
        
    end
    
end

%%%
all_p_values_eachgroup_days=ones(num_genes_test,length(day_comp)-1,2);
for i_g=1:(length(day_comp)-1)
    i_g
    aux_all_expr1=all_days_expr{1,i_g}.log{1,3};
    aux_all_expr2=all_days_expr{1,i_g+1}.log{1,3};
    genes1=all_days_genes{1,i_g};
    genes2=all_days_genes{1,i_g+1};
    
    [list_of_genes,ind_inters1,ind_inters2]=intersect(genes1,genes2);

    aux_all_expr1=aux_all_expr1(ind_inters1,:);
    aux_all_expr2=aux_all_expr2(ind_inters2,:);
    labels1=all_days_labels{1,i_g};
    labels2=all_days_labels{1,i_g+1};

    for i_gene=1:num_genes_test
        
        gene_p=struct_genes{i_gene,1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        if ~isempty(ind_gene)
            
            pval_gene1=permtest(aux_all_expr1(ind_gene,(labels1==0)),aux_all_expr2(ind_gene,(labels2==0)),num_perm,'conservative');
            pval_gene2=permtest(aux_all_expr1(ind_gene,(labels1==1)),aux_all_expr2(ind_gene,(labels2==1)),num_perm,'conservative');
            all_p_values_eachgroup_days(i_gene,i_g,1)=pval_gene1;
            all_p_values_eachgroup_days(i_gene,i_g,2)=pval_gene2;
            
        end
        
    end
    
end

stop=1;
elseif 1
   load test_hist.mat 
end


threshold_p_corr=threshold_p/num_genes_test;

all_p_values_days_val=sum((all_p_values_days<threshold_p_corr),2);
all_p_values_groups_val=sum((all_p_values_groups<threshold_p_corr),2);
all_p_values_eachgroup_days_val=sum(sum((all_p_values_eachgroup_days<threshold_p_corr),2),3);

% Plot histograms

%%%%%
if 1
ind_plot=find(all_p_values_days_val>=min_tests_rej);
for i_gene=1:length(ind_plot)
    i_g
    FigHandle = figure('Position', [100, 100, 1400, 600]);    
    for i_exp=1:(length(day_comp)-1)
        subplot(1,(length(day_comp)-1),i_exp);
        aux_all_expr1=all_days_expr{1,i_exp}.log{1,3};
        aux_all_expr2=all_days_expr{1,i_exp+1}.log{1,3};
        genes1=all_days_genes{1,i_exp};
        genes2=all_days_genes{1,i_exp+1};
        [list_of_genes,ind_inters1,ind_inters2]=intersect(genes1,genes2);
        aux_all_expr1=aux_all_expr1(ind_inters1,:);
        aux_all_expr2=aux_all_expr2(ind_inters2,:);
        gene_p=struct_genes{ind_plot(i_gene),1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        if ~isempty(ind_gene)
            histogram(aux_all_expr1(ind_gene,:),'Normalization','probability');
            hold on
            histogram(aux_all_expr2(ind_gene,:),'Normalization','probability');
            list_file_genes=struct_genes{ind_plot(i_gene),2};
            list_file_genes(strfind(list_file_genes,'_'))=' ';            
            title(['Day ' num2str(day_comp(i_exp)) ' vs Day ' num2str(day_comp(i_exp+1))...
                ' Gene ' gene_p ', Pval ' num2str(all_p_values_days(ind_plot(i_gene),i_exp))...
                ' ' list_file_genes]);
            legend(['Day ' num2str(day_comp(i_exp))],['Day ' num2str(day_comp(i_exp+1))])
        end
        
    end

end
end

%%%%%
if 1
ind_plot=find(all_p_values_groups_val>=min_tests_rej);
for i_gene=1:length(ind_plot)
    i_g
    FigHandle = figure('Position', [100, 100, 1400, 600]);        
    for i_exp=1:length(day_comp)
        subplot(1,length(day_comp),i_exp);
        aux_all_expr=all_days_expr{1,i_exp}.log{1,3};
        list_of_genes=all_days_genes{1,i_exp};
        labels=all_days_labels{1,i_exp};
        ind_g1=find(labels==0);
        ind_g2=find(labels==1);
                
        gene_p=struct_genes{ind_plot(i_gene),1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        if ~isempty(ind_gene)
            histogram(aux_all_expr(ind_gene,ind_g1),'Normalization','probability');
            hold on
            histogram(aux_all_expr(ind_gene,ind_g2),'Normalization','probability');
            list_file_genes=struct_genes{ind_plot(i_gene),2};
            list_file_genes(strfind(list_file_genes,'_'))=' ';
            title(['MutantVSControl Day ' num2str(day_comp(i_exp))...
                ' Gene ' gene_p ', Pval ' num2str(all_p_values_groups(ind_plot(i_gene),i_exp))...
                ' ' list_file_genes]);
            legend('Mutant','Control')
        end
        
    end

end
end

%%%%%
if 1
ind_plot=find(all_p_values_eachgroup_days_val>=2*min_tests_rej);
for i_gene=1:length(ind_plot)
    i_g
    FigHandle = figure('Position', [100, 100, 1400, 1000]);        
    for i_exp=1:(length(day_comp)-1)
        
        aux_all_expr1=all_days_expr{1,i_exp}.log{1,3};
        aux_all_expr2=all_days_expr{1,i_exp+1}.log{1,3};
        genes1=all_days_genes{1,i_exp};
        genes2=all_days_genes{1,i_exp+1};
        [list_of_genes,ind_inters1,ind_inters2]=intersect(genes1,genes2);
        aux_all_expr1=aux_all_expr1(ind_inters1,:);
        aux_all_expr2=aux_all_expr2(ind_inters2,:);
        labels1=all_days_labels{1,i_exp};
        labels2=all_days_labels{1,i_exp+1};
        
        gene_p=struct_genes{ind_plot(i_gene),1};
        ind_gene=cellfun(@(s) (strcmp(s,gene_p)), list_of_genes);
        ind_gene=find(ind_gene==1);
        subplot(1,length(day_comp)-1,i_exp);
        if ~isempty(ind_gene)
            histogram(aux_all_expr1(ind_gene,(labels1==0)),'Normalization','probability');
            hold on
            histogram(aux_all_expr2(ind_gene,(labels2==0)),'Normalization','probability');
            list_file_genes=struct_genes{ind_plot(i_gene),2};
            list_file_genes(strfind(list_file_genes,'_'))=' ';            
            title(['Per group Day ' num2str(day_comp(i_exp)) 'vs Day ' num2str(day_comp(i_exp+1))...
                ' Gene ' gene_p ', Pval ' num2str(all_p_values_eachgroup_days(ind_plot(i_gene),i_exp,1))...
                ' ' list_file_genes]);
        end
        subplot(2,length(day_comp)-1,(length(day_comp)-1)+i_exp);
        if ~isempty(ind_gene)
            histogram(aux_all_expr1(ind_gene,(labels1==1)),'Normalization','probability');
            hold on
            histogram(aux_all_expr2(ind_gene,(labels2==1)),'Normalization','probability');
            list_file_genes=struct_genes{ind_plot(i_gene),2};
            list_file_genes(strfind(list_file_genes,'_'))=' ';            
            title(['Per group Day ' num2str(day_comp(i_exp)) 'vs Day ' num2str(day_comp(i_exp+1))...
                ' Gene ' gene_p ', Pval ' num2str(all_p_values_eachgroup_days(ind_plot(i_gene),i_exp,2))...
                ' ' list_file_genes]);
            legend(['Day ' num2str(day_comp(i_exp))],['Day ' num2str(day_comp(i_exp+1))])
        end
        
    end

end
end
 
