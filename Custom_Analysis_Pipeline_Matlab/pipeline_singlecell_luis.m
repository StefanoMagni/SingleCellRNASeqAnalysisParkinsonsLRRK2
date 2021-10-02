
% All the functions for the pipeline of the single cell

function pipeline_singlecell_luis(folder_files,day_comp,day_comp_index,...
    name_groups,numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,...
    flag_save,flag_useTSNE,flag_comp_log,folder_res,max_clusters,all_genes_to_explore,...
    thres_pval,max_num_genes)

num_perm=10000;

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
        figure
        scatter_colorCoded(A_expr_sca_tsne,colorsA,A_expr_max(ind_gene,:),num_dim_tsne)
        colorbar
        colormap('jet')
    end
    if flag_save
        print([folder_res 'Mutant_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
            num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
    end
    
    
    if flag_plot
        figure
        scatter_colorCoded(B_expr_sca_tsne,colorsB,B_expr_max(ind_gene,:),num_dim_tsne)
        colorbar
        colormap('jet')
    end
    if flag_save
        print([folder_res 'Control_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
            num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
    end
    
    if flag_plot
        figure
        scatter_colorCoded(Both_expr_sca_tsne,colorsBoth,...
            cat(2,A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)),num_dim_tsne);
        colorbar
        colormap('jet')
    end
    % Color code by label
    %scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
    if flag_save
        print([folder_res 'Both_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
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


each_num_comp=zeros(2,3);
all_structures_gm=cell(max_clusters,3);
for i_g=1:3
    data_group=all_expr_sca_cluster{i_g};
    AIC=zeros(1,max_clusters);
    for k=1:max_clusters
        all_structures_gm{k,i_g}= fitgmdist(data_group,k);
        AIC(k)= all_structures_gm{k,i_g}.AIC;
    end
    [minAIC,numComponents]=min(AIC);
    each_num_comp(:,i_g)=[minAIC;numComponents];
end

%% Get the assignment to each cluster

all_assignments=cell(1,3);
for i_g=1:3
    data_group=all_expr_sca_cluster{i_g};
    obj=all_structures_gm{each_num_comp(2,i_g),i_g};
    idx=cluster(obj,data_group);
    all_assignments{1,i_g}=idx;
end

%% Plot the assignments

if flag_plot
    
    for i_g=1:3
        data_group=all_expr_sca_cluster{i_g};
        indexes=all_assignments{1,i_g};
        obj=all_structures_gm{each_num_comp(2,i_g),i_g};
        figure
        h1 = gscatter(data_group(:,1),data_group(:,2),indexes);
        plot_gaussians(data_group,obj,indexes,h1)
        if(i_g==1); group_val='mutant'; elseif(i_g==2); group_val='control'; elseif(i_g==3); group_val='both'; end
        title(['Cluster with Gaussian to 99% for ' group_val])
        if flag_save
            print([folder_res 'Clusters_Day' num2str(day_comp) '_' group_val '_k' num2str(each_num_comp(1,i_g)) '_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
        end
    end
end

%% Now, using these cluster, explore all possibilities of comparison between clusters
% For this, we will use non-parametric testing, and we will proceed gene by
% gene. For those most significative genes below a threshold, we will then 
else
    load('alldata_test.mat')
end

all_pvalues_and_comp=cell(2,3);
indexes_all_genes=zeros(1,size(all_genes_to_explore,1));

for i_g=1:3
    
    number_comb=0;
    for i=1:floor(each_num_comp(2,i_g)/2)
        if i==each_num_comp(2,i_g)/2
            number_comb=number_comb+(nchoosek(each_num_comp(2,i_g),i)/2);
        else
            number_comb=number_comb+nchoosek(each_num_comp(2,i_g),i);
        end
    end
%     all_perms=perms(1:each_num_comp(2,i_g));
%     all_perms=all_perms(end:-1:1,:);
    
    table_pvalues=zeros(number_comb,size(all_genes_to_explore,1));
    table_comparison_and_ind=cell(number_comb,3);
    
    if flag_comp_log
        data_group=all_expr_log{1,i_g};
    else
        data_group=all_expr_sca_cluster{1,i_g};
    end
    indexes_assign=all_assignments{1,i_g};
    
    for i_gene=1:size(all_genes_to_explore,1)
       
       all_list_gene=all_genes_to_explore{i_gene};
       ind_gene=cellfun(@(s) (strcmp(s,all_list_gene)), list_of_genes);
       ind_gene=find(ind_gene==1);       
       
       if ~isempty(ind_gene)
           indexes_all_genes(i_gene)=ind_gene;
           
           count_save=1;
           for i_comb1=1:floor(each_num_comp(2,i_g)/2)
               
               comp_to_do=combntns(1:each_num_comp(2,i_g),i_comb1);
               if i_comb1==floor(each_num_comp(2,i_g)/2)
                  comp_to_do=comp_to_do(1:(size(comp_to_do,1)/2),:); 
               end
               comp_to_do=cat(2,comp_to_do,zeros(size(comp_to_do,1),each_num_comp(2,i_g)-size(comp_to_do,2)));
               for i_aux=1:size(comp_to_do,1)
                   comp_to_do(i_aux,(i_comb1+1):each_num_comp(2,i_g))=setdiff(1:each_num_comp(2,i_g),...
                       comp_to_do(i_aux,1:i_comb1));
                   
                   table_comparison_and_ind{count_save,1}=comp_to_do(i_aux,1:i_comb1);
                   table_comparison_and_ind{count_save,2}=comp_to_do(i_aux,(i_comb1+1):end);
                   
                   max_ind=max(indexes_assign);
                   aux_indexes_assign=changem(indexes_assign,kron(max_ind+1,ones(1,i_comb1)),comp_to_do(i_aux,1:i_comb1));
                   aux_indexes_assign=changem(aux_indexes_assign,kron(max_ind+2,ones(1,each_num_comp(2,i_g)-i_comb1)),...
                       comp_to_do(i_aux,(i_comb1+1):end));
                   aux_indexes_assign(aux_indexes_assign==(max_ind+1))=1;
                   aux_indexes_assign(aux_indexes_assign==(max_ind+2))=2;
                   
                   table_comparison_and_ind{count_save,3}=aux_indexes_assign;
                   
                   vec_comp1=data_group(ind_gene,(aux_indexes_assign==1));
                   vec_comp2=data_group(ind_gene,(aux_indexes_assign==2));
                   
                   %table_pvalues(count_save,i_gene)=ranksum(vec_comp1,vec_comp2);
                   table_pvalues(count_save,i_gene)=permtest(vec_comp1,vec_comp2,num_perm,'conservative');
                   
                   count_save=count_save+1;
               end
           end
           
       end
       
    end
    
    all_pvalues_and_comp{1,i_g}=table_pvalues;
    all_pvalues_and_comp{2,i_g}=table_comparison_and_ind;
    
end

%% Now, we just choose those cases with p_values below some threshold and plot them

if flag_plot

for i_g=1:3
   
    table_pvalues=all_pvalues_and_comp{1,i_g};
    table_comparison_and_ind=all_pvalues_and_comp{2,i_g};
    
    [min_val,indexes_min]=min(table_pvalues);
    ind_cols=find(((min_val<thres_pval)==1));
    row_to_plot=indexes_min(ind_cols);
    col_to_plot=ind_cols;
    
%     [p_val_asc,ind_asc]=sort(table_pvalues(:),'ascend');
%     ind_thres=max(find((p_val_asc<thres_pval)==1));
%     number_plot=min(ind_thres,max_num_genes);
%     ind_to_plot=ind_asc(1:number_plot);
%     [row_to_plot,col_to_plot]=ind2sub(size(table_pvalues),ind_to_plot);        
    
    for i_p=1:length(row_to_plot)
       
        p_val_plot=table_pvalues(row_to_plot(i_p),col_to_plot(i_p));
        ind_gene=indexes_all_genes(col_to_plot(i_p));
        indexes_comp=table_comparison_and_ind{row_to_plot(i_p),3};
        data_group_expr=all_expr{1,i_g};
        data_group=all_expr_sca_cluster{1,i_g};
        data_group_expr_log=all_expr_log{1,i_g};
        
        if ind_gene>0
            %figure;
            FigHandle = figure('Position', [100, 100, 1200, 900]);
            subplot(2,2,1)
            gscatter(data_group(:,1),data_group(:,2),indexes_comp);
            subplot(2,2,2)
            colorsA=jet(max(data_group_expr(ind_gene,:))+1);
            scatter_colorCoded(data_group,colorsA,data_group_expr(ind_gene,:),2)
            if(i_g==1); group_val='mutant'; elseif(i_g==2); group_val='control'; elseif(i_g==3); group_val='both'; end
            title(['Comparison for gene ' all_genes_to_explore{col_to_plot(i_p)} ' for group ' group_val ...
                ', p value of ' num2str(p_val_plot)])
            colorbar;
            colormap('jet');
            aux_log_gene=data_group_expr_log(ind_gene,:);
            subplot(2,2,3)
            histogram(aux_log_gene(indexes_comp==1),20,'Normalization','probability');
            hold on;
            histogram(aux_log_gene(indexes_comp==2),20,'Normalization','probability');
            subplot(2,2,4)
            boxplot(aux_log_gene,indexes_comp);
            if flag_save
                print([folder_res 'Comp_Day' num2str(day_comp) '_' group_val '_Gene' all_genes_to_explore{col_to_plot(i_p)}...
                    '_Comp' row_to_plot(i_p) 'Pval' num2str(p_val_plot) '.png'], '-dpng', '-r300')
            end
        end
        
    end
    close all
end

end



%%%%%%% STILL SOMETHING FOR OVERLAPPING THE CONTROLS VS MUTANTS!!!


    function plot_gaussians(X0,gmfit,index,h1)
        
%         mahalDist = mahal(gmfit,X0);
%         threshold=sqrt(chi2inv(0.99,2));
%         k=max(unique(index));
%         for m = 1:k
%             idx = mahalDist(:,m)<=threshold;
%             Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
%             h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
%             uistack(h2,'bottom');
%             hold on
%         end
        maxminXY=get(gca,{'XLim','YLim'});
        xlim=maxminXY{1};
        ylim=maxminXY{2};
        %fcontour(@(x1,x2)pdf(gmfit,[x1 x2]),[xlim(1) ylim(1) xlim(2) ylim(2)])
        hold on;
        ezcontour(@(x1,x2)pdf(gmfit,[x1 x2]),xlim,ylim)
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        
    end

end

