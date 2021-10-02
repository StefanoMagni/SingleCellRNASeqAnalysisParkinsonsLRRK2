
    
number_pcs=20;

min_expr_genes=10;

num_dim_tsne=3;

folder_res='./Images_results/Replication/';

thres_pval=0.001;
max_num_genes=10;
th_prop_cells=10;
folder_pheno='./Drop_Seq/';
list_pheno={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt','list_ph_mDANeurons.txt',...
    'list_stem_2.txt'};
list_pheno_Linn={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt','list_ph_mDANeurons_Linnarson2.txt',...
    'list_stem_2.txt'};
num_perm=20000;
number_plot=1;

% Flags to do stuff
flag_plot=1;
flag_save=0;
if flag_save
    flag_plot=1;
end
flag_useTSNE=1;
flag_comp_log=1;
flag_min_pval=1;
flag_comb_cumgenes=2; % 1 for sum, 2 for median, 3 for mean
flag_perc_cells=1; % For running in pipeline_reduce_more the percentage of cells in the most highly expressed cluster
flag_comp_pheno_cell=1; % For running in pipeline_reduce_more the computation of phenotypes
flag_log_cum=1; % To use the log of the cumulative expression for computing the p value

colormap_val='gray';
alpha_val=0.2;
flag_legend=0;
flag_nothing=1;

%%
% Open the file
day_comp=0;
    
folder_in={'./LinnarsonData/'};
filename={'Human_Embryo_fulldataset_week_6.txt'};
numb_cells={'all'};

input_param.folder_in=folder_in;
input_param.filename=filename;
input_param.numb_cells=numb_cells;
input_param.min_expr=0;

% Open, getting common genes and normalizing
out_struct=open_preprocessGeneExpr(input_param);

cell_type_ord=out_struct.names_cells{1,1};
table_expr_ord=out_struct.all_expr_mat{1,1};
table_expr_ord_sca=out_struct.all_expr_mat{3,1};
list_of_genes=out_struct.names_genes{1,1};

% For knowing the cells that they claim as mDA neurons
labels=zeros(1,size(table_expr_ord,2));
name_search='Prog';
ind_celltype=[];
for i_c=1:length(labels)
    if strfind(cell_type_ord{1,i_c},name_search)
        ind_celltype=[ind_celltype,i_c];
    end;
end;
labels(ind_celltype)=1;

%% Let's apply PCA
%%%% Both groups
[~,score,~,~,explained,~]=pca(table_expr_ord_sca');
norm_eig_s=cumsum(explained);
ind_minBoth=min(find(norm_eig_s>75));
table_expr_ord_sca_trans=score';
fprintf('Variance captured = %f \n',norm_eig_s(number_pcs))

% TSNE
%%%% Both group
% If this line fails is because it is configured the wrong tsne function.
% It should be the Matlab tsne
[table_expr_ord_sca_tsne]=tsne(table_expr_ord_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
 
colors_label=[1,0,0;0,0,1];
% Plot for groups
scatter_colorCoded(table_expr_ord_sca_tsne,colors_label,labels,num_dim_tsne);


% Plot for a gene
gene_colorcode='TH';
ind_gene=cellfun(@(s) (strcmp(s,gene_colorcode)), list_of_genes);
ind_gene=find(ind_gene==1);
colorsA=jet(max(table_expr_ord(ind_gene,:))+1);
scatter_colorCoded(table_expr_ord_sca_tsne,colorsA,table_expr_ord(ind_gene,:),num_dim_tsne);
colorbar
colormap('jet')



% Obtain phenotypical expression per cell
if flag_comp_pheno_cell
    
    all_expr{1,1}=table_expr_ord;
    all_expr_log{1,1}=log10(table_expr_ord+1);
    all_assignments=ones(1,size(table_expr_ord,2));
    data_group=table_expr_ord_sca_tsne;
    az=-15.9000;
    el=61.2000;
    
    title_eachfig=['Phenotypes with cumulative expressions for each cell, Day ' num2str(day_comp)];
    name_save=[folder_res 'DiffPheno_PerCell_CumExpr_Day' num2str(day_comp) '.pdf'];
    flag_comb_cumgenes=3;
    colors_pheno=hsv(length(list_pheno));
    colors_pheno=cat(1,[0,0,0],colors_pheno);
    [table_expr_all_phenos,all_genes_to_check,index_group_genes,ind_max_pheno]=...
        func_assign_phenotypes_percell_noClustering(list_pheno,folder_pheno,list_of_genes,all_expr,all_expr_log,...
        flag_comb_cumgenes,flag_plot,num_dim_tsne,data_group,...
        colors_pheno,title_eachfig,az,el,flag_save,name_save,labels,...
        alpha_val,colormap_val,flag_log_cum,day_comp,flag_nothing);
    title('Our list of genes, and Linnarson data for Human embryo 6 weeks')
    [table_expr_all_phenos_Linn,all_genes_to_check_Linn,index_group_genes_Linn,ind_max_pheno_Linn]=...
        func_assign_phenotypes_percell_noClustering(list_pheno_Linn,folder_pheno,list_of_genes,all_expr,all_expr_log,...
        flag_comb_cumgenes,flag_plot,num_dim_tsne,data_group,...
        colors_pheno,title_eachfig,az,el,flag_save,name_save,labels,...
        alpha_val,colormap_val,flag_log_cum,day_comp,flag_nothing);    
    title('Our list of genes, except for mDA list of genes, and Linnarson data for Human embryo 6 weeks')
    
    % Now, we obtain a cumulative expression for the genes indicated in the
% Linnarson paper, and then we plot it, also depicting the cells that we have marked 
% as the ones corresponding to mDA neurons
    
pheno_sel=8;
labels_aux=labels;
data_venn=[sum(labels_aux),sum(labels_aux.*(ind_max_pheno==pheno_sel)),...
    sum(ind_max_pheno==pheno_sel),sum((ind_max_pheno_Linn==pheno_sel).*(ind_max_pheno==pheno_sel)),...
    sum(ind_max_pheno_Linn==pheno_sel),sum((ind_max_pheno_Linn==pheno_sel).*labels_aux),...
    sum(labels_aux.*(ind_max_pheno==pheno_sel).*(ind_max_pheno_Linn==pheno_sel))];
fprintf('Correctly detected with our list: %f - FP: %f which is %f of total labeled \n',...
    sum(labels_aux.*(ind_max_pheno==pheno_sel))/sum(labels_aux),...
    sum(not(labels_aux).*(ind_max_pheno==pheno_sel)),...
    sum(not(labels_aux).*(ind_max_pheno==pheno_sel))/sum(ind_max_pheno==pheno_sel));
fprintf('Correctly detected with Linnarson list: %f - FP: %f which is %f of total labeled \n',...
    sum(labels_aux.*(ind_max_pheno_Linn==pheno_sel))/sum(labels_aux),...
    sum(not(labels_aux).*(ind_max_pheno_Linn==pheno_sel)),...
    sum(not(labels_aux).*(ind_max_pheno_Linn==pheno_sel))/sum(ind_max_pheno_Linn==pheno_sel));

figure;vennX(data_venn,1/10)

end