
% Just a function to iteratively run the pipeline

folder_files='./jonas/';
day_comp=[0,10,14,42]; % [0,10,14,42]
day_comp_index=1:4; % 1:4
%day_comp_index=2;
day_comp=day_comp(day_comp_index);
name_groups={'A','B'}; % A is mutated and B control
% The columns are A and B, and the rows 0, 10, 14 or 42 days
numb_cells=[500,500;250,250;505,350;400,300];
list_genes_plot={'SOX2','MAP2'};
list_genes_plot={};
gene_comp_pheno={'LRRK2'};
number_pcs=20;

az_vals=[-155.9,14.8;-161.1,17.2;40.9,13.2;-5.5,9.2];

num_dim_tsne=3;

%folder_res='./Images_results/Exp1/';
folder_res='./Results/Sims06_06_20dimPCA_MeanCum_Pval0c01_Paper1/';
max_clusters=10;

file_genes={'./Drop_Seq/list_stem_2.txt','./Drop_Seq/list_CellCycle.txt',...
    './Drop_Seq/list_ph_mDANeurons.txt','./Drop_Seq/list_Mito.txt'}; % list_stem / list_CellCycle / list_total
%file_genes={'./Drop_Seq/list_Mito.txt'};
all_genes_to_explore={};
for i_list_genes=1:length(file_genes)
    aux_all_genes_to_explore=write_struct_fromtxt(file_genes{1,i_list_genes});
    empty_cols=sum(not(cellfun(@isempty,aux_all_genes_to_explore)));
    all_genes_to_explore{1,i_list_genes}=aux_all_genes_to_explore(:,find(empty_cols>0));
    [~,file_name,~]=fileparts(file_genes{1,i_list_genes});
    file_name(strfind(file_name,'_'))=' ';
    all_genes_to_explore{2,i_list_genes}=file_name;
end

thres_pval=0.01;
max_num_genes=10;
th_prop_cells=10;
folder_pheno='./Drop_Seq/';
list_pheno={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt','list_ph_mDANeurons.txt',...
    'list_stem_2.txt'};
% list_pheno={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
%     'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
%     'list_ph_Microglia.txt','list_ph_Endothelial.txt','list_ph_mDANeurons_Linnarson2.txt',...
%     'list_stem_2.txt'};
list_correlation={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt','list_ph_mDANeurons.txt',...
    'list_stem_2.txt','list_CellCycle.txt','list_Mito.txt'};
num_perm=20000;
number_plot=1;
min_perc_exp=0.1;
th_expr=1;

% Flags to do stuff
flag_plot=1;
flag_save=0;
if flag_save
    flag_plot=1;
end
flag_useTSNE=1;
flag_comp_log=1;
flag_min_pval=1;
flag_comb_cumgenes=2; % 1 for sum, cl2 for median, 3 for mean
flag_perc_cells=0; % For running in pipeline_reduce_more the percentage of cells in the most highly expressed cluster
flag_comp_pheno=0; % For running in pipeline_reduce_more the computation of phenotypes
flag_comp_pheno_cell=0;
flag_log_cum=1; % To use the log of the cumulative expression for computing the p value
flag_corr_genes=0;
flag_hist_forPheno=1;
if flag_hist_forPheno
    flag_comp_pheno_cell=1;
end
index_chosen_pheno=7;

az_vals=[];
flag_save_first=0; % To run or not the first part, where we open all the data, and compute PCA,
% tSNE, and obtain the clusters and assignments
flag_choose_ncomp=[1,8];

for i_d=1:length(day_comp)

    pipeline_reduce_more(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
        numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
        flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,th_prop_cells,...
        flag_min_pval,flag_comb_cumgenes,num_perm,number_plot,list_pheno,list_correlation,...
        folder_pheno,flag_comp_pheno,flag_comp_pheno_cell,flag_perc_cells,flag_log_cum,min_perc_exp,th_expr,...
        flag_save_first,flag_choose_ncomp,flag_corr_genes,az_vals,flag_hist_forPheno,...
        index_chosen_pheno,gene_comp_pheno);

end