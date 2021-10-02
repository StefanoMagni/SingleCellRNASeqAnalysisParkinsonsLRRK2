
% Just a function to iteratively run the pipeline

folder_files='./jonas/';
day_comp=[0,10,14,42]; % [0,10,14,42]
day_comp_index=1:3; % 1:4
day_comp_index=3;
day_comp=day_comp(day_comp_index);
name_groups={'A','B'}; % A is mutated and B control
% The columns are A and B, and the rows 0, 10, 14 or 42 days
numb_cells=[500,500;250,250;505,350;400,300];
list_genes_plot={'SOX2','MAP2'};
number_pcs=40;

num_dim_tsne=3;

folder_res='./Images_results/Exp1/';
max_clusters=6;

file_genes='./Drop_Seq/list_stem.txt'; % list_stem / list_CellCycle / list_total
all_genes_to_explore=write_struct_fromtxt(file_genes);
empty_cols=sum(not(cellfun(@isempty,all_genes_to_explore)));
all_genes_to_explore=all_genes_to_explore(:,find(empty_cols>0));
thres_pval=0.001;
max_num_genes=10;
th_prop_cells=10;
folder_pheno='./Drop_Seq/';
list_pheno={'list_ph_Oligodendrocyte.txt','list_ph_Neuron.txt',...
    'list_ph_OPC.txt','list_ph_Astrocyte.txt',...
    'list_ph_Microglia.txt','list_ph_Endothelial.txt'};
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
flag_comb_cumgenes=1; % 1 for sum, 2 for median, 3 for mean
flag_perc_cells=1; % For running in pipeline_reduce_more the percentage of cells in the most highly expressed cluster
flag_comp_pheno=1; % For running in pipeline_reduce_more the computation of phenotypes
flag_log_cum=1; % To use the log of the cumulative expression for computing the p value


for i_d=1:length(day_comp)
    
    %         pipeline_singlecell_luis(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
    %             numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
    %             flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,max_num_genes);
    pipeline_bothgroups_diffexp(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
        numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
        flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,th_prop_cells,...
        flag_min_pval,flag_comb_cumgenes,num_perm,number_plot);
%     pipeline_reduce_more(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
%         numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
%         flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,th_prop_cells,...
%         flag_min_pval,flag_comb_cumgenes,num_perm,number_plot,list_pheno,...
%         folder_pheno,flag_comp_pheno,flag_perc_cells,flag_log_cum)
end