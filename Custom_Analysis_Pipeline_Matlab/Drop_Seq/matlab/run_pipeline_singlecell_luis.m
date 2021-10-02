
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

file_genes='./Drop_Seq/list_CellCycle.txt'; % list_stem / list_CellCycle / list_total
all_genes_to_explore=write_struct_fromtxt(file_genes);
empty_cols=sum(not(cellfun(@isempty,all_genes_to_explore)));
all_genes_to_explore=all_genes_to_explore(:,find(empty_cols>0));
thres_pval=0.001;
max_num_genes=10;
th_prop_cells=10;

% Flags to do stuff
flag_plot=1;
flag_save=0;
if flag_save
    flag_plot=1;
end
flag_useTSNE=1;
flag_comp_log=1;
flag_min_pval=1;


for i_d=1:length(day_comp)
    
%         pipeline_singlecell_luis(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
%             numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
%             flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,max_num_genes);
    pipeline_bothgroups_diffexp(folder_files,day_comp(i_d),day_comp_index(i_d),name_groups,...
        numb_cells,list_genes_plot,number_pcs,num_dim_tsne,flag_plot,flag_save,flag_useTSNE,...
        flag_comp_log,folder_res,max_clusters,all_genes_to_explore,thres_pval,th_prop_cells,...
        flag_min_pval);
end