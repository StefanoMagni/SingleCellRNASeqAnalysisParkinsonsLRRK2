

addpath './CodeDef'
addpath(genpath('./GeneralFunctions'))
rmpath(genpath('./GeneralFunctions/Data_analysis/tSNE/'))
all_day_comp=[0,10,14,42];
week_embryo=6;

for i_day=1:length(all_day_comp)
    
    groups_use=[2,3];
    
    day_comp=all_day_comp(i_day);
    
    folder_in={'./Tables_Phenotypes/','./Tables_Phenotypes/','./LinnarsonData/'};
    filename={['A_' num2str(day_comp) '_Red_GenExprMat.txt'],...
        ['B_' num2str(day_comp) '_Red_GenExprMat.txt'],...
        ['Human_Embryo_fulldataset_week_' num2str(week_embryo) '.txt']};
    % 'Human_Embryo_fulldataset_week_6.txt' , iPSC_fulldataset_day_63
    
    numb_cells={'all','all','all'};
    
    clims=[0.3,0.7];
    
    input_param.folder_in=folder_in;
    input_param.filename=filename;
    input_param.numb_cells=numb_cells;
    input_param.min_expr=0;
    
    % Open, getting common genes and normalizing
    out_struct=open_preprocessGeneExpr(input_param);
    
    %%
    
    all_vector_to_order={};
    % Open the phenotypes obtained
    folder_in={'./Tables_Phenotypes/','./Tables_Phenotypes/'};
    filename={['A_' num2str(day_comp) '_Highest_Pheno.txt'],['B_' num2str(day_comp) '_Highest_Pheno.txt']};
    
    for i=1:length(filename)
        table=readtable([folder_in{1,i} '/' filename{1,i}]);
        all_vector_to_order{1,i}=table2array(table);
    end;
    
    % Remove final suffix from the cell type
    aux_vector_to_order={};
    aux_names_cells=out_struct.names_cells{1,3};
    for i=1:length(aux_names_cells)
        ind_str=strfind(aux_names_cells{1,i},'_');
        aux_vector_to_order{1,i}=aux_names_cells{1,i}(1:(ind_str(end)-1));
    end;
    all_vector_to_order{1,3}=aux_vector_to_order;
    
    %
    
    %%
    % Order matrices considering some criteria
    % First, order each matrix independently
    all_expr_mat_ord={};
    for i=1:length(all_vector_to_order)
        vector_to_order=all_vector_to_order{1,i};
        expr_mat=out_struct.all_expr_mat{1,i};
        [expr_mat_ord,names_order,index_order,labels_order]=order_GeneExprMat(expr_mat,vector_to_order);
        all_expr_mat_ord{1,i}=expr_mat_ord;
        all_expr_mat_ord{2,i}=index_order;
        %labels_type=vector_to_order(index_order);
        all_expr_mat_ord{3,i}=labels_order;
        all_expr_mat_ord{4,i}=names_order;
    end;
        
    
    %% Obtain correlation matrices using cumulative expression for all cells of one tipe
    %cov_mat=corr(zscore(all_expr_mat_ord{1,1}')',zscore(all_expr_mat_ord{1,3}')');
    
    phenos={'Oligodendrocyte','Neuron','OPC','Astrocyte','Microglia',...
        'Endothelial','mDA','Stemness'};
    
    all_expr_mat_redtype={};
    all_expr_mat_redtype_sum={};
    all_expr_mat_together=[];
    vec_label_all_expr=[];
    for i=1:size(all_expr_mat_ord,2)
        mattored=all_expr_mat_ord{1,i};
        mattored_norm=zscore(mattored);
        indextored=all_expr_mat_ord{3,i};
        mat_red=zeros(size(mattored,1),length(unique(indextored)));
        mat_red_sum=zeros(size(mattored,1),length(unique(indextored)));
        for i_un=unique(indextored)
            mat_red(:,i_un)=mean(mattored(:,(indextored==i_un)),2);
            mat_red_sum(:,i_un)=mean(mattored_norm(:,(indextored==i_un)),2);
        end;
        all_expr_mat_redtype{1,i}=mat_red;
        all_expr_mat_redtype_sum{1,i}=mat_red_sum;
    end;
    
    %cov_mat=corr(all_expr_mat_redtype{1,groups_use(1)},all_expr_mat_redtype{1,groups_use(2)});
    cov_mat=corr(all_expr_mat_redtype_sum{1,groups_use(1)},all_expr_mat_redtype_sum{1,groups_use(2)});
    
    aux_mat=[all_expr_mat_redtype{1,groups_use(1)},all_expr_mat_redtype{1,groups_use(2)}];
    dist_mat=squareform(pdist(zscore(aux_mat)','euclidean'));
    
    figure(1);
    subplot(2,2,i_day);
    imagesc(cov_mat)
    xticks(1:length(all_expr_mat_ord{4,groups_use(2)}))
    xticklabels(all_expr_mat_ord{4,groups_use(2)})
    yticks(1:length(all_expr_mat_ord{4,groups_use(1)}))
    yticklabels(phenos(all_expr_mat_ord{4,groups_use(1)}))
    title(['CumCorr Day ' num2str(day_comp) ' and Week ' num2str(week_embryo)])
    colorbar
    caxis(clims)
    
    %% In this other method, we only pick the most expressed X genes for both
    % for each group
    
    num_genes_sel=10000;
    
    corr_mat_selgenes=zeros(length(all_expr_mat_ord{4,groups_use(1)}),...
        length(all_expr_mat_ord{4,groups_use(2)}));
    mutinfo_selgenes=zeros(length(all_expr_mat_ord{4,groups_use(1)}),...
        length(all_expr_mat_ord{4,groups_use(2)}));
    
    auxmat1=all_expr_mat_redtype{1,groups_use(1)};
    auxmat2=all_expr_mat_redtype{1,groups_use(2)};
    
    auxmat1_n=all_expr_mat_redtype_sum{1,groups_use(1)};
    auxmat2_n=all_expr_mat_redtype_sum{1,groups_use(2)};
    
    for i1=1:size(auxmat1,2)
        
        for i2=1:size(auxmat2,2)
            
            all_genes_exp=abs(auxmat1(:,i1))+abs(auxmat2(:,i2));
            [val_genes,ind_max]=sort(all_genes_exp,'descend');
            
            ind_max=ind_max(1:num_genes_sel);
            corr_val=corr(auxmat1_n(ind_max,i1),auxmat2_n(ind_max,i2));
            mutinfo_val=mutualinfo(auxmat1_n(ind_max,i1),auxmat2_n(ind_max,i2));
            %corr_val=var(auxmat1(ind_max,i1),auxmat2(ind_max,i2));
            
            corr_mat_selgenes(i1,i2)=corr_val;
            mutinfo_selgenes(i1,i2)=mutinfo_val;
            
        end;
        
    end;
    
    figure(2);
    subplot(2,2,i_day);
    imagesc(corr_mat_selgenes);
    xticks(1:length(all_expr_mat_ord{4,groups_use(2)}))
    xticklabels(all_expr_mat_ord{4,groups_use(2)})
    yticks(1:length(all_expr_mat_ord{4,groups_use(1)}))
    yticklabels(phenos(all_expr_mat_ord{4,groups_use(1)}))
    title(['CumCorr SelGenes Day ' num2str(day_comp) ' and Week ' num2str(week_embryo)])
    colorbar
    caxis(clims)
end;





