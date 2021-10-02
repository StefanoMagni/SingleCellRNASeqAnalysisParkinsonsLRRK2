

if 1
    
folder_files='./LinnarsonData/';
filename='iPSC_fulldataset.txt';

numb_cells=[500,500;250,250;505,350;400,300];
list_genes_plot={'SOX2','MAP2'};
number_pcs=20;

min_expr_genes=10;

num_dim_tsne=3;

folder_res='./Images_results/Replication/';
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

table=readtable([folder_files filename]);
table_cell=table2cell(table);

table_expr=str2double(table_cell(6:end,3:end));

end;

list_of_genes=table_cell(6:end,1);
cell_type=table_cell(3,3:end);
day_expr=table_cell(4,3:end);

% Ordered by the level of expression, and then select the number indicated
% in numb_cells
numb_cells=size(table_expr,2);

aux_exp=sum(table_expr,1);
[val_ord,ind_ord]=sort(aux_exp,'descend');

table_expr_ord=table_expr(:,ind_ord(1:numb_cells));
cell_type_ord=cell_type(ind_ord(1:numb_cells));
day_expr_ord=day_expr(ind_ord(1:numb_cells));

% Also, reduce the tables by only considering those genes expressed
genes_exp=sum(table_expr_ord,2);
ind_genes_exp=find(genes_exp>min_expr_genes);
list_of_genes=list_of_genes(ind_genes_exp);
table_expr_ord=table_expr_ord(ind_genes_exp,:);


% Then, we perform a normalization to obtain for each cell the number
% of transcripts per 10.000 (in python P.append(10000*M[:,i]/sum(M[:,i]))),
% and also obtain the zscore, i.e., normalize to have mean 0 and variance 1
[table_expr_ord_norm,table_expr_ord_sca]=get_normalization(table_expr_ord);

day_comp='42'; % 42 or 63
ind_day=find(cellfun(@(s) (strcmp(s,['day_' day_comp])), day_expr_ord));
%labels=zeros(1,size(table_expr_ord,2));
%labels(ind_day)=1;

% For knowing the cells that they claim as mDA neurons
labels=zeros(1,size(table_expr_ord,2));
name_search='DA';
ind_DA=[];
for i_c=1:length(labels)
    if strfind(cell_type_ord{1,i_c},name_search)
        ind_DA=[ind_DA,i_c];
    end;
end;
labels(ind_DA)=1;


table_expr_ord_aux=table_expr_ord(:,ind_day);
table_expr_ord_sca_aux=table_expr_ord_sca(:,ind_day);
labels_aux=labels(ind_day);

%% Let's apply PCA
%%%% Both groups
[~,score,~,~,explained,~]=pca(table_expr_ord_sca_aux');
norm_eig_s=cumsum(explained);
ind_minBoth=min(find(norm_eig_s>75));
table_expr_ord_sca_trans=score';
fprintf('Variance captured = %f \n',norm_eig_s(number_pcs))

% TSNE
%%%% Both group
[table_expr_ord_sca_tsne]=tsne(table_expr_ord_sca_trans(1:number_pcs,:)','NumDimensions',num_dim_tsne);
 

colors_label=[1,0,0;0,0,1];
% Plot for groups
scatter_colorCoded(table_expr_ord_sca_tsne,colors_label,labels_aux,num_dim_tsne);


% Plot for a gene
gene_colorcode='TH';
ind_gene=cellfun(@(s) (strcmp(s,gene_colorcode)), list_of_genes);
ind_gene=find(ind_gene==1);
colorsA=jet(max(table_expr_ord_aux(ind_gene,:))+1);
scatter_colorCoded(table_expr_ord_sca_tsne,colorsA,table_expr_ord_aux(ind_gene,:),num_dim_tsne);
colorbar
colormap('jet')



% Obtain phenotypical expression per cell
if flag_comp_pheno_cell
    
    all_expr{1,1}=table_expr_ord_aux;
    all_expr_log{1,1}=log10(table_expr_ord_aux+1);
    all_assignments=ones(1,size(table_expr_ord_aux,2));
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
        colors_pheno,title_eachfig,az,el,flag_save,name_save,labels_aux,...
        alpha_val,colormap_val,flag_log_cum,day_comp,flag_nothing);
    
    [table_expr_all_phenos_Linn,all_genes_to_check_Linn,index_group_genes_Linn,ind_max_pheno_Linn]=...
        func_assign_phenotypes_percell_noClustering(list_pheno_Linn,folder_pheno,list_of_genes,all_expr,all_expr_log,...
        flag_comb_cumgenes,flag_plot,num_dim_tsne,data_group,...
        colors_pheno,title_eachfig,az,el,flag_save,name_save,labels_aux,...
        alpha_val,colormap_val,flag_log_cum,day_comp,flag_nothing);    
  
    % Now, we obtain a cumulative expression for the genes indicated in the
% Linnarson paper, and then we plot it, also depicting the cells that we have marked 
% as the ones corresponding to mDA neurons
    
end



if 0 
    
    % Here we are using the Human embryo dataset, in order to correlate it
    % with the data of Day 0, to see to which group we are more correlated
    filename_embryo='Human_Embryo_fulldataset.txt';
    
    table_embryo=readtable([folder_files filename_embryo]);
    table_embryo_cell=table2cell(table_embryo);
    
    table_embryo_expr=str2double(table_embryo_cell(6:end,3:end));
    
    list_of_genes_embryo=table_embryo_cell(6:end,1);
    cell_type_embryo=table_embryo_cell(3,3:end);
    day_expr_embryo=table_embryo_cell(4,3:end);
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if 0
if flag_save_first
     
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
            figure(1)
            scatter_colorCoded(Both_expr_sca_tsne,colorsBoth,...
                cat(2,A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)),num_dim_tsne);
            colorbar
            colormap('jet')
        end
        % Color code by label
        %scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
        if flag_save
            figure(1)
            print([folder_res 'Both_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            figure(2)
            print([folder_res 'Both_Day' num2str(day_comp) '_Coloured_by_Group_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            
        end
        
        if num_dim_tsne==2
            close all
        end
        
    end
    
    if flag_plot
        figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape')
        scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
        
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        box on
        %grid off
        %axis off
        if isempty(az_vals)
            [az,el]=view;
        else
            view(az_vals(day_comp_index,:));
            [az,el]=view;
        end
        title(['Results of PCA and tSNE for Day ' num2str(day_comp) ', Az ' num2str(az) ', El ' num2str(el)])
        print([folder_res 'DimensReduction_Day' num2str(day_comp) '.pdf'], '-dpdf')
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
    
    each_num_comp=zeros(2,1);
    all_structures_gm=cell(max_clusters,1);
    i_g=3;
    
    data_group=all_expr_sca_cluster{i_g};
    AIC=zeros(1,max_clusters);
    for k=1:max_clusters
        all_structures_gm{k,1}= fitgmdist(data_group,k,'Replicates',50);
        AIC(k)= all_structures_gm{k,1}.AIC;
    end
    [minAIC,numComponents]=min(AIC);
    each_num_comp(:,1)=[minAIC;numComponents];
    figure('Position', [100, 100, 600, 400])
    plot(AIC,'LineWidth',6)
    title(['AIC for 10 clusters, Day ' num2str(day_comp)])
    print([folder_res 'AICcurve_Day' num2str(day_comp) '.pdf'], '-dpdf')
    
%     save([folder_res 'all_data_test_both_Day_' num2str(day_comp) '.mat'],'AIC','all_structures_gm',...
%         'each_num_comp','all_expr','all_expr_log','all_expr_sca_cluster','labels',...
%         'list_of_genes','az','el');
        
else
    
    load([folder_res 'all_data_test_both_Day_' num2str(day_comp) '.mat']);
    
    if flag_plot
        data_plot=all_expr_sca_cluster{1,3};
        colors_label=[1,0,0;0,0,1];
        figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape')
        scatter_colorCoded(data_plot,colors_label,labels,num_dim_tsne);
        
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        box on
        %grid off
        %axis off
        view(az,el);
        title(['Results of PCA and tSNE for Day ' num2str(day_comp) ', Az ' num2str(az) ', El ' num2str(el)])
        if flag_nothing
            grid off
            axis off
            box off
            title('')
        end
        print([folder_res 'DimensReduction_Day' num2str(day_comp) '.pdf'], '-dpdf')
    end
    
    A_expr_max=all_expr{1,1};
    B_expr_max=all_expr{1,2};
    Both_expr_sca_tsne=all_expr_sca_cluster{1,3};
    for i_g=1:length(list_genes_plot)
        
        gene_colorcode=list_genes_plot{i_g};
        
        ind_gene=cellfun(@(s) (strcmp(s,gene_colorcode)), list_of_genes);
        ind_gene=find(ind_gene==1);
        
        colorsA=jet(max(A_expr_max(ind_gene,:))+1);
        colorsB=jet(max(B_expr_max(ind_gene,:))+1);
        colorsBoth=jet(max([A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)])+1);
        
        if flag_plot
            figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape')
            scatter_colorCoded(Both_expr_sca_tsne,colorsBoth,...
                cat(2,A_expr_max(ind_gene,:),B_expr_max(ind_gene,:)),num_dim_tsne);
            colorbar
            colormap('jet')
            xticklabels(''); yticklabels('');zticklabels('');
            xticks([]);yticks([]);zticks([]);
            axis square
            box on
            %grid off
            %axis off
            view(az,el);
        end
        % Color code by label
        %scatter_colorCoded(Both_expr_sca_tsne,colors_label,labels,num_dim_tsne);
        if flag_save
            figure(1)
            print([folder_res 'Both_Day' num2str(day_comp) '_Gene' gene_colorcode '_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            figure(2)
            print([folder_res 'Both_Day' num2str(day_comp) '_Coloured_by_Group_PCAdim'...
                num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
            
        end
        
        if num_dim_tsne==2
            close all
        end
        
    end
    
    
end    

if 1
    
    if flag_choose_ncomp(1)
        numComponents=flag_choose_ncomp(2);
        each_num_comp(:,1)=[AIC(flag_choose_ncomp(2)),numComponents];
    end
    
    %% Get the assignment to each cluster
    i_g=3;
    all_assignments=cell(1,1);
    
    data_group=all_expr_sca_cluster{i_g};
    obj=all_structures_gm{each_num_comp(2,1),1};
    idx=cluster(obj,data_group);
    all_assignments{1,1}=idx;
    
    % Measure of interdistances
    vec_ind_comp=[1:each_num_comp(2,1)]';
    vec_ind_comp=cat(2,kron(vec_ind_comp,ones(each_num_comp(2,1),1)),...
        kron(ones(each_num_comp(2,1),1),vec_ind_comp));
    ind_remove=each_num_comp(2,1)*((1:each_num_comp(2,1))-1)+(1:each_num_comp(2,1));
    vec_ind_comp=vec_ind_comp(setdiff(1:size(vec_ind_comp,1),ind_remove),:);
    all_distances=zeros(1,size(vec_ind_comp,1));
    
    mean_obj=obj.mu;
    for i_distance=1:size(vec_ind_comp,1)
       
        all_distances(i_distance)=sqrt(sum((mean_obj(vec_ind_comp(i_distance,1),:)-...
            mean_obj(vec_ind_comp(i_distance,2),:)).^2));
        
    end
    fprintf('Distances for Day %d: Sum: %f / Mean: %f \n',day_comp,sum(all_distances),mean(all_distances));
    fprintf('AIC for Day %d, %f \n',day_comp,obj.AIC);
    
    %% Plot the assignments
    
    data_group=all_expr_sca_cluster{i_g};
    indexes=all_assignments{1,1};
    obj=all_structures_gm{each_num_comp(2,1),1};
    
    if flag_plot        

        colorsBoth=hsv(length(unique(indexes)));
        
        h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
        % To obtain a first plot of single points to create the legend
        leg_cluster={};
        for i_leg=1:numComponents
            ind_cl=find(indexes==i_leg);
            aux_point=data_group(ind_cl(1),:);
            scatter3(aux_point(1),aux_point(2),aux_point(3),1,colorsBoth(i_leg,:),'filled')
            hold on
            leg_cluster{1,i_leg}=['Cluster ' num2str(i_leg)];
        end
        
        plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne,alpha_val,colormap_val);
        hold on
        if num_dim_tsne==2
            h1 = gscatter(data_group(:,1),data_group(:,2),indexes);
        elseif num_dim_tsne==3
            %gscatter3_custom(data_group,indexes,h1);            
            scatter_colorCoded(data_group,colorsBoth,indexes-1,num_dim_tsne);
        end
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        box on
        view(az,el)
        xlim_val=xlim;
        ylim_val=ylim;
        zlim_val=zlim;
        all_lims=[xlim_val;ylim_val;zlim_val];
        
        if(i_g==1); group_val='mutant'; elseif(i_g==2); group_val='control'; elseif(i_g==3); group_val='both'; end
        title(['Cluster with Gaussians for Day ' num2str(day_comp)])
        
        legend(leg_cluster)
        if flag_nothing
            grid off
            axis off
            box off
            title('')
            legend off
        end
        if flag_save
        %if 0
            print([folder_res 'Clusters_3D_withGaussian_Day' num2str(day_comp) '.pdf'],'-dpdf','-fillpage')
        end
        
        % NOW, with different markers for mutant and control
        h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
              
        % To obtain a first plot of single points to create the legend
        plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne,alpha_val,colormap_val);
        hold on
        
        perc_group_per_clus=zeros(3,numComponents);
        for i_leg=1:numComponents
            ind_cl=find(indexes==i_leg);
            ind_cl0=ind_cl((labels(ind_cl)==0));
            ind_cl1=ind_cl((labels(ind_cl)==1));
            perc_group_per_clus(:,i_leg)=[length(ind_cl0),length(ind_cl1),length(ind_cl)]';
            plot3(data_group(ind_cl0,1),data_group(ind_cl0,2),...
                data_group(ind_cl0,3),'s','MarkerEdgeColor',colorsBoth(i_leg,:),...
                'MarkerFaceColor',[1,1,1],'MarkerSize',6)
            hold on
            plot3(data_group(ind_cl1,1),data_group(ind_cl1,2),...
                data_group(ind_cl1,3),'o','MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',colorsBoth(i_leg,:),'MarkerSize',6)
        end
        fprintf(['Proportions of cells per group and cluster Day ' num2str(day_comp) '\n'])
        fprintf('%f ',perc_group_per_clus(1,:)./perc_group_per_clus(3,:))
        fprintf('\n');
        fprintf('%f ',perc_group_per_clus(2,:)./perc_group_per_clus(3,:))
        fprintf('\n');
        
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        box on
        view(az,el)
        xlim(all_lims(1,:))
        ylim(all_lims(2,:))
        zlim(all_lims(3,:))
        
        if(i_g==1); group_val='mutant'; elseif(i_g==2); group_val='control'; elseif(i_g==3); group_val='both'; end
        title(['Cluster with Gaussians for Day ' num2str(day_comp) ', Square Mutand and X control'])
        
        if flag_legend
            legend(leg_cluster)
        end
        if flag_nothing
            grid off
            %axis off
            %box off
            title('')
            legend off
        end
        if flag_save
        %if 0
            print([folder_res 'Clusters_andGroups_3D_withGaussian__Day' num2str(day_comp) '.pdf'],'-dpdf','-fillpage')
        end 
        
        % Only ellipsoids
        h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
        plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne,1,'gray');
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        view(az,el)
        xlim(all_lims(1,:))
        ylim(all_lims(2,:))
        zlim(all_lims(3,:))
        grid off
        axis off
        if flag_save
            %if 0
            print([folder_res 'OnlyEllipsoids_Day' num2str(day_comp) '.pdf'],'-dpdf','-fillpage')
        end
%         if flag_save
%             print([folder_res 'Clusters_Day' num2str(day_comp) '_' group_val '_k' num2str(each_num_comp(1,1)) '_PCAdim'...
%                 num2str(number_pcs) '_tSNRdim' num2str(num_dim_tsne) '.png'], '-dpng', '-r300')
%         end
        
    end
    
    %% Obtain some statistics of how many of each group are there in each cluster
    
    percent_group_per_cluster=zeros(2,numComponents);
    number_group_per_cluster=zeros(2,numComponents);
    for i=1:max(indexes)
        aux_ind=find(indexes==i);
        percent=sum(labels(aux_ind)==0)/length(aux_ind);
        number_group_per_cluster(:,i)=[sum(labels(aux_ind)==0);sum(labels(aux_ind)==1)];
        percent_group_per_cluster(:,i)=[percent;1-percent]*100;
    end
    


%% Now, using these cluster, explore all possibilities of comparison between clusters
% For this, we will use non-parametric testing, and we will proceed gene by
% gene. For those most significative genes below a threshold, we will then

number_comb=0;
for i=1:floor(each_num_comp(2,1)/2)
    if i==each_num_comp(2,1)/2
        number_comb=number_comb+(nchoosek(each_num_comp(2,1),i)/2);
    else
        number_comb=number_comb+nchoosek(each_num_comp(2,1),i);
    end
end
%     all_perms=perms(1:each_num_comp(2,i_g));
%     all_perms=all_perms(end:-1:1,:);

end

if flag_corr_genes
    
    title_eachfig=['Phenotypes with cumulative expressions for each cell, Day ' num2str(day_comp)];
    name_save=[folder_res 'DiffPheno_PerCell_CumExpr_Day' num2str(day_comp) '.pdf'];
    flag_comb_cumgenes=2;
    colors_pheno=hsv(length(list_pheno));
    colors_pheno=cat(1,[0,0,0],colors_pheno);
    [both_cum_expression_all_corr,data_group_log_red,all_genes_to_check,index_group_genes]=...
        func_correlation_lists(list_correlation,folder_pheno,list_of_genes,all_expr_log,...
        all_assignments,flag_comb_cumgenes,flag_plot,folder_res,day_comp,flag_save);

end

if flag_comp_pheno_cell
    
    title_eachfig=['Phenotypes with cumulative expressions for each cell, Day ' num2str(day_comp)];
    name_save=[folder_res 'DiffPheno_PerCell_CumExpr_Day' num2str(day_comp) '.pdf'];
    flag_comb_cumgenes=2;
    colors_pheno=hsv(length(list_pheno));
    colors_pheno=cat(1,[0,0,0],colors_pheno);
    [both_cum_expression_all_phenos,all_genes_to_check,index_group_genes,ind_max_pheno]=...
        func_assign_phenotypes_percell(list_pheno,folder_pheno,list_of_genes,all_expr,all_expr_log,...
        all_assignments,flag_comb_cumgenes,flag_plot,num_dim_tsne,data_group,...
        colors_pheno,title_eachfig,az,el,obj,flag_save,all_lims,name_save,labels,...
        alpha_val,colormap_val,flag_log_cum,day_comp,flag_nothing);
    
end

if flag_comp_pheno
    if ~exist([folder_res 'Res_ClusterPheno_Day' num2str(day_comp) '.mat'],'file')
        %     all_phenotypes_assignment=func_assign_phenotypes(list_pheno,folder_pheno,list_of_genes,all_expr_log,...
        %         all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes,flag_comp_log,...
        %         flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
        %         flag_plot,th_prop_cells,num_dim_tsne,min_perc_exp,th_expr);
        flag_comb_cumgenes_pheno=2;
        flag_per_gene=0;
        [prob_cluster_expressed,all_genes_to_check,index_group_genes,index_group_genes_cluster,...
            pvalues_all_clusters,pval_and_cum_cluster_expressed]=...
            func_assign_phenotypes_percluster(list_pheno,folder_pheno,list_of_genes,all_expr_log,...
            all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes_pheno,flag_comp_log,...
            flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
            flag_plot,th_prop_cells,num_dim_tsne,num_perm,th_expr,flag_per_gene);
        save([folder_res 'Res_ClusterPheno_Day' num2str(day_comp) '.mat'],'prob_cluster_expressed',...
            'number_comb','percent_group_per_cluster','number_group_per_cluster',...
            'data_group','indexes','obj','all_genes_to_check','index_group_genes',...
            'index_group_genes_cluster','pvalues_all_clusters','pval_and_cum_cluster_expressed',...
            'flag_comb_cumgenes_pheno');
        
    else
        load([folder_res 'Res_ClusterPheno_Day' num2str(day_comp) '.mat']);
        colors_pheno=hsv(length(list_pheno));
        colors_pheno=cat(1,[0,0,0],colors_pheno);
        h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');
        % To obtain a first plot of single points to create the legend
        aux_point=data_group;
        for i_leg=1:size(colors_pheno,1)
            ind_cl=find(indexes==i_leg);
            aux_point=data_group(1,:);
            scatter3(aux_point(1,1),aux_point(1,2),aux_point(1,3),1,colors_pheno(i_leg,:),'filled')
            hold on
            if i_leg==1
                leg_cluster={'No particular phenotype'};
            else
                aux_list=list_pheno{i_leg-1};
                aux_list(strfind(aux_list,'_'))=' ';
                leg_cluster{1,i_leg}=aux_list;
            end
        end
        view(az,el)
        leg_ob=legend(leg_cluster);
        for i_cl=1:size(prob_cluster_expressed,2)
            ind_el_cluster=find(indexes==i_cl);
            aux_val_cl=prob_cluster_expressed(:,i_cl);
            [max_val,max_ind]=max(aux_val_cl);
            if max_val==0
                group_belong=1; % Doesn't belong to any group
            else
                group_belong=max_ind+1;
            end
            aux_point=data_group(ind_el_cluster,:);
            scatter3(aux_point(:,1),aux_point(:,2),aux_point(:,3),30,colors_pheno(group_belong,:),'filled')
        end
        xticklabels(''); yticklabels('');zticklabels('');
        xticks([]);yticks([]);zticks([]);
        axis square
        xlim(all_lims(1,:))
        ylim(all_lims(2,:))
        zlim(all_lims(3,:))
        box on
        title(['Different phenotypes, Day ' num2str(day_comp)])
        leg_ob.String=leg_ob.String(1:(length(list_pheno)+1));
        if flag_nothing
            grid off
            axis off
            box off
            title('')
            legend off
        end
        if flag_save
           print([folder_res 'DifferentPhenotypes_Day' num2str(day_comp) '.pdf'], '-dpdf','-fillpage') 
        end
        % Now, another figure with the histograms of probability
        title_eachfig='Probability of expressing for cluster ';
        name_save=[folder_res 'ProbabilitiesEachPheno_Day' num2str(day_comp) '.pdf'];
        plot_barplots(prob_cluster_expressed,colors_pheno,leg_cluster,title_eachfig,flag_save,name_save)
        
        %%%% Now, plot also for the cumulative effect
        % To one all those where the expression of the cluster is smaller
                
        title_eachfig=['Different phenotypes using the cumulative expressions, Day ' num2str(day_comp)];
        name_save=[folder_res 'DifferentPhenotypes_CumulativeExpr_Day' num2str(day_comp) '.pdf'];
        plot_phenotypes(data_group,colors_pheno,indexes,title_eachfig,list_pheno,...
            az,el,obj,num_dim_tsne,flag_save,pval_and_cum_cluster_expressed,thres_pval,...
            all_lims,name_save,[],alpha_val,colormap_val,flag_nothing);

        % For those that are also differentially expressed, plot the
        % cumulative value
        title_eachfig='Cumulative expression for cluster ';
        name_save=[folder_res 'CumulativeExprEachPheno_Day' num2str(day_comp) '.pdf'];
        pval_and_cum_cluster_expressed_aux=pval_and_cum_cluster_expressed(:,:,2);
        pval_and_cum_cluster_expressed_aux(pval_and_cum_cluster_expressed(:,:,1)>...
            (thres_pval/size(pval_and_cum_cluster_expressed_aux,1)))=0;
        pval_and_cum_cluster_expressed_aux(pval_and_cum_cluster_expressed(:,:,2)<...
            pval_and_cum_cluster_expressed(:,:,3))=0;        
        plot_barplots(pval_and_cum_cluster_expressed_aux,colors_pheno,leg_cluster,title_eachfig,flag_save,name_save)
        
        %%%% Now, again the phenotypes using the cumulative effect, but
        %%%% also now with a different marker for each group
                
        title_eachfig=['Phenotypes with cumulative expressions, Day ' num2str(day_comp) ', Sq mutant, o control'];
        name_save=[folder_res 'DifferentPheno_andGroups_CumulativeExpr_Day' num2str(day_comp) '.pdf'];
        plot_phenotypes(data_group,colors_pheno,indexes,title_eachfig,list_pheno,...
            az,el,obj,num_dim_tsne,flag_save,pval_and_cum_cluster_expressed,thres_pval,...
            all_lims,name_save,labels,alpha_val,colormap_val,flag_nothing);
        
    end
end

if flag_perc_cells
    p_vals_perc_cells_cells={};
    for i_list_g=1:size(all_genes_to_explore,2)
        flag_save=1;
        flag_comb_cumgenes=2;
        p_vals_perc_cells=func_perc_cells_cluster(all_genes_to_explore{1,i_list_g},list_of_genes,all_expr_log,...
            all_expr_sca_cluster,all_assignments,all_expr,flag_comb_cumgenes,flag_comp_log,...
            flag_log_cum,each_num_comp,number_comb,number_group_per_cluster,thres_pval,...
            flag_plot,number_plot,flag_min_pval,flag_save,labels,...
            th_prop_cells,num_dim_tsne,obj,num_perm,all_genes_to_explore{2,i_list_g},az,el,...
            folder_res,day_comp,all_lims,flag_nothing);
        p_vals_perc_cells_cells{1,i_list_g}=p_vals_perc_cells;
        flag_save=0;
    end
    save([folder_res 'Res_CellsCluster_Day' num2str(day_comp) '.mat'],'p_vals_perc_cells')
end

if flag_hist_forPheno
        
    Both_expr_max=all_expr{1,3};
    aux_gene=gene_comp_pheno{1,1};
    ind_gene=cellfun(@(s) (strcmp(s,aux_gene)), list_of_genes);
    ind_gene=find(ind_gene==1);
    Both_expr_max_gene=Both_expr_max(ind_gene,:);
    cells_chosen_pheno=find(ind_max_pheno==index_chosen_pheno);
    cells_other=setdiff(1:length(Both_expr_max_gene),cells_chosen_pheno);
    figure
    histogram(Both_expr_max_gene(cells_chosen_pheno),'Normalization','probability');
    hold on;
    histogram(Both_expr_max_gene(cells_other),'Normalization','probability');
    xlim([-1,max(Both_expr_max_gene)+2])
    legend('mDA neurons','Rest of Cells')
    title(['Histogram for gene ' aux_gene ', Pheno ' num2str(index_chosen_pheno) ', Day ' num2str(day_comp)]);
    print(['Histogram_Pheno' num2str(index_chosen_pheno) '_Gene' aux_gene '_Day' num2str(day_comp) '.pdf'],'-dpdf')
    
end
end;
