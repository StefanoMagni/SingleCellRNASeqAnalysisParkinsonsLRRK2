

function plot_phenotypes(data_group,colors,indexes,title_eachfig,list_pheno,...
    az,el,obj,num_dim_tsne,flag_save,pval_and_cum,thres_pval,...
    all_lims,name_save,labels,alpha_val,colormap_val,flag_nothing)

flag_legend=0;

pval_and_cum_aux=pval_and_cum(:,:,1);
pval_and_cum_aux(pval_and_cum(:,:,3)>pval_and_cum(:,:,2))=1;

% Generate the figure
h1=figure('Position', [100, 100, 800, 500],'PaperOrientation','landscape');

% To obtain a first plot of single points to create the legend
for i_leg=1:size(colors,1)
    ind_cl=find(indexes==i_leg);
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
plot_gaussians(data_group,obj,indexes,h1,num_dim_tsne,alpha_val,colormap_val);
% Plot the figure for the cumulative effect
% HERE PERHAPS WE CAN CHOOSE THE ONE WITH P VALUE THAT SATISFIES,
% AND WITH THE HIGHEST CUMULATIVE EXPRESSION
if flag_legend
    leg_ob=legend(leg_cluster);
end
for i_cl=1:size(pval_and_cum_aux,2)
    ind_el_cluster=find(indexes==i_cl);
    %aux_val_cl=prob_cluster_expressed(:,i_cl);
    ind_pval_sat=find(pval_and_cum_aux(:,i_cl,1)<...
        (thres_pval/size(pval_and_cum_aux,1)));
    %[min_pval,ind_pval]=min(pval_and_cum_cluster_expressed_aux(:,i_cl));
    %if min_pval>(thres_pval/size(pval_and_cum_cluster_expressed_aux,1))
    if isempty(ind_pval_sat)
        group_belong=1; % Doesn't belong to any group
    else
        val_cum_exp=pval_and_cum(ind_pval_sat,i_cl,2);
        [max_val,ind_maxval]=max(val_cum_exp);
        ind_pval=ind_pval_sat(ind_maxval);
        group_belong=ind_pval+1;
    end
    aux_point=data_group(ind_el_cluster,:);
    if isempty(labels)
        scatter3(aux_point(:,1),aux_point(:,2),aux_point(:,3),30,colors(group_belong,:),'filled')
    else
        ind_cl0=ind_el_cluster((labels(ind_el_cluster)==0));
        ind_cl1=ind_el_cluster((labels(ind_el_cluster)==1));
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
xlim(all_lims(1,:))
ylim(all_lims(2,:))
zlim(all_lims(3,:))
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