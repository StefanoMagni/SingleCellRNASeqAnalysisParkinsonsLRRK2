
if 0
    folder_files='./jonas/';
    day_comp=[0,10,14,42]; % [0,10,14,42]
    day_comp_index=1:4; % 1:4
    day_comp=day_comp(day_comp_index);
    name_groups={'A','B'}; % A is mutated and B control
    % The columns are A and B, and the rows 0, 10, 14 or 42 days
    numb_cells=[500,500;250,250;505,350;400,300];
    
    
    A_knee_plots=[];
    B_knee_plots=[];
    
    A_knee_plots_inv=[];
    B_knee_plots_inv=[];
    
    all_numbers_cells=[];
    for i_days=1:length(day_comp)
        
        files_to_compare={[name_groups{1,1} '_' num2str(day_comp(i_days)) '.txt'],...
            [name_groups{1,2} '_' num2str(day_comp(i_days)) '.txt']};
        
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
        A_expr_max=get_top_cells_matrix(A_expr,size(A_expr,2));
        A_expr_max_red=get_top_cells_matrix(A_expr,numb_cells(i_days,1));
        
        B_expr=all_arrays_comp_red{1,2};
        B_expr_max=get_top_cells_matrix(B_expr,size(B_expr,2));
        B_expr_max_red=get_top_cells_matrix(A_expr,numb_cells(i_days,1));
        
        % Also, reduce the tables by only considering those genes expressed in both
        % groups
        genes_expA=sum(A_expr_max,2);
        genes_expB=sum(B_expr_max,2);
        
        ind_genes_expA=find(genes_expA>0);
        ind_genes_expB=find(genes_expB>0);
        ind_genes_exp_tot=intersect(ind_genes_expA,ind_genes_expB);
        
        A_expr_max=A_expr_max(ind_genes_exp_tot,:);
        B_expr_max=B_expr_max(ind_genes_exp_tot,:);
        
        aux_n_cells=[size(A_expr_max,2);size(B_expr_max,2)];
        all_numbers_cells=[all_numbers_cells,aux_n_cells];
        if i_days>1
            A_knee_plots=fit_size_mat(sum(A_expr_max),A_knee_plots);
            B_knee_plots=fit_size_mat(sum(B_expr_max),B_knee_plots);
        else
            A_knee_plots=sum(A_expr_max);
            B_knee_plots=sum(B_expr_max);
        end
        
        
        % Now, for the inverted order, i.e., we first filter by genes, and then
        % by number of cells
        A_expr_red=A_expr(ind_genes_exp_tot,:);
        B_expr_red=B_expr(ind_genes_exp_tot,:);
        
        A_expr_max_red_inv=get_top_cells_matrix(A_expr_red,size(A_expr_red,2));
        B_expr_max_red_inv=get_top_cells_matrix(B_expr_red,size(B_expr_red,2));
        
        if i_days>1
            A_knee_plots_inv=fit_size_mat(sum(A_expr_max_red_inv),A_knee_plots_inv);
            B_knee_plots_inv=fit_size_mat(sum(B_expr_max_red_inv),B_knee_plots_inv);
        else
            A_knee_plots_inv=sum(A_expr_max_red_inv);
            B_knee_plots_inv=sum(B_expr_max_red_inv);
        end
        
    end
    
    A_knee_plots_cum=A_knee_plots./kron(sum(A_knee_plots,2),ones(1,size(A_knee_plots,2)));
    B_knee_plots_cum=B_knee_plots./kron(sum(B_knee_plots,2),ones(1,size(B_knee_plots,2)));
    
    A_knee_plots_cum_inv=A_knee_plots_inv./kron(sum(A_knee_plots_inv,2),ones(1,size(A_knee_plots_inv,2)));
    B_knee_plots_cum_inv=B_knee_plots_inv./kron(sum(B_knee_plots_inv,2),ones(1,size(B_knee_plots_inv,2)));
    
else
    load('All_loaded_kneeplots.mat')
    A_knee_plots_cum=cumsum(A_knee_plots_cum,2);
    B_knee_plots_cum=cumsum(B_knee_plots_cum,2);
    A_knee_plots_cum_inv=cumsum(A_knee_plots_cum_inv,2);
    B_knee_plots_cum_inv=cumsum(B_knee_plots_cum_inv,2);
end

%% Plots
colors=[0.133 0.47 0.71;1 0.5 0.06;0.27 0.64 0.2;0.84 0.15 0.157];
all_numbers_cells_prev=all_numbers_cells;
all_numbers_cells=ones(2,4)*999;

var_plot{1,1}=A_knee_plots;
var_plot{1,2}=B_knee_plots;
var_plot{1,3}=A_knee_plots_cum;
var_plot{1,4}=B_knee_plots_cum;
all_titles={'G2019S','Control','G2019S Cum Exp','G2019S Cum Exp'};

figure('Position', [100, 100, 1400, 1000],'PaperOrientation','landscape')

values_chosen_cells=zeros(8,3);
values_chosen_cells(:,1)=reshape(numb_cells,[],1);
for i_sp=1:4
    subplot(2,2,i_sp);
    hold on
    var_p=var_plot{1,i_sp};
    if mod(i_sp,2)==0
        ind_g=2;
    else
        ind_g=1;
    end
    for i_p=1:4
        plot([0 0],[0 0],'Color',colors(i_p,:));
    end
    leg=legend(['Day 0, ' num2str(numb_cells(1,ind_g))],['Day 10, ' num2str(numb_cells(2,ind_g))],...
        ['Day 14, ' num2str(numb_cells(3,ind_g))],['Day 42, ' num2str(numb_cells(4,ind_g))]);
    
    for i_p=1:4        
        val_y=var_p(i_p,numb_cells(i_p,ind_g));
        plot(1:all_numbers_cells(ind_g,i_p),var_p(i_p,1:all_numbers_cells(ind_g,i_p)),'Color',colors(i_p,:),'LineWidth',4);
        plot([0 all_numbers_cells(ind_g,i_p)],[val_y val_y],'Color',colors(i_p,:));
        plot([numb_cells(i_p,ind_g) numb_cells(i_p,ind_g)],[0 30000],'Color',colors(i_p,:));
        xlim([0,all_numbers_cells(ind_g,i_p)]);
        
        if i_sp<3
            values_chosen_cells(4*(ind_g-1)+i_p,2)=val_y;
        else
            values_chosen_cells(4*(ind_g-1)+i_p,3)=val_y;
        end
        if i_sp<3
            ylim([0,30000]);
        else
            ylim([0,1]);
        end
    end
    title(all_titles{1,i_sp});
    leg.String=leg.String(1:4);
    
end


var_plot{1,1}=A_knee_plots_inv;
var_plot{1,2}=B_knee_plots_inv;
var_plot{1,3}=A_knee_plots_cum_inv;
var_plot{1,4}=B_knee_plots_cum_inv;
all_titles={'G2019S Inv','Control Inv','G2019S Inv Cum Exp ','G2019S Inv Cum Exp '};

figure('Position', [100, 100, 1400, 1000],'PaperOrientation','landscape')

values_chosen_cells_inv=zeros(8,3);
values_chosen_cells_inv(:,1)=reshape(numb_cells,[],1);
for i_sp=1:4
    subplot(2,2,i_sp);
    hold on
    var_p=var_plot{1,i_sp};
    if mod(i_sp,2)==0
        ind_g=2;
    else
        ind_g=1;
    end
    for i_p=1:4
        plot([0 0],[0 0],'Color',colors(i_p,:));
    end
    leg=legend(['Day 0, ' num2str(numb_cells(1,ind_g))],['Day 10, ' num2str(numb_cells(2,ind_g))],...
        ['Day 14, ' num2str(numb_cells(3,ind_g))],['Day 42, ' num2str(numb_cells(4,ind_g))]);    
    
    for i_p=1:4        
        val_y=var_p(i_p,numb_cells(i_p,ind_g));
        plot(1:all_numbers_cells(ind_g,i_p),var_p(i_p,1:all_numbers_cells(ind_g,i_p)),'Color',colors(i_p,:),'LineWidth',4);
        plot([0 size(var_p,2)],[val_y val_y],'Color',colors(i_p,:));
        plot([0 size(var_p,2)],[val_y val_y],'Color',colors(i_p,:));   
        plot([numb_cells(i_p,ind_g) numb_cells(i_p,ind_g)],[0 30000],'Color',colors(i_p,:));
        xlim([0,all_numbers_cells(ind_g,i_p)]);
        
        if i_sp<3
            values_chosen_cells_inv(4*(ind_g-1)+i_p,2)=val_y;
        else
            values_chosen_cells_inv(4*(ind_g-1)+i_p,3)=val_y;
        end        
        if i_sp<3
            ylim([0,30000]);
        else
            ylim([0,1]);
        end
    end
    title(all_titles{1,i_sp});
    leg.String=leg.String(1:4);
    
end

group_day={'Mutant, Day 0, ','Mutant, Day 10, ','Mutant, Day 14, ','Mutant, Day 42, ',...
    'Control, Day 0, ','Control, Day 10, ','Control, Day 14, ','Control, Day 42, '};

for i=1:8
   
    fprintf([group_day{1,i} num2str(values_chosen_cells(i,2)) ', ' ...
        num2str(values_chosen_cells(i,3)) ' and Inv, ' ...
        num2str(values_chosen_cells_inv(i,2)) ', ' ...
        num2str(values_chosen_cells_inv(i,3))])
    fprintf('\n')

end

aux_values_chosen_cells=cat(2,values_chosen_cells(1:4,2:3),values_chosen_cells(5:8,2:3));
for i=1:4
   
    fprintf([num2str(values_chosen_cells(i,1)) ' ' group_day{1,i} num2str(aux_values_chosen_cells(i,1)) ', ' ...
        num2str(aux_values_chosen_cells(i,2)) ' // ' num2str(values_chosen_cells(i+4,1)) ' ' group_day{1,4+i} ...
        num2str(aux_values_chosen_cells(i,3)) ', ' ...
        num2str(aux_values_chosen_cells(i,4))])
    fprintf('\n')

end

%%

function X_knee_plots=fit_size_mat(X_expr_max,X_knee_plots)

if length(X_expr_max)>size(X_knee_plots,2)
    X_knee_plots=cat(2,X_knee_plots,zeros(size(X_knee_plots,1),length(X_expr_max)-size(X_knee_plots,2)));
    X_knee_plots=[X_knee_plots;X_expr_max];
else
    X_expr_max=cat(2,X_expr_max,zeros(1,size(X_knee_plots,2)-length(X_expr_max)));
    X_knee_plots=[X_knee_plots;X_expr_max];
end

end


