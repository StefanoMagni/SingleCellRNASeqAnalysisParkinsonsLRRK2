

function out_struct=open_preprocessGeneExpr(input_param)

% File to open and do some quick standard preprocessing on the gene
% expression matrices. We are doing the following processes:
% 1 - Open and storing the tables with only those genes where there is some
% expression higher than min_expr
% 2 - Getting the names of the genes and intersecting them to get a common
% set
% 3 - Getting the arrays, order cells from those more express to less, pick
% a specific number of them, and then, get normalizations.
%

% Default values
min_expr=0;
flag_common_genes=1;

%% Assigning all parameters from input structure
all_inputs=fieldnames(input_param);
for i_in=1:length(all_inputs)
    eval([all_inputs{i_in} '=input_param.' all_inputs{i_in}]);
end;

if exist('filename','var')
    
    if ~iscell(filename)
        aux_filename=filename;
        filename={};
        filename{1,1}=aux_filename;
        aux_folder_in=folder_in;
        folder_in={};
        folder_in{1,1}=aux_folder_in;
    elseif ischar(folder_in)
        folder_in_aux=folder_in;
        folder_in={};
        folder_in{1,1}=folder_in_aux;
        folder_in=repmat(folder_in,1,length(filename));
    end;
    
    all_tables_comp={};
    for i=1:length(filename)
        table=readtable([folder_in{1,i} '/' filename{1,i}],'Delimiter','tab','ReadRowNames',1,'ReadVariableNames',1);
        if isempty(table)
            table=readtable([folder_in{1,i} '/' filename{1,i}],'Delimiter',',','ReadRowNames',1,'ReadVariableNames',1);
        end;
        aux=table2array(table);
        aux_cum=sum(aux,2);    
        ind_genes=find(aux_cum>min_expr);
        
        all_tables_comp{1,i}=table(ind_genes,:);
    end;

    % Get names of genes
    names_genes={};
    names_cells={};
    for i=1:length(all_tables_comp)
        names_genes{1,i}=all_tables_comp{1,i}.Properties.RowNames;
        names_cells{1,i}=all_tables_comp{1,i}.Properties.VariableNames;
    end
    
    % Get only common genes and indexes
    if flag_common_genes
        if length(all_tables_comp)>1
            [common_genes,ind_c1,ind_c2]=intersect(names_genes{1,1},names_genes{1,2});
            ind_common_genes{1,1}=ind_c1;
            ind_common_genes{1,2}=ind_c2;
            for i_int=3:length(names_genes)
                [common_genes,ind_c1,ind_c2]=intersect(common_genes,names_genes{1,i_int});
                for i_u=1:(i_int-1)
                    aux_ind=ind_common_genes{1,i_u};
                    ind_common_genes{1,i_u}=aux_ind(ind_c1);
                end;
                ind_common_genes{1,i_int}=ind_c2;
            end;
        else
            ind_common_genes{1,1}=1:size(all_tables_comp{1,1});
            common_genes=names_genes{1,i};
        end;
    else
        for i=1:length(all_tables_comp)
            ind_common_genes{1,i}=1:size(all_tables_comp{1,i},1);
            common_genes{1,i}=names_genes{1,i};
        end
    end;
    
    % Getting the order gene expression matrices only for the number of
    % specified cells, and their normalizations
    labels=[];
    for i=1:length(all_tables_comp)
        aux_expr=table2array(all_tables_comp{1,i});
        aux_expr=aux_expr(ind_common_genes{1,i},:);
        if exist('numb_cells','var')
            if strcmp(numb_cells{1,i},'all')
                aux_expr_max=get_top_cells_matrix(aux_expr,size(aux_expr,2));
            else
                aux_expr_max=get_top_cells_matrix(aux_expr,numb_cells{1,i});
            end;
        else
            aux_expr_max=get_top_cells_matrix(aux_expr,size(aux_expr,2));
        end;
        [aux_expr_norm,aux_expr_sca]=get_normalization(aux_expr_max);
        labels=cat(2,labels,(i-1)*ones(1,size(aux_expr_norm,2)));
        
        all_expr_mat{1,i}=aux_expr_max;
        all_expr_mat{2,i}=aux_expr_norm;
        all_expr_mat{3,i}=aux_expr_sca;
    end;
    
    out_struct.folder_in=folder_in;
    out_struct.filename=filename;
    out_struct.all_expr_mat=all_expr_mat;
    out_struct.common_genes=common_genes;
    out_struct.labels=labels;
    out_struct.ind_common_genes=ind_common_genes;
    out_struct.names_genes=names_genes;
    out_struct.names_cells=names_cells;
    
else
    
    out_struct.error=1;
    
end;


