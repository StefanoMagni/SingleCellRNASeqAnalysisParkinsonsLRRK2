
function output_names=convert_ceftotxt(folder_files,filename,input_param)

% Simply, to convert the cef files into txt, splitting by days, and
% indicating in the name of the cell its type, if provided. The cef has to
% be copied already in a txt, and the files are separated by days. The name
% for the variables are the cell types

% Default values
range_table=[6,3];
ind_type=[3,3];
ind_day=[4,3];
ind_genes=[6,1];


%% Assigning all parameters from input structure
if exist('input_param','var')
    all_inputs=fieldnames(input_param);
    for i_in=1:length(all_inputs)
        eval([all_inputs{i_in} '=input_param.' all_inputs{i_in}]);
    end;
end;

%
if ~iscell(filename)
    aux_filename=filename;
    filename={};
    filename{1,1}=aux_filename;
    aux_folder_in=folder_files;
    folder_files={};
    folder_files{1,1}=aux_folder_in;
elseif ischar(folder_files)
    folder_in_aux=folder_files;
    folder_files={};
    folder_files{1,1}=folder_in_aux;
    folder_files=repmat(folder_files,1,length(filename));
end;

% Opening file
output_names={};
cnames=1;
for i_f=1:length(filename)
    
    table=readtable([folder_files{1,i_f} '/' filename{1,i_f}]);
    filename_aux=filename{1,i_f}(1:(end-4));
    
    table_cell=table2cell(table);
    
    table_expr=str2double(table_cell(range_table(1):end,range_table(2):end));
    
    cell_type=table_cell(ind_type(1),ind_type(2):end);
    day_expr=table_cell(ind_day(1),ind_day(2):end);
    list_of_genes=table_cell(ind_genes(1):end,ind_genes(2));
    
    diff_days=unique(day_expr);
    
    for i_d=1:length(diff_days)
        
        ind_day=find(cellfun(@(s) (strcmp(s,diff_days{1,i_d})), day_expr));
        table_red=table_expr(:,ind_day);
        cell_type_red=cell_type(1,ind_day);
        
        name_file=[filename_aux '_' diff_days{1,i_d} '.txt'];
        
        cell_type_red_ext={};
        for i_name=1:length(cell_type_red)
            cell_type_red_ext{1,i_name}=[cell_type_red{1,i_name} '_' num2str(i_name)];
        end;
        table_aux=array2table(table_red,'RowNames',list_of_genes,'VariableNames',cell_type_red_ext);
        writetable(table_aux,[folder_files{1,i_f} '/' name_file],'WriteRowNames',1,'WriteVariableNames',1)
        
        output_names{1,cnames}=[folder_files{1,i_f} '/' name_file];
        cnames=cnames+1;
        
    end;
    
end;