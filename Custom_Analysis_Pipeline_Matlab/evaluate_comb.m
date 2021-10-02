

function [table_pvalues,table_comparison_and_ind]=evaluate_comb(ind_gene,...
    each_num_comp,indexes_assign,data_group,number_comb,flag_cum,number_group_per_cluster)

if flag_cum
    table_pvalues=zeros(number_comb,3);
else
    table_pvalues=zeros(number_comb,1);
end
table_comparison_and_ind=cell(number_comb,4);

if ~isempty(ind_gene)
    
    count_save=1;
    for i_comb1=1:floor(each_num_comp(2,1)/2)
        
        comp_to_do=combntns(1:each_num_comp(2,1),i_comb1);
        if i_comb1==floor(each_num_comp(2,1)/2)
            comp_to_do=comp_to_do(1:(size(comp_to_do,1)/2),:);
        end
        comp_to_do=cat(2,comp_to_do,zeros(size(comp_to_do,1),each_num_comp(2,1)-size(comp_to_do,2)));
        for i_aux=1:size(comp_to_do,1)
            comp_to_do(i_aux,(i_comb1+1):each_num_comp(2,1))=setdiff(1:each_num_comp(2,1),...
                comp_to_do(i_aux,1:i_comb1));
            
            table_comparison_and_ind{count_save,1}=comp_to_do(i_aux,1:i_comb1);
            table_comparison_and_ind{count_save,2}=comp_to_do(i_aux,(i_comb1+1):end);
            
            max_ind=max(indexes_assign);
            aux_indexes_assign=changem(indexes_assign,kron(max_ind+1,ones(1,i_comb1)),comp_to_do(i_aux,1:i_comb1));
            aux_indexes_assign=changem(aux_indexes_assign,kron(max_ind+2,ones(1,each_num_comp(2,1)-i_comb1)),...
                comp_to_do(i_aux,(i_comb1+1):end));
            aux_indexes_assign(aux_indexes_assign==(max_ind+1))=1;
            aux_indexes_assign(aux_indexes_assign==(max_ind+2))=2;
            
            table_comparison_and_ind{count_save,3}=aux_indexes_assign;
            
            vec_comp1=data_group(ind_gene,(aux_indexes_assign==1));
            vec_comp2=data_group(ind_gene,(aux_indexes_assign==2));
            
            if flag_cum                
                table_pvalues(count_save,1)=ranksum(data_group(1,(aux_indexes_assign==1)),...
                    data_group(1,(aux_indexes_assign==2)));
                aux_c1=table_comparison_and_ind{count_save,1};
                aux_c2=table_comparison_and_ind{count_save,2};
                table_pvalues(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2))';
                %table_pvalues_cum_exp(count_save,2:3)=(sum(number_group_per_cluster(:,aux_c1),2)./sum(number_group_per_cluster,2))';
                %table_pvalues_cum_exp(count_save,3)=sum(number_group_per_cluster(:,aux_c2),2)./sum(number_group_per_cluster,2);
            else
                table_pvalues(count_save,1)=ranksum(vec_comp1,vec_comp2);
            end
            
            count_save=count_save+1;
        end
    end
end