

function [matrix_ord,names_order,index_order,labels_order]=order_GeneExprMat(matrix,vector_to_order)

% In this case, we just receive as input a vector of either strings or
% indexes, and we reorder the matrix according to them

unique_el=unique(vector_to_order);
names_order=unique_el;

index_order=[];
labels_order=[];
clab=1;
for i=unique_el
    if iscell(unique_el)
        ind_g=find(cellfun(@(s) (strcmp(s,i)), vector_to_order));
    else
        ind_g=find(vector_to_order==i);
    end;
    index_order=[index_order,ind_g];
    labels_order=[labels_order,clab*ones(1,length(ind_g))];
    clab=clab+1;
end;

matrix_ord=matrix(:,index_order);