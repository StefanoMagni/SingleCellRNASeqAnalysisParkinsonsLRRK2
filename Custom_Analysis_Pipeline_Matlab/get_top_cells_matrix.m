

function array_red=get_top_cells_matrix(array_cells,number_cells)
    
    aux_exp=sum(array_cells,1);
    [val_ord,ind_ord]=sort(aux_exp,'descend');
    array_red=array_cells(:,ind_ord(1:number_cells)); 

end