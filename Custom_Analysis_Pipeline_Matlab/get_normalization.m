
function [array_cells_norm,array_cells_sca]=get_normalization(array_cells)
    
    array_cells_norm=10000*(array_cells./kron(sum(array_cells,1),ones(size(array_cells,1),1)));
    array_cells_sca=zscore(array_cells,[],2);

end