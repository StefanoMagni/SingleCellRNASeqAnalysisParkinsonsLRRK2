
function scatter_colorCoded(array_plot,colors,array_max,num_dim_tsne)

if num_dim_tsne==3
    scatter3(array_plot(:,1),array_plot(:,2),array_plot(:,3),30,colors(array_max+1,:),'filled')
elseif num_dim_tsne==2
    scatter(array_plot(:,1),array_plot(:,2),30,colors(array_max+1,:),'filled');
end
