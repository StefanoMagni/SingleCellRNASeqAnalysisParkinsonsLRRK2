

function h1=gscatter3_custom(data,indexes,h1)

[~,~,id] = unique(indexes);
colors = 'rgbmckyrgbmckyrgbmcky';
markers = 'osdosdosdosdosdosd';
for idx = 1:max(id)
    data_aux=data(indexes == idx,:);
    plot3(data_aux(:,1), data_aux(:,2), data_aux(:,3), [colors(idx) markers(idx)]);
    hold on;
end
grid

