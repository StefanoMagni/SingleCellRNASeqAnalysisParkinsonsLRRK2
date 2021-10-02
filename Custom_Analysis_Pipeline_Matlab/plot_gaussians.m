

function plot_gaussians(X0,gmfit,index,h1,num_dim_tsne,alpha_val,colormap_val)

if nargin==6
   colormap_val='copper';
end

%fcontour(@(x1,x2)pdf(gmfit,[x1 x2]),[xlim(1) ylim(1) xlim(2) ylim(2)])
if nargin==5
    alpha_val=0.1;
    colormap_val='copper';
end
hold on;
if num_dim_tsne==2
    maxminXY=get(gca,{'XLim','YLim'});
    xlim=maxminXY{1};
    ylim=maxminXY{2};
    ezcontour(@(x1,x2)pdf(gmfit,[x1 x2]),xlim,ylim)
    plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
elseif num_dim_tsne==3
    for i_k=1:size(gmfit.mu,1)
        hold on;
        
        plot_ellipsoid(gmfit.Sigma(:,:,i_k),gmfit.mu(i_k,:),alpha_val,colormap_val);
    end
    plot3(gmfit.mu(:,1),gmfit.mu(:,2),gmfit.mu(:,3),'kx','LineWidth',2,'MarkerSize',10)
end