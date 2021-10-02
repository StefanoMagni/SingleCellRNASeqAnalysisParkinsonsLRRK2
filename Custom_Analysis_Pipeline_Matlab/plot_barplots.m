

function plot_barplots(data_vals,colors,leg_cluster,title_eachfig,flag_save,name_save)

h1=figure('Position', [100, 100, 1400, 1400],'PaperOrientation','landscape');
num_x=ceil((size(data_vals,2)+1)/3);
for i_cl=1:size(data_vals,2)
    subplot(3,num_x,i_cl);
    %aHand = axes('parent', h1);
    hold on
    for i=1:size(data_vals,1)
        bar(i,data_vals(i,i_cl), 'facecolor', colors(i+1,:));
    end
    %set(gca, 'XTick', 1:size(prob_cluster_expressed,2), 'XTickLabel', leg_cluster(2:end))
    %TH=rotateticklabel(gca,30);
    title([title_eachfig num2str(i_cl)]);
    if sum(data_vals(:,i_cl))==0
       ylim([0,1]) 
    end
end
subplot(3,num_x,size(data_vals,2)+1);
for i=1:size(data_vals,1)
    bar(i,0, 'facecolor', colors(i+1,:));
    hold on
end
grid off
legend(leg_cluster(2:end))
if flag_save
    print(name_save, '-dpdf','-fillpage')
end