function plotting_histogram_delay

figure('Position',[0,0,1200,1900],'Color','w');
b=bar([30,3;13,0;20,0]);
b(1).FaceColor='y';
b(2).FaceColor='c';

% bar([1,3,5],[30,13,20],'FaceColor','g');
% hold on;
% bar([2,4,6],[3,0,0],'FaceColor','c');
box on;
xlabel('Strain');
ylabel('Fold Change');
set(gca,'Linewidth',3,'FontSize',20,'XTicklabel',{'WT','ECF41','UN'});
a=axis;
axis([a(1),a(2),a(3),35]);
legend({'YFP','CFP'});
title('Induction after 2h');