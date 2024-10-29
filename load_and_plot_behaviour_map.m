!function load_and_plot_behaviour_map

in_data=load('test.mat');

%bg=flipud(in_data.behaviours);
bg=in_data.behaviours;
S=in_data.S_grid;
D=in_data.D_grid;
[x,y]=meshgrid(S,D);

color_code_in=load("color_codes.mat");
color_code=color_code_in.color_codes;

un_b=fieldnames(color_code);
un_b_used=unique(bg);
[~,id]=ismember(un_b,un_b_used);
f=find(id);

cmap=zeros(3,length(f));
ind=0;
for i=1:length(f)
    ind=ind+1;
    cmap(:,id(f(i)))=color_code.(un_b{f(i)});
end
cmap=cmap';

%
bg_out=zeros(size(bg));
for i=1:length(un_b)
    for j=1:size(bg,1)
        for k=1:size(bg,2)
            bg_out(j,k)=find(strcmp(un_b_used,bg{j,k}));
        end
    end
end

figure;
p=surf(x,y,bg_out,'EdgeColor','none');
set(gca,'xscale','log','yscale','log','zscale','log');
set(gca,'TickDir','both');
set(gca,'Linewidth',3,'ylim',[1,100]);
view(0,90);
box on;

colormap(cmap);
xlabel('D','FontSize',20);
ylabel('S','Fontsize',20);

%imagesc(bg_out);


