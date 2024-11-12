%!function load_and_plot_behaviour_map

% in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_3.0-v0_0.01-n_2.0-nu_0.05';
% data_name='behavior_map_tau_3.0-v0_0.01-n_2.0-nu_0.05.mat';
in_path={'/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05',...
    '/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_3.0-v0_0.01-n_2.0-nu_0.05'};
data_name={'behavior_map_tau_0.2-v0_0.01-n_2.0-nu_0.05.mat',...
    'behavior_map_tau_3.0-v0_0.01-n_2.0-nu_0.05.mat'};
% in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05_grid';
% data_name='behavior_map_tau_0.2-v0_0.01-n_2.0-nu_0.05.mat';
%in_data=load('test.mat');
figure('Position',[1,1,1920/2,800],'color','w');
ax=
for now=flip(1:length(in_path))
    in_data=load([in_path{now},'/',data_name{now}]);
    
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
    
    ax(now)=subplot(1,2,now);
    p=surf(x,y,bg_out,'EdgeColor','none');
    set(gca,'xscale','log','yscale','log','zscale','log');
    set(gca,'TickDir','both');
    set(gca,'Linewidth',3,'ylim',[1,100]);
    set(gca,'Layer','top');
    view(0,90);
    box on;
    
    colormap(ax(now),cmap);
    set(gca,'FontSize',16,'FontWeight','bold');
    xlabel('D [au]','FontSize',20);
    ylabel('S [au]','Fontsize',20);
end

% %Getting points where plotted
% %in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05';
% Di=dir([in_path,'/S*']);
% S=zeros(1,length(Di));
% D=zeros(1,length(Di));
% ind=1;
% do=[4,7,9,10,11,15,5];
% for i=do
%     f=strfind(Di(i).name,'-');
%     name_now=Di(i).name(1:f(2)-1);
%     f2=strfind(Di(i).name,'_');
%     S(ind)=str2num(name_now(f2(1)+1:f(1)-1));
%     D(ind)=str2num(name_now(f2(2)+1:f(2)-1));
%     ind=ind+1;
% end
% t_name=Di(i).name(f(2)+1:end-4);
% title(strrep(strrep(t_name,'_',': '),'-',' '));
% 
% %goodones=D==5|D==6|D>=15;
% %hold on; plot3(D(goodones),S(goodones),ones(1,sum(goodones))*10,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white','Linewidth',2);
% 
% hold on; plot3(D,S,ones(1,length(D))*10,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white','Linewidth',2);
% % [D,ind]=sort();
% % D=D(ind);
% fig;
% 
% %imagesc(bg_out);


