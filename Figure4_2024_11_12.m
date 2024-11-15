function Figure3_2024_11_12
%Figure 3 Version 1
%Plotting different behaviour maps
%Figure has panel A and B

%Options:
%Function to plot single cell trajcejotries
plot_all_trajectories=1;
plot_mean_activity_and_hetero=0;
save_fig_as_pdf=1;
do=[4,7,9,10,11,15,5];

%Path to data
in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05_grid';
data_name='behavior_map_tau_0.2-v0_0.01-n_2.0-nu_0.05.mat';
t_name='tau=0.2, nu=2.0, v_{0}=0.01, n=2';
out_path='Matlab_Fig';
current_directory=pwd;

%Loading color code map
color_code_in=load("color_codes.mat");
color_code=color_code_in.color_codes;
un_b=fieldnames(color_code);

%Peparing data to plot
%figure('Position',[1,1,1920,700],'color','w');
%figure('centimeters',[1,1,1920,700],'color','w');

fig_para.Linewidth=2;
fig_para.FontSize=8;
fig_para.Fontweight='bold';
fig_para.FontName='Arial';
fig_para.Position=[1, 1, 35, 12];
fig_para.hline_width=1;
fig_para.PlotLinewidth=1;
fig_para.gap_h=0.1;
fig_para.margin_h=[0.15,0.025];


figure('Color','w');
set(gcf, 'Units', 'centimeters', ...
    'Position', fig_para.Position, 'PaperUnits', 'centimeters', 'PaperType','a4',...
    'PaperPosition',[1,1,22,12]);


ax=zeros(1,length(in_path));
%for now=1:length(in_path)
    %Loading Data
    in_data=load([in_path,'/',data_name]);

    %Getting S and D
    bg=in_data.behaviours;
    S=in_data.S_grid;
    D=in_data.D_grid;
    [x,y]=meshgrid(S,D);
    
    %Finding used colors
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


    %mapping behaviour map point (text) to color (number)
    bg_out=zeros(size(bg));
    for i=1:length(un_b)
        for j=1:size(bg,1)
            for k=1:size(bg,2)
                bg_out(j,k)=find(strcmp(un_b_used,bg{j,k}));
            end
        end
    end
    
    %Plotting and making figure nice
    %ax(now)=subplot(1,2,now);
    p=surf(x,y,bg_out,'EdgeColor','none');
    set(gca,'xscale','log','yscale','log','zscale','log');
    set(gca,'TickDir','both');
    set(gca,'Linewidth',2,'ylim',[1,100]);
    set(gca,'Layer','top');
    view(0,90);
    box on;
    
    colormap(cmap);
    set(gca,'FontSize',16,'FontWeight','bold','color',[1,1,1],'FontName','Verdana');
    xlabel('D [au]','FontSize',20);
    ylabel('S [au]','Fontsize',20);
    title(t_name);
%end
set(gcf, 'InvertHardcopy', 'off')
% if save_fig_as_pdf==1
%     exportgraphics(gcf,[current_directory,'/', out_path,'/Figure3.pdf'],'ContentType','vector');
%     %fig
%     exportgraphics(gcf,[current_directory,'/', out_path,'/Figure3.eps'],'ContentType','vector');
% end


%Getting points where plotted
%in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05';
Di=dir([in_path,'/S*']);
S=zeros(1,length(Di));
D=zeros(1,length(Di));
ind=1;
do=[4,7,9,10,11,15,5];
for i=do
    f=strfind(Di(i).name,'-');
    name_now=Di(i).name(1:f(2)-1);
    f2=strfind(Di(i).name,'_');
    S(ind)=str2num(name_now(f2(1)+1:f(1)-1));
    D(ind)=str2num(name_now(f2(2)+1:f(2)-1));
    ind=ind+1;
end
t_name=Di(i).name(f(2)+1:end-4);
title(strrep(strrep(t_name,'_',': '),'-',' '));

%goodones=D==5|D==6|D>=15;
%hold on; plot3(D(goodones),S(goodones),ones(1,sum(goodones))*10,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white','Linewidth',2);

hold on; plot3(D,S,ones(1,length(D))*10,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white','Linewidth',2);
% [D,ind]=sort();
% D=D(ind);
fig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


%getting used conditions
S_pre=cell(1,7);
D_pre=cell(1,7);
S_pre_sort=nan(1,7);
D_pre_sort=nan(1,7);


ind=1;

for i=do
    f=strfind(Di(i).name,'-');
    name_now=Di(i).name(1:f(2)-1);
    f2=strfind(Di(i).name,'_');
    %S_pre_sort(i)=str2num(name_now(f2(1)+1:f(1)-1));
    %D_pre_sort(i)=str2num(name_now(f2(2)+1:f(2)-1));
    S_pre{ind}=name_now(f2(1)+1:f(1)-1);
    D_pre{ind}=name_now(f2(2)+1:f(2)-1);
    ind=ind+1;
end
S_pre_u=unique(S_pre);
D_pre_u=unique(D_pre);
[D_y,s_sort]=sort(cellfun(@(a) str2num(a),S_pre_u));
[D_x,d_sort]=sort(cellfun(@(a) str2num(a),D_pre_u));

S=S_pre_u(s_sort);
D=D_pre_u(d_sort);



tvnn=Di(i).name(f(2):end);



if plot_all_trajectories==1
    figure('Color','w');
    set(gcf, 'Units', 'centimeters', ...
        'Position', fig_para.Position, 'PaperUnits', 'centimeters', 'PaperType','a4',...
        'PaperPosition',[1,1,22,12]);
    
    ind=1;
    sgtitle(strrep(strrep(tvnn(1:end-4),'_',':'),'-',' '));
    data_out=cell(length(S),length(D));
    for S_ind=flip(1:length(S))
        for D_ind=1:length(D)
            t_name=['S_',num2str(S{S_ind}),'-D_',num2str(D{D_ind})];
            data=load([in_path,'/S_',num2str(S{S_ind}),'-D_',num2str(D{D_ind}),tvnn]);
            plot_now=squeeze(data.traces(:,1,:));
            subplot(length(S),length(D),ind);
    
            plot_now_x=plot_now(:,1);
            plot_now_y=plot_now(:,2:end);
            plot(plot_now_x,plot_now_y);
    
            if ind==1
                data_x=plot_now_x;
            end
            data_out{S_ind,D_ind}=plot_now_y;
    
            title(strrep(strrep(t_name,'_',':'),'-',' '));
            
            a=axis;
            axis([-50,a(2),0,1.2])
            if mod(ind,length(D))==1
                ylabel('Sigma Factor [au]');
            end
            if ind>length(D)*(length(S)-2)
                xlabel('Time [au]');
            end
            ind=ind+1;
            box on;
            set(gca,'Linewidth',2,'Fontweight','bold','Fontsize',12);

        end
    end

end



