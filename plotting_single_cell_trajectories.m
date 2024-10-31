%Function to plot single cell trajcejotries
plot_all_trajectories=1;
plot_mean_activity_and_hetero=1;


%Getting file loction
p = mfilename('fullpath');
f=strfind(p,'/');
data_path_main=[p(1:f(end-1)),'Simulation_Data/','Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05_grid/'];
%data_path_main=[p(1:f(end-1)),'Simulation_Data/','Behaviour_map_tau_3.0-v0_0.01-n_2.0-nu_0.05/'];
%data_path_main=[p(1:f(end-1)),'Simulation_Data/','Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05/'];

%getting used conditions
Di=dir([data_path_main,'S*']);
S_pre=cell(1,length(Di));
D_pre=cell(1,length(Di));
S_pre_sort=nan(1,length(Di));
D_pre_sort=nan(1,length(Di));

for i=1:length(Di)
    f=strfind(Di(i).name,'-');
    name_now=Di(i).name(1:f(2)-1);
    f2=strfind(Di(i).name,'_');
    %S_pre_sort(i)=str2num(name_now(f2(1)+1:f(1)-1));
    %D_pre_sort(i)=str2num(name_now(f2(2)+1:f(2)-1));
    S_pre{i}=name_now(f2(1)+1:f(1)-1);
    D_pre{i}=name_now(f2(2)+1:f(2)-1);
end
S_pre_u=unique(S_pre);
D_pre_u=unique(D_pre);
[D_y,s_sort]=sort(cellfun(@(a) str2num(a),S_pre_u));
[D_x,d_sort]=sort(cellfun(@(a) str2num(a),D_pre_u));

S=S_pre_u(s_sort);
D=D_pre_u(d_sort);

%Filtering out conditions
%goodones=D_x>15;
% goodones=D_x==6|D_x==5;
% D=D(goodones);

tvnn=Di(i).name(f(2):end);



if plot_all_trajectories==1
    figure('Position',[0,0,1900,1200]);
    
    ind=1;
    sgtitle(strrep(strrep(tvnn(1:end-4),'_',':'),'-',' '));
    data_out=cell(length(S),length(D));
    for S_ind=flip(1:length(S))
        for D_ind=1:length(D)
            t_name=['S_',num2str(S{S_ind}),'-D_',num2str(D{D_ind})];
            data=load([data_path_main,'S_',num2str(S{S_ind}),'-D_',num2str(D{D_ind}),tvnn]);
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
        end
    end
    set(gcf,'color','w');
end

if plot_mean_activity_and_hetero==1

    %Calculating what to plot
    for i=1:size(data_out,1);
        for j=1:size(data_out,2);
            %calulating mean
            m(i,j)=mean(data_out{i,j}(end,:));
            %calculating activation
            data_pre=data_out{i,j}(1:250,:);
            m_pre=mean(data_pre(:));
            s_pre=std(data_pre(:));
            activation(i,j)=sum(data_out{i,j}(end,:)>m_pre+3*s_pre)/size(data_pre,2);
        end
    end

    %plotting mean
    figure('Position',[0,0,1200,1900],'Color','w');
    subplot(2,1,1);
    plot(D_x,m','Linewidth',3);
    xlabel('D [au]');
    ylabel('Mean steady state');
    title('Mean steady state');
    set(gca,'Linewidth',3,'Fontsize',20);
    legend(cellfun(@(a) ['S: ',a],S,'UniformOutput',false));

    %plotting activation
    subplot(2,1,2);
    plot(D_x,activation','Linewidth',3);
    xlabel('D [au]');
    ylabel('Fraction activated [au]');
    title('Activated fraction');
    set(gca,'Linewidth',3,'Fontsize',20);
    legend(cellfun(@(a) ['S: ',a],S,'UniformOutput',false));

    figure('Position',[0,0,1200,1900],'Color','w');
    for i=1:3
        subplot(2,2,i);
        plot(m(i,:),activation(i,:),'x','Linewidth',3,'MarkerSize',10);
        xlabel('Mean steady state');
        ylabel('Activated fraction');
        set(gca,'Linewidth',3,'FontSize',20,'Fontweight','bold');
        title(['S: ',S{i}]);
        axis([0,1.1,0,1.1])
        %names_scatter=cellfun(@(a) ['\leftarrowD:' a],D,'UniformOutput',false);
        %text(m(i,:),activation(i,:),names_scatter,'Fontsize',16);
    end

end

