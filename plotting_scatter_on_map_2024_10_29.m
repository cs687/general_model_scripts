

in_path='/Users/cs687/Documents/General Model Article/Simulation_Data/Behaviour_map_tau_0.2-v0_0.01-n_2.0-nu_0.05';
Di=dir([in_path,'/S*']);
S=zeros(1,length(Di));
D=zeros(1,length(Di));
for i=1:length(Di)
    f=strfind(Di(i).name,'-');
    name_now=Di(i).name(1:f(2)-1);
    f2=strfind(Di(i).name,'_');
    S(i)=str2num(name_now(f2(1)+1:f(1)-1));
    D(i)=str2num(name_now(f2(2)+1:f(2)-1));
end

hold on; plot3(D,S,ones(1,48)*10,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white','Linewidth',2);
% [D,ind]=sort();
% D=D(ind);
fig;
