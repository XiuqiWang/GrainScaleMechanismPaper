function [ParticleID,E,VX,VExVector,VEzVector,Mass_tot,VEx,VEz,ME,MD] = storeParticleIDData(ID_Particle,Z,X,Vxstored,Vzstored,Vp,coe_h,dt,N_inter)
ParticleID=[];VExVector=[];VEzVector=[];ME_tot=0;MD_tot=0;VX=[];VZ=[];E=[];
VEx = cell(1,N_inter);VEz = VEx;ME = zeros(1,N_inter);MD=ME;
for i=1:length(ID_Particle)
    g=9.81;Lx=0.00025*100;Ly=2*0.00025;A=Lx*Ly;
    z=Z(:,i);
    
    %calculate average particle velocities from position differences/dt
    Vx=[0; (X(2:end,i)-X(1:end-1,i))/dt];
    Vz=[0; (Z(2:end,i)-Z(1:end-1,i))/dt];
    
    ke=0.5*Vz.^2;%(sqrt(Vxstored(:,i).^2+Vz.^2)).^2;
    pe=g*z;
    e=ke+pe;
    d_h=coe_h*0.00025;
    thre_e=g*d_h;
    [ID_Ei,ID_Di]=OutputID(e,thre_e);
    ParticleID{i}={ID_Ei,ID_Di};
    
    %correct the Vx when an ejected particle crosses the boundary in
    %positive direction
    Index_neg = find(Vx<-Lx*0.25/dt);
    Vx(Index_neg) = Vx(Index_neg)+Lx/dt;
%     Index_pos = find(Vx>Lx*0.5/dt);
%     Vx(Index_pos) = Vx(Index_pos)-Lx/dt;
    
    VX=[VX;Vx'];
    VZ=[VZ;Vz'];
    E=[E;e'];
    VExi=Vx(ID_Ei+1);%Vxstored(ID_Ei+1,i);%
    VEzi=Vz(ID_Ei+1);%Vzstored(ID_Ei+1,i);
    mE=Vp(ID_Ei,i)*2650;
    mD=Vp(ID_Di,i)*2650;
    
    VExVector = [VExVector; VExi];
    VEzVector = [VEzVector; VEzi];
    ME_tot = ME_tot+sum(mE)/5/A;
    MD_tot = MD_tot+sum(mD)/5/A;
    
    %determine which interval the values should go
    for j = 1:length(ID_Ei)
    % 找到值的范围
    idx = ceil((ID_Ei(j)+1) / (size(X,1)/N_inter)); % 每个cell表示一个(size(X,1)/Ninter个元素的范围
    if idx > N_inter
        idx = N_inter; % 防止溢出
    end
    % 将值添加到对应的cell中
    VEx{idx} = [VEx{idx}, VExi(j)];
    VEz{idx} = [VEz{idx}, VEzi(j)];
    ME(idx) = ME(idx)+mE(j)/(5/N_inter)/A;
    end
    %deposition rate versus time
    for j = 1:length(ID_Di)
    % 找到值的范围
    idx = ceil((ID_Di(j)-1) / (size(X,1)/N_inter)); % 每个cell表示一个(size(X,1)/Ninter个元素的范围
    if idx > N_inter
        idx = N_inter; % 防止溢出
    end
    MD(idx) = MD(idx)+mD(j)/(5/N_inter)/A;
    end
end
Mass_tot = [ME_tot;MD_tot];
end

function [ID_E,ID_D]=OutputID(e,thre_e)
ID_E = [];ID_D = [];t=linspace(5/length(e),5,length(e));
g=9.81;
condition_indices = find(e<=thre_e);%find indices with low energy and near bed surface
%disp(['condition_indices:',num2str(condition_indices)]);
segments = [];
% 遍历找到的索引，找出每一段的第一个和最后一个索引
if ~isempty(condition_indices)
    start_idx = condition_indices(1);
    for i = 2:length(condition_indices)
        %the interval between every E and D should be long enough to be
        %considered saltation
        if t(condition_indices(i))-t(condition_indices(i-1))>sqrt((thre_e/g-12*0.00025)*2/g)%condition_indices(i) ~= condition_indices(i-1) + 1
            end_idx = condition_indices(i-1);
            %set an interval between every ID_D and ID_E to avoid counting rebounds
            %for now use 2*0.01s from the coarsest output steps
            if t(end_idx)-t(start_idx) > 0.03
               segments = [segments; start_idx, end_idx];
            end
            start_idx = condition_indices(i);
        end
    end
    % 添加最后一段
    if t(condition_indices(end))-t(start_idx) > 0.03
       segments = [segments; start_idx, condition_indices(end)];
    end
end
if length(segments)> 1
    ID_E = segments(1:end,2);%count all the erosions
    ID_D = segments(2:end,1);%exclude the first deposition to avoid overcounting (happens when a particle has been eroded at the first time step)
%delete the ID_D if it appears at the first time step and the ID_E if it
%appears at the last time step
ID_D = ID_D(ID_D>1);
ID_E = ID_E(ID_E<size(e,1));
end
end






