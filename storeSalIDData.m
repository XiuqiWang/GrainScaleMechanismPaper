function [Par,VX,ez,exz,Vxr,ez_mean_t,RIM] = storeSalIDData(ID_Particle,VXstored,X,Z,Vp,coe_h,dt,N_inter)
Par={};VX=[];ez=[];exz=[];
Vxr=[];
ez_t = cell(1,N_inter);RIM = zeros(1,N_inter);
for i=1:length(ID_Particle)
    g=9.81;
    Lx=100*0.00025;
    Ly=2*0.00025;
    A=Lx*Ly;
    z=Z(:,i);
    %calculate average particle velocities from position differences/dt
    Vz=[0; (Z(2:end,i)-Z(1:end-1,i))/dt];
    Vx=[0; (X(2:end,i)-X(1:end-1,i))/dt];
    %correct the Vx when an ejected particle crosses the boundary in
    %positive direction
    Index_neg = find(Vx<-Lx*0.1/dt);
    F_amp = ceil(0.5*(VXstored(Index_neg,i)+VXstored(Index_neg,i))*dt/Lx);
    Vx(Index_neg) = Vx(Index_neg)+F_amp*Lx/dt;
    
    ke=0.5*Vz.^2;
    pe=g*z;
    e=ke+pe;
    d_h=coe_h*0.00025;
    thre_e=g*d_h;
    thres_Zsal=(coe_h-12)*0.00025;
    mp=Vp(1,i)*2650;
    [Moms_col, IDvzri_vec, ezi, exzi, Vxr_i]=findSaltationID(e,Vx,Vz,VXstored(:,i),thre_e,dt);
    %collect the global vectors
    Par{i}={Moms_col, IDvzri_vec};
    VX=[VX;Vx'];
    ez = [ez,ezi];
    exz = [exz,exzi];
    Vxr = [Vxr,Vxr_i];
    
    %time-varying: determine which interval the values should go
    for j = 1:size(IDvzri_vec,1)
    % 按照rebound的时间找到值的范围
    idx = ceil((IDvzri_vec(j,2)+1) / (size(X,1)/N_inter)); % 每个cell表示一个(size(X,1)/Ninter个元素的范围
    if idx > N_inter
        idx = N_inter; % 防止溢出
    end
    % 将值添加到对应的cell中
    ez_t{idx} = [ez_t{idx}, ezi(j)];
    RIM(idx) = RIM(idx)+mp/(5/N_inter)/A;
    end
end
ez_mean_t=NaN(1,N_inter);
for i=1:size(ez_mean_t,2)
    if ~isempty(ez_t{i})
    ez_mean_t(3,i)=mean(ez_t{i});
    end
end
end

%%找出mobile的ID数组，遍历找saltation
function [Moments_col, IDvzri_vec, ez_vector, exz_vector, Vxr_vec]=findSaltationID(e,Vxi,Vzi,Vxstoredi,thre_e,dt)
Moments_col = [];ez_vector = [];exz_vector = [];
IDvzri_vec=[]; 
Vxr_vec=[];
%get the mobile intervals 
IntervalMobile=OutputMobileInterval(e,thre_e);
t=linspace(dt,5,5/dt);
%detect saltations in the mobile intervals 
% 找到峰值的索引
for i=1:size(IntervalMobile,1)
Vz_sal=Vzi(IntervalMobile(i,1):IntervalMobile(i,2));
%time vector
ti=t(IntervalMobile(i,1):IntervalMobile(i,2));

% % 初始化存储谷值的索引数组
% troughIndices = [];
% 找到从负到正穿越的索引
crossing_indices = find(Vz_sal(1:end-1) < 0 & Vz_sal(2:end) >= 0);
% GlobalI=linspace(IntervalMobile(i,1),IntervalMobile(i,2),IntervalMobile(i,2)-IntervalMobile(i,1)+1);
% CI=GlobalI(crossing_indices);
% 计算每次穿越X轴的时刻
ratioST = Vz_sal(crossing_indices)./(Vz_sal(crossing_indices+1) - Vz_sal(crossing_indices));
crossing_times = ti(crossing_indices) - dt*ratioST';

%locate the next ids (global) of the crossing moments
ID_next=ceil(crossing_times/dt);
Vxzi = sqrt(Vxi.^2+Vzi.^2);
IDvzri=Findez(ID_next,Vzi,Vxzi,thre_e);

ez = abs(Vzi(IDvzri(:,2))./Vzi(IDvzri(:,1)));
exz = Vxzi(IDvzri(:,2))./Vxzi(IDvzri(:,1));

Moments_col = [Moments_col, crossing_times];
ez_vector = [ez_vector, ez'];
exz_vector = [exz_vector, exz'];
IDvzri_vec = [IDvzri_vec; IDvzri];
vxri = Vxi(IDvzri(:,2));
Vxr_vec = [Vxr_vec, vxri'];
%end
end
end

function segments=OutputMobileInterval(e,thre_e)
t=linspace(5/length(e),5,length(e));
g=9.81;
condition_indices = find(e>thre_e);%find indices with low energy and near bed surface
%disp(['condition_indices:',num2str(condition_indices)]);
segments = [];
% 遍历找到的索引，找出每一段的第一个和最后一个索引
if ~isempty(condition_indices)
    start_idx = condition_indices(1);
    for i = 2:length(condition_indices)
        %set a limit to every static interval to avoid counting rebounds
        %for now use 3*0.01s from the coarsest output steps
        if t(condition_indices(i))-t(condition_indices(i-1)) > 0.03%condition_indices(i) ~= condition_indices(i-1) + 1
           end_idx = condition_indices(i-1);
           %set a limit to every mobile interval so that the saltation is long enough
           if t(end_idx)-t(start_idx) > sqrt((thre_e/g-12*0.00025)*2/g)
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
end

function [IDvzri]=Findez(ID_next,Vzi,Vxzi,thre_e)
IDvzr=[];IDvzi=[];IDvzri=[];
for i=1:length(ID_next)-1
    ID_range = linspace(ID_next(i)-1,ID_next(i+1),ID_next(i+1)-ID_next(i)+2);
    if length(ID_range)>2
    local_maxima = Vzi(ID_range(2):ID_range(end-1))>Vzi(ID_range(1):ID_range(end-2)) & Vzi(ID_range(2):ID_range(end-1))>Vzi(ID_range(3):ID_range(end));
    local_minimu = Vzi(ID_range(2):ID_range(end-1))<Vzi(ID_range(1):ID_range(end-2)) & Vzi(ID_range(2):ID_range(end-1))<Vzi(ID_range(3):ID_range(end));
%     disp('find local_maxima');
%     disp(find(local_maxima));
    end
    IDvzr=[IDvzr, ID_range(find(local_maxima)+1)];
    IDvzi=[IDvzi, ID_range(find(local_minimu)+1)];
%     Vz_re=Vzi(ID_range(local_maxima));
%     Vz_im=Vzi(ID_range(local_minimu));
end
%exclude the first and the last id in each mobile interval
IDvzri=[IDvzi(1:end-1)', IDvzr(2:end)'];
%only include the different-sign pairs to avoid some local spikes in the Vz profile
product=Vzi(IDvzri(:,1)).*Vzi(IDvzri(:,2));
IDvzri=IDvzri(product<0,:);
%the rebound velocity should be high enough to reach (coe_h-12)D
vz_crit = sqrt((thre_e/9.81-12*0.00025)*2*9.81);
IDvzri=IDvzri(Vzi(IDvzri(:,2))>vz_crit,:);
%exclude the ones with exz > 1
exz=Vxzi(IDvzri(:,2))./Vxzi(IDvzri(:,1));
IDvzri=IDvzri(exz<1,:);
end
