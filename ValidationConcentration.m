%validate the slope of concentration profile/with two numbers of output step 
close all
% clear
% data4=load('matData\storedfromData\S005M20LBIni.mat');
% data5=load('matData\storedfromData\FTSS005M20LBIni.mat');
% data6=load('matData\storedfromData\FiftyTSS005M20LBIni.mat');
% 
g=9.81;
for i=4:6
    num_p=2725;
eval(['ID_Particle',num2str(i),'=linspace(num_p-10-299,num_p-10,300);']);
eval(['[X',num2str(i),', Y',num2str(i),', Z',num2str(i),', Vx',num2str(i),', Vy',num2str(i),', Vz',num2str(i),', Vp',num2str(i),',W',num2str(i),']=getDryParticleInfoFromData(data',num2str(i),',',num2str(i),',ID_Particle',num2str(i),');']);
%eval(['KE',num2str(i),'=0.5*(sqrt(Vx',num2str(i),'.^2 + Vy',num2str(i),'.^2 + Vz',num2str(i),'.^2)).^2;']);
end

D=0.00025;
%t=linspace(0, 5, size(Z0,1));
coe_h = 15;%critial height for a mobile particle to reach
N_inter = 250;%number of output timesteps for erosion and deposition properties
thres_Zbed = 13*D;
dt4=5/501;
dt5=5/5004;
dt6=5/50803;

%access all the ids of all the particle motions
for i=4:6
eval(['[Par',num2str(i),',E',num2str(i),',VXCal',num2str(i),',VExVector',num2str(i),',VEzVector',num2str(i),',Mass',num2str(i),',VEx',num2str(i),',VEz',num2str(i),',ME',num2str(i),',MD',num2str(i),'] = storeParticleIDData(ID_Particle',num2str(i),',Z',num2str(i),',X',num2str(i),',Vx',num2str(i),',Vz',num2str(i),',Vp',num2str(i),',coe_h,dt',num2str(i),',N_inter);']);
end

%verification
t=linspace(0, 5, size(Z4,1));t5=linspace(0, 5, size(Z5,1));t6=linspace(0, 5, size(Z6,1));
colors = [0 0 0;1 0 0;0 0 1;0.4660 0.6740 0.1880];
id_p=294;
figure
subplot(3,1,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t,Z4(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t(Par4{id_p}{1}), Z4(Par4{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t(Par4{id_p}{2}), Z4(Par4{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
%plot(t(ParS4{id_p}{1}), Z4(ParS4{id_p}{1},id_p),'diamond','MarkerSize',8,'DisplayName','Saltation rebound');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$Z_\mathrm{p}$ [m]','Interpreter','Latex','fontsize',16);
legend('fontsize',12);
subplot(3,1,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t5,Z5(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t5(Par5{id_p}{1}), Z5(Par5{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t5(Par5{id_p}{2}), Z5(Par5{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$Z_\mathrm{p}$ [m]','Interpreter','Latex','fontsize',16);
subplot(3,1,3)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t6,Z6(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t6(Par6{id_p}{1}), Z6(Par6{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t6(Par6{id_p}{2}), Z6(Par6{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$Z_\mathrm{p}$ [m]','Interpreter','Latex','fontsize',16);



figure
subplot(3,1,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t,E4(id_p,:),'-','marker','.','DisplayName','Particle trajectory');
hold on;
% plot(t(Par4{id_p}{1}), E4(id_p,Par4{id_p}{1}),'o','MarkerSize',10,'DisplayName','Erosion');
% plot(t(Par4{id_p}{2}), E4(id_p,Par4{id_p}{2}),'x','MarkerSize',10,'DisplayName','Deposition');
plot(t,ones(size(t))*coe_h*g*D,'--','DisplayName','Energy criteria');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$E_\mathrm{p}$ []','Interpreter','Latex','fontsize',16);
legend('fontsize',12);
subplot(3,1,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t5,E5(id_p,:),'-','marker','.','DisplayName','Particle trajectory');
hold on;
% plot(t5(Par5{id_p}{1}), E5(id_p,Par5{id_p}{1}),'o','MarkerSize',10,'DisplayName','Erosion');
% plot(t5(Par5{id_p}{2}), E5(id_p,Par5{id_p}{2}),'x','MarkerSize',10,'DisplayName','Deposition');
plot(t5,ones(size(t5))*coe_h*g*D,'--','DisplayName','Energy criteria');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$E_\mathrm{p}$ []','Interpreter','Latex','fontsize',16);
subplot(3,1,3)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t6,E6(id_p,:),'-','marker','.','DisplayName','Particle trajectory');
hold on;
% plot(t6(Par6{id_p}{1}), E6(id_p,Par6{id_p}{1}),'o','MarkerSize',10,'DisplayName','Erosion');
% plot(t6(Par6{id_p}{2}), E6(id_p,Par6{id_p}{2}),'x','MarkerSize',10,'DisplayName','Deposition');
plot(t6,ones(size(t6))*coe_h*g*D,'--','DisplayName','Energy criteria');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$E_\mathrm{p}$ []','Interpreter','Latex','fontsize',16);


figure
subplot(3,1,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t,Vz4(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t(Par4{id_p}{1}), Vz4(Par4{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t(Par4{id_p}{2}), Vz4(Par4{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$U_\mathrm{p,z}$ [m/s]','Interpreter','Latex','fontsize',16);
legend('fontsize',12);
subplot(3,1,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t5,Vz5(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t5(Par5{id_p}{1}), Vz5(Par5{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t5(Par5{id_p}{2}), Vz5(Par5{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$U_\mathrm{p,z}$ [m/s]','Interpreter','Latex','fontsize',16);
subplot(3,1,3)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t6,Vz6(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t6(Par6{id_p}{1}), Vz6(Par6{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t6(Par6{id_p}{2}), Vz6(Par6{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$U_\mathrm{p,z}$ [m/s]','Interpreter','Latex','fontsize',16);


% data54=readMercuryCG('..\WilletViscousTransport\stat\LBIniTP\S005M20LBIni.stat');
% data_ini=readMercuryCG('..\WilletViscousTransport\stat\LBIniTP\S005M20LBIniFirstTimeStep.stat');
% dz=0.00025;
% for i=4:4
%     eval(['[Q5',num2str(i),', Upx5', num2str(i), ', Cvx5', num2str(i),', H5', num2str(i),'] = CharTester(data5', num2str(i),',5',',D, dz);']); 
% end
% [Qini,Upini,Cvxini,Hini]=CharTester(data_ini,5,D,dz);
% Cini=Cvxini(1);

%calculate the slope
t_cg=linspace(0.01,5,501);
n=4;slopes=zeros(5,(length(Cvx54)-1)/n);
for i=4:4
    eval(['slopes(',num2str(i+1),',:) = CalSlope(n,t_cg,Cvx5',num2str(i),');']);
end

n_dis=size(X4,1)/N_inter;
t_inter = linspace((t_cg(1)+t_cg(1+n_dis))*0.5,(t_cg(end)+t_cg(end-n_dis))*0.5,(length(t_cg)-1)/n_dis);%Discrete t
t_plot=linspace((t_cg(1)+t_cg(1+n))*0.5,(t_cg(end)+t_cg(end-n))*0.5,(length(t_cg)-1)/n);%CG t
% figure
% for i=0:4
%     subplot(2,3,i+1)
%     plot([0, t_plot],[0, slopes(i+1,:)],'k--');
%     hold on
%     eval(['plot([0, t_inter], [0, Me',num2str(i),'-Md',num2str(i),'],"r-");']);
%     xlabel('time [s]');ylabel('$\frac{\partial C_\mathrm{sal}}{\partial t}$ [kg/m$^2$/s]','Interpreter','Latex');
% end
% legend('CG concentration growth rate','Discrete net erosion rate');

%compare 500 output steps and 5000
t_cg5 = linspace(0.001,5,5003);
n_dis5 = round(size(X5,1)/N_inter);
t_cg6 = linspace(0.0001,5,50802);
n_dis6 = round(size(X6,1)/N_inter);
t_inter5 = linspace((t_cg5(1)+t_cg5(1+n_dis5))*0.5,(t_cg5(end)+t_cg5(end-n_dis5))*0.5,N_inter);
t_inter6 = linspace((t_cg6(1)+t_cg6(1+n_dis6))*0.5,(t_cg6(end)+t_cg6(end-n_dis6))*0.5,N_inter);
figure
subplot(3,1,1)
plot([0, t_inter], [0, ME4],'k-');
hold on;
plot([0, t_inter5], [0, ME5],'r-');
plot([0, t_inter6], [0, ME6],'b-');
xlabel('time [s]');ylabel('Erosion [kg/m^2/s]');
subplot(3,1,2)
plot([0, t_inter], [0, MD4],'k-');
hold on;
plot([0, t_inter5], [0, MD5],'r-');
plot([0, t_inter6], [0, MD6],'b-');
xlabel('time [s]');ylabel('Deposition [kg/m^2/s]');
subplot(3,1,3)
plot([0, t_inter], [0, ME4-MD4],'k-');
hold on;
plot([0, t_inter5], [0, ME5-MD5],'r-');
plot([0, t_inter6], [0, ME6-MD6],'b-');
plot([0, t_plot],[0, slopes(5,:)],'g--');
xlabel('time [s]');ylabel('Net erosion [kg/m^2/s]');
legend('500 outputs','5000 outputs','50000 outputs','CG concentration growth rate');

colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
figure
subplot(1,2,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
        eval(['[x,y]=getXYfromVar(VExVector',num2str(i),',0.05);']);
        stairs(x, y, 'LineWidth', 1);
        hold on;
end
%xlim([-1.5 2]);
box on;grid on;
xlabel('$U_\mathrm{E,x}$ [m/s]','Interpreter','Latex');
ylabel('Frequency');
subplot(1,2,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
        eval(['[x,y]=getXYfromVar(VEzVector',num2str(i),',0.05);']);
        stairs(x, y, 'LineWidth', 1);
        hold on;
end
%xlim([-1.5 2]);
box on;grid on;
xlabel('$U_\mathrm{E,z}$ [m/s]','Interpreter','Latex');
legend('500 outputs','5000 outputs','50000 outputs');


function slopes=CalSlope(n,t,C)
dx = t(n+1:n:end)-t(1:n:end-n);
dy = C(n+1:n:end)-C(1:n:end-n);
% 计算每个点的斜率
slopes = dy ./ dx;
end