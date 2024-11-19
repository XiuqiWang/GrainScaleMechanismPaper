close all
clear
% data40=load('matData\storedfromData\S004DryLBIni.mat');
% data41=load('matData\storedfromData\S004M1LBR1Ini.mat');
% data42=load('matData\storedfromData\S004M5LBIni.mat');
% data43=load('matData\storedfromData\S004M10LBIni.mat');
% data44=load('matData\storedfromData\S004M20LBIni.mat');
data50=load('matData\storedfromData\S005DryLBIni.mat');
data51=load('matData\storedfromData\S005M1LBIni.mat');
data52=load('matData\storedfromData\S005M5LBR1Ini.mat');
data53=load('matData\storedfromData\S005M10LBR1Ini.mat');
data54=load('matData\storedfromData\S005M20LBIni.mat');
% data60=load('matData\storedfromData\S006DryLBIni.mat');
% data61=load('matData\storedfromData\S006M1LBIni.mat');
% data62=load('matData\storedfromData\S006M5LBIni.mat');
% data63=load('matData\storedfromData\S006M10LBIni.mat');
% data64=load('matData\storedfromData\S006M20LBIni.mat');

for j=5:5
for i=0:4
    if (j==5 && (i == 2 || i == 3)) || (j==4 && i==1)
        num_p=2714;
    else
        num_p=2725;
    end
eval(['ID_Particle',num2str(j),num2str(i),'=linspace(num_p-10-299,num_p-10,300);']);
eval(['[X',num2str(j),num2str(i),', Y',num2str(j),num2str(i),', Z',num2str(j),num2str(i),', Vx',num2str(j),num2str(i),', Vy',num2str(j),num2str(i),', Vz',num2str(j),num2str(i),', Vp',num2str(j),num2str(i),',W',num2str(j),num2str(i),']=getDryParticleInfoFromData(data',num2str(j),num2str(i),',',num2str(i),',ID_Particle',num2str(j),num2str(i),');']);
end
end

D=0.00025;
dt=0.01;
coe_h = 20;%critial height for a mobile particle to reach
N_inter = 100;%number of output timesteps for erosion and deposition properties
omega = [0 1 5 10 20];

%access all the ids of all the particle motions
for j=5:5
for i=0:4
eval(['[Par',num2str(j),num2str(i),',E',num2str(j),num2str(i),',VXCal',num2str(j),num2str(i),',VExVector',num2str(j),num2str(i),',VEzVector',num2str(j),num2str(i),',Mass',num2str(j),num2str(i),',VEx',num2str(j),num2str(i),',VEz',num2str(j),num2str(i),',ME',num2str(j),num2str(i),',MD',num2str(j),num2str(i),'] = storeParticleIDData(ID_Particle',num2str(j),num2str(i),',Z',num2str(j),num2str(i),',X',num2str(j),num2str(i),',Vx',num2str(j),num2str(i),',Vz',num2str(j),num2str(i),',Vp',num2str(j),num2str(i),',coe_h,dt,N_inter);']);
end
end

%verification
colors = [0 0 0;1 0 0;0 0 1;0.4660 0.6740 0.1880];
t=linspace(0, 5, 500);
id_p=294;
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(t,Z54(:,id_p),'-','marker','.','DisplayName','Particle trajectory');
hold on;
plot(t(Par54{id_p}{1}), Z54(Par54{id_p}{1},id_p),'o','MarkerSize',10,'DisplayName','Erosion');
plot(t(Par54{id_p}{2}), Z54(Par54{id_p}{2},id_p),'x','MarkerSize',10,'DisplayName','Deposition');
%plot(t(ParS4{id_p}{1}), Z4(ParS4{id_p}{1},id_p),'diamond','MarkerSize',8,'DisplayName','Saltation rebound');
xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$Z_\mathrm{p}$ [m]','Interpreter','Latex','fontsize',16);
legend('fontsize',12);
% 
% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% plot(t,VXCal4(id_p,:),'-','marker','.','DisplayName','Calculated');
% hold on;
% plot(t,Vx4(:,id_p),'-','marker','.','DisplayName','Stored');
% plot(t(Par4{id_p}{1}), VXCal4(id_p,Par4{id_p}{1}),'o','MarkerSize',10,'DisplayName','Erosion');
% plot(t(Par4{id_p}{2}), VXCal4(id_p,Par4{id_p}{2}),'x','MarkerSize',10,'DisplayName','Deposition');
% xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$U_\mathrm{p,x}$ [m/s]','Interpreter','Latex','fontsize',16);
% legend('fontsize',12);
% 
% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% plot(t,E4(id_p,:),'-','marker','.','DisplayName','Energy');
% hold on;
% plot(t(Par4{id_p}{1}), E4(id_p,Par4{id_p}{1}),'o','MarkerSize',10,'DisplayName','Erosion');
% plot(t(Par4{id_p}{2}), E4(id_p,Par4{id_p}{2}),'x','MarkerSize',10,'DisplayName','Deposition');
% xlabel('$time$ [s]','Interpreter','Latex','fontsize',16);ylabel('$E_\mathrm{p}$ [-]','Interpreter','Latex','fontsize',16);
% legend('fontsize',12);


%distribution
% colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
% figure
% for j=4:6
% subplot(2,3,j-3)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
%         eval(['[x,y]=getXYfromVar(VExVector',num2str(j),num2str(i),',0.05);']);
%         stairs(x, y, 'LineWidth', 1);
%         hold on;
% end
% box on;grid on;
% xlim([-1 2.5]);ylim([0 0.12]);
% xlabel('$U_\mathrm{E,x}$ [m/s]','Interpreter','Latex');
% ylabel('Normalized frequency');
% eval(['title("$\Theta=0.0',num2str(j),'$","Interpreter","Latex");']);
% end
% for j=4:6
% subplot(2,3,j)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
%         eval(['[x,y]=getXYfromVar(VEzVector',num2str(j),num2str(i),',0.05);']);
%         stairs(x, y, 'LineWidth', 1);
%         hold on;
% end
% box on;grid on;
% xlim([0 2]);ylim([0 0.4]);
% xlabel('$U_\mathrm{E,z}$ [m/s]','Interpreter','Latex');
% ylabel('Normalized frequency');
% end 
% legend('$\Omega$=0','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex','fontsize',10);


%mean and std
VxEskew=[];VzEskew=[];VxEkurt=[];VzEkurt=[];
VxEmean=[];VzEmean=[];VxEstd=[];VzEstd=[];
VxEQ1=[];VxEQ2=[];VxEQ3=[];VzEQ1=[];VzEQ2=[];VzEQ3=[];
for j=5:5
for i=0:4
eval(['VxEmean(',num2str(j-3),',',num2str(i),'+1)=mean(VExVector',num2str(j),num2str(i),');']);
eval(['VxEstd(',num2str(j-3),',',num2str(i),'+1)=std(VExVector',num2str(j),num2str(i),');']);
eval(['VzEmean(',num2str(j-3),',',num2str(i),'+1)=mean(VEzVector',num2str(j),num2str(i),');']);
eval(['VzEstd(',num2str(j-3),',',num2str(i),'+1)=std(VEzVector',num2str(j),num2str(i),');']);
eval(['VxEQ1(',num2str(j-3),',',num2str(i),'+1)=prctile(VExVector',num2str(j),num2str(i),',25);']);
eval(['VxEQ2(',num2str(j-3),',',num2str(i),'+1)=prctile(VExVector',num2str(j),num2str(i),',50);']);
eval(['VxEQ3(',num2str(j-3),',',num2str(i),'+1)=prctile(VExVector',num2str(j),num2str(i),',75);']);
eval(['VzEQ1(',num2str(j-3),',',num2str(i),'+1)=prctile(VEzVector',num2str(j),num2str(i),',25);']);
eval(['VzEQ2(',num2str(j-3),',',num2str(i),'+1)=prctile(VEzVector',num2str(j),num2str(i),',50);']);
eval(['VzEQ3(',num2str(j-3),',',num2str(i),'+1)=prctile(VEzVector',num2str(j),num2str(i),',75);']);
eval(['VxEskew(',num2str(j-3),',',num2str(i),'+1)=skewness(VExVector',num2str(j),num2str(i),');']);
eval(['VzEskew(',num2str(j-3),',',num2str(i),'+1)=skewness(VEzVector',num2str(j),num2str(i),');']);
eval(['VxEkurt(',num2str(j-3),',',num2str(i),'+1)=kurtosis(VExVector',num2str(j),num2str(i),');']);
eval(['VzEkurt(',num2str(j-3),',',num2str(i),'+1)=kurtosis(VEzVector',num2str(j),num2str(i),');']);

% eval(['Me(',num2str(i),'+1)=Mass',num2str(i),'(1);']);
% eval(['Md(',num2str(i),'+1)=Mass',num2str(i),'(2);']);
% eval(['Med(',num2str(i),'+1)=Mass',num2str(i),'(1)-Mass',num2str(i),'(2);']);
end
end

% figure;
% subplot(1,2,1);
% for i=1:3
% errorbar(omega,VxEmean(i,:),VxEstd(i,:),'o--','LineWidth', 1.5, 'MarkerSize', 8);
% hold on;
% end
% box on;
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\bar{U}_\mathrm{E,x} \quad \& \quad \sigma_{U_\mathrm{E,x}}$ [m/s]','Interpreter','Latex','fontsize',12);
% subplot(1,2,2);
% for i=1:3
% errorbar(omega,VzEmean(i,:),VzEstd(i,:),'o--','LineWidth', 1.5, 'MarkerSize', 8);
% hold on;
% end
% box on;
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\bar{U}_\mathrm{E,z} \quad \& \quad \sigma_{U_\mathrm{E,z}}$ [m/s]','Interpreter','Latex','fontsize',12);
% sgtitle('Mean and std');
% legend('$\Theta=0.04$','$\Theta=0.05$','$\Theta=0.06$','Interpreter','Latex','fontsize',10);

%%central mark: median; edges: 25% and 75%; whiskers|: extreme data not
%%considered outliers; '+' marks: outliers
for i=4:6
eval(['VExMatrix',num2str(i),'=[VExVector',num2str(i),'0;VExVector',num2str(i),'1;VExVector',num2str(i),'2;VExVector',num2str(i),'3;VExVector',num2str(i),'4];']);
eval(['VEzMatrix',num2str(i),'=[VEzVector',num2str(i),'0;VEzVector',num2str(i),'1;VEzVector',num2str(i),'2;VEzVector',num2str(i),'3;VEzVector',num2str(i),'4];']);
eval(['omegamatrix',num2str(i),'=[omega(1)*ones(length(VExVector',num2str(i),'0),1);omega(2)*ones(length(VExVector',num2str(i),'1),1);omega(3)*ones(length(VExVector',num2str(i),'2),1);omega(4)*ones(length(VExVector',num2str(i),'3),1);omega(5)*ones(length(VExVector',num2str(i),'4),1)];']);
end
omegalabels={'0','1','5','10','20'};

figure
for i=1:3
subplot(1,3,i);
eval(['boxplot(VExMatrix',num2str(i+3),',omegamatrix',num2str(i+3),',"Labels",omegalabels);']);
%errorbar(omega, VxEQ2(i,:),VxEQ2(i,:)-VxEQ1(i,:),VxEQ3(i,:)-VxEQ2(i,:),'o--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
box on;ylim([-1 2.5]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_\mathrm{E,x}$ [m/s]','Interpreter','Latex','fontsize',12);
eval(['title("$\Theta=0.0',num2str(i+3),'$","Interpreter","Latex");']);
end
figure
for i=1:3
subplot(1,3,i);
    eval(['boxplot(VEzMatrix',num2str(i+3),',omegamatrix',num2str(i+3),',"Labels",omegalabels);']);
%errorbar(omega, VzEQ2(i,:),VzEQ2(i,:)-VzEQ1(i,:),VzEQ3(i,:)-VzEQ2(i,:),'o--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
box on;ylim([0 2]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_\mathrm{E,z}$ [m/s]','Interpreter','Latex','fontsize',12);
eval(['title("$\Theta=0.0',num2str(i+3),'$","Interpreter","Latex");']);
end
%legend('$\Theta=0.04$','$\Theta=0.05$','$\Theta=0.06$','Interpreter','Latex','fontsize',10);


%erosion rate and deposition rate
colors = [1 0 0;0 0 1;0 0 0;0.9290 0.6940 0.1250];
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
plot(omega,Me,'+--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(omega,Md,'o--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(omega,Med,'*--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
box on;
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$R_\mathrm{mass}$ [kg/m$^2$/s]','Interpreter','Latex','fontsize',12);
%legend('$R_\mathrm{mass,E}$','$R_\mathrm{mass,D}$','$R_\mathrm{mass,E}-R_\mathrm{mass,D}$','Interpreter','Latex'); 

%%values versus time!
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
t_cg=linspace(0.01,5,501);
n_dis=500/N_inter;
t_inter = linspace((t_cg(1)+t_cg(1+n_dis))*0.5,(t_cg(end)+t_cg(end-n_dis))*0.5,(length(t_cg)-1)/n_dis);
% figure
% subplot(2,1,1)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot(t_inter, Emean',num2str(i),'(1,:));']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{U}_\mathrm{E,x}$ [m/s]','Interpreter','Latex','fontsize',12);
% subplot(2,1,2)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot(t_inter, Emean',num2str(i),'(2,:));']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{U}_\mathrm{E,z}$ [m/s]','Interpreter','Latex','fontsize',12);
% legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex');

% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot(t_inter, Emean',num2str(i),'(5,:));']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{\Omega}_\mathrm{p,E}$ [$\%$]','Interpreter','Latex','fontsize',12);
% legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);

% figure
% subplot(3,1,1)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot([0, t_inter], [0, ME',num2str(i),']);']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{R}_\mathrm{E}$ [kg/m$^2$/s]','Interpreter','Latex','fontsize',12);
% subplot(3,1,2)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot([0, t_inter], [0, MD',num2str(i),']);']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{R}_\mathrm{D}$ [kg/m$^2$/s]','Interpreter','Latex','fontsize',12);
% subplot(3,1,3)
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot([0, t_inter], [0, ME',num2str(i),'-MD',num2str(i),']);']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{R}_\mathrm{E}-\bar{R}_\mathrm{D}$ [kg/m$^2$/s]','Interpreter','Latex','fontsize',12);
% legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);

function SEM = sem(data)
    SEM = std(data) / sqrt(length(data));
end