close all
% clear
% data0=load('matData\storedfromData\S005DryLBIni.mat');
% data1=load('matData\storedfromData\S005M1LBIni.mat');
% data2=load('matData\storedfromData\S005M5LBR1Ini.mat');
% data3=load('matData\storedfromData\S005M10LBR1Ini.mat');
% data4=load('matData\storedfromData\S005M20LBIni.mat');
% data5=load('matData\storedfromData\FTSS005M20LBIni.mat');
% data6=load('matData\storedfromData\FiftyTSS005M20LBIni.mat');

% g=9.81;
% for i=4:6
%     if i == 2 || i == 3
%         num_p=2714;
%     else
%         num_p=2725;
%     end
%     %num_p=2725;
% eval(['ID_Particle',num2str(i),'=linspace(num_p-309,num_p,310);']);
% eval(['[X',num2str(i),', Y',num2str(i),', Z',num2str(i),', Vx',num2str(i),', Vy',num2str(i),', Vz',num2str(i),', Vp',num2str(i),',W',num2str(i),']=getDryParticleInfoFromData(data',num2str(i),',',num2str(i),',ID_Particle',num2str(i),');']);
% end

D=0.00025;
dt4=5/size(Z4,1);dt0=dt4;dt1=dt4;dt2=dt4;dt3=dt4;
dt5=5/size(Z5,1);
dt6=5/size(Z6,1);
t4=linspace(dt4, 5, size(Z4,1));
t5=linspace(dt5, 5, size(Z5,1));
t6=linspace(dt6, 5, size(Z6,1));
coe_h = 20;%critial height for a mobile particle to reach
N_inter = 100;%number of output timesteps for erosion and deposition properties
omega = [0 1 5 10 20];

%access the sal ids of all the particle jumps
for i=4:6
eval(['[ParS',num2str(i),',VXCal',num2str(i),',ez',num2str(i),',exz',num2str(i),',Vxr',num2str(i),',ez_t',num2str(i),',RIM',num2str(i),'] = storeSalIDData(ID_Particle',num2str(i),',Vx',num2str(i),',X',num2str(i),',Z',num2str(i),',Vp',num2str(i),',coe_h,dt',num2str(i),',N_inter);']);
end

%verification
id_p=310;
% figure
% plot(t,Z3(:,id_p),'-','marker','.');
% hold on;
% plot(ParS3{id_p}{1}, Z3(ParS3{id_p}{1},id_p),'x','DisplayName','Saltation');
% xlabel('t [s]');ylabel('Z [m]');
% legend();

figure
subplot(3,1,1)
plot(t4,VZCal4(id_p,:),'-','marker','.');
hold on;
plot(ParS4{id_p}{1}, zeros(length(ParS4{id_p}{1}),1),'x');
plot(t4(ParS4{id_p}{2}(:,1)), VZCal4(id_p,ParS4{id_p}{2}(:,1)),'diamond');
plot(t4(ParS4{id_p}{2}(:,2)), VZCal4(id_p,ParS4{id_p}{2}(:,2)),'o');
xlabel('t [s]');ylabel('Uz [m/s]');
subplot(3,1,2)
plot(t5,VZCal5(id_p,:),'-','marker','.');
hold on;
plot(ParS5{id_p}{1}, zeros(length(ParS5{id_p}{1}),1),'x');
plot(t5(ParS5{id_p}{2}(:,1)), VZCal5(id_p,ParS5{id_p}{2}(:,1)),'diamond');
plot(t5(ParS5{id_p}{2}(:,2)), VZCal5(id_p,ParS5{id_p}{2}(:,2)),'o');
xlabel('t [s]');ylabel('Uz [m/s]');
subplot(3,1,3)
plot(t6,VZCal6(id_p,:),'-','marker','.');
hold on;
plot(ParS6{id_p}{1}, zeros(length(ParS6{id_p}{1}),1),'x');
plot(t6(ParS6{id_p}{2}(:,1)), VZCal6(id_p,ParS6{id_p}{2}(:,1)),'diamond');
plot(t6(ParS6{id_p}{2}(:,2)), VZCal6(id_p,ParS6{id_p}{2}(:,2)),'o');
xlabel('t [s]');ylabel('Uz [m/s]');
legend('Particle velocity','Collision moments','impact','rebound');


% figure
% plot(t(1:200),X3(1:200,id_p),'-','marker','.');
% hold on;
% plot(t(Par3{id_p}{1}), X3(Par3{id_p}{1},id_p),'x','DisplayName','Saltation');
% xlabel('t [s]');ylabel('X [m]');
% legend();
% 
% figure
% plot(t,Vx3(:,id_p),'-','marker','.');
% hold on;
% plot(t(ParS3{id_p}{1}), Vx3(ParS3{id_p}{1},id_p),'x','DisplayName','Saltation');
% xlabel('t [s]');ylabel('Ux [m/s]');
% legend();

%distribution of ezs
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
        eval(['[x,y]=getXYfromVar(ez',num2str(i),',0.05);']);
        stairs(x, y, 'LineWidth', 1);
        hold on;
end
box on;grid on;
xlabel('$e_\mathrm{z}$ [-]','Interpreter','Latex');
ylabel('Normalized frequency');
legend('500','5000','50000');

figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
        eval(['[x,y]=getXYfromVar(exz',num2str(i),',0.05);']);
        stairs(x, y, 'LineWidth', 1);
        hold on;
end
box on;grid on;
xlabel('$e_\mathrm{xz}$ [-]','Interpreter','Latex');
ylabel('Normalized frequency');
legend('500','5000','50000');

figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
        eval(['[x,y]=getXYfromVar(Vxr',num2str(i),',0.05);']);
        stairs(x, y, 'LineWidth', 1);
        hold on;
end
box on;grid on;
xlabel('$U_\mathrm{xcal}$ [m/s]','Interpreter','Latex');
ylabel('Normalized frequency');
legend('500','5000','50000');


% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
%         eval(['[x,y]=getXYfromVar(ez',num2str(i),',0.05);']);
%         stairs(x, y, 'LineWidth', 1);
%         hold on;
% end
% box on;grid on;
% xlabel('$e_\mathrm{z}$ [-]','Interpreter','Latex');
% ylabel('Normalized frequency');
% title('$\Theta=0.05$','Interpreter','Latex');
% legend('$\Omega$=0','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex','fontsize',10);

% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
%         eval(['[x,y]=getXYfromVar(exz',num2str(i),',0.05);']);
%         stairs(x, y, 'LineWidth', 1);
%         hold on;
% end
% box on;grid on;
% xlabel('$e_\mathrm{xz}$ [-]','Interpreter','Latex');
% ylabel('Normalized frequency');
% legend('$\Omega$=0','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex','fontsize',10);



%time-varying values
for i=4:6
eval(['n_dis=ceil(length(t',num2str(i),')/N_inter);']);
eval(['t_inter',num2str(i),' = linspace((t',num2str(i),'(1)+t',num2str(i),'(n_dis))*0.5,(t',num2str(i),'(end)+t',num2str(i),'(end-n_dis+1))*0.5,N_inter);']);
end
% figure
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% for i=0:4
% eval(['plot(t_inter, ez_t',num2str(i),');']);
% hold on;
% end
% xlabel('t [s]');
% ylabel('$\bar{e}_\mathrm{z}$ [$\%$]','Interpreter','Latex','fontsize',12);
% legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);

figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=4:6
eval(['plot(t_inter',num2str(i),', RIM',num2str(i),');']);
hold on;
end
xlabel('t [s]');
ylabel('$\bar{R}_\mathrm{impact}$ [kg/m$^2$/s]','Interpreter','Latex','fontsize',12);
legend('500','5000','50000');
%legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);


%Restitution coefficient and particle liquid content change
ezmean = [];exmean = [];dLmean = [];
for i=0:4
    eval(['[ez',num2str(i),',ex',num2str(i),',dL',num2str(i),',Vim',num2str(i),',thetaim',num2str(i),',thetare',num2str(i),',Eim',num2str(i),',Smean',num2str(i),']=portDataToVectorSal(Par',num2str(i),',Vx',num2str(i),',Vz',num2str(i),',X',num2str(i),',Z',num2str(i),',Vp',num2str(i),',LCm',num2str(i),',N_inter);']);
end

for i=0:4
eval(['ezmean(',num2str(i),'+1)=mean(ez',num2str(i),');']);
eval(['exmean(',num2str(i),'+1)=mean(ex',num2str(i),');']);
eval(['exstd(',num2str(i),'+1)=std(ex',num2str(i),');']);
eval(['ezstd(',num2str(i),'+1)=std(ez',num2str(i),');']);
eval(['dLmean(',num2str(i),'+1)=mean(dL',num2str(i),');']);
eval(['dLstd(',num2str(i),'+1)=std(dL',num2str(i),');']);
eval(['Vimmean(',num2str(i),'+1)=mean(Vim',num2str(i),');']);
eval(['Vimstd(',num2str(i),'+1)=std(Vim',num2str(i),');']);
eval(['thetaimmean(',num2str(i),'+1)=mean(thetaim',num2str(i),');']);
eval(['thetaimstd(',num2str(i),'+1)=std(thetaim',num2str(i),');']);
eval(['thetaremean(',num2str(i),'+1)=mean(thetare',num2str(i),');']);
eval(['Eimmean(',num2str(i),'+1)=mean(Eim',num2str(i),');']);
eval(['Eimstd(',num2str(i),'+1)=std(Eim',num2str(i),');']);
% eval(['Elossmean(',num2str(i),'+1)=mean(Eim',num2str(i),num2str(j),'.*(1-ex',num2str(i),num2str(j),'));']);
% eval(['Elossstd(',num2str(j),',',num2str(i),'+1)=std(Eim',num2str(i),num2str(j),'.*(1-ex',num2str(i),num2str(j),'));']);
end

%distribution
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
titles = {'$\Omega$=0','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$'};
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=1:5
    eval(['[x,y]=getXYfromVar(ex',num2str(i-1),',0.05);']);
    stairs(x, y, 'LineWidth', 1);
    hold on;
    %xlim([-1 2]);
    %ylim([0 0.4]);
    box on;grid on;
    xlabel('$e_\mathrm{xz}$ [-]','Interpreter','Latex');
end
legend('$\Omega$=0$','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex');


figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=1:5
    eval(['[x,y]=getXYfromVar(thetaim',num2str(i-1),',5);']);
    stairs(x, y, 'LineWidth', 1);
    hold on;
%     xlim([0 80]);
%     ylim([0 0.5]);
    box on;grid on;
    xlabel('$\theta_\mathrm{im}$ [$\circ$]','Interpreter','Latex');
end
legend('$\Omega$=0$','$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex');


colors = [0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
titles = {'Eroding phase','Depositing phase','Steady state'};
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=1:4
    eval(['[x,y]=getXYfromVar(dL',num2str(i),');']);
    stairs(x, y, 'LineWidth', 1);
    hold on;
end
xlim([-20 10]);ylim([0 0.8]);
box on;grid on;
xlabel('$\Omega_\mathrm{p,r}-\Omega_\mathrm{p,i}$ [$\%$]','Interpreter','Latex','fontsize',12);
ylabel('Normalized frequency');
legend('$\Omega$=1$\%$','$\Omega$=5$\%$','$\Omega$=10$\%$','$\Omega$=20$\%$','Interpreter','Latex');

%mean and std
colors = [1 0 0;0 0 1;0 0 0];

figure
subplot(1,2,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    %plot(omega,thetaremean(j,:),'o--');
errorbar(omega,Vimmean(j,:),Vimstd(j,:),'^--');
hold on;
end
box on;
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);%ylabel('$\bar{\theta}_\mathrm{re}$ [$\circ$]','Interpreter','Latex','fontsize',12);
ylabel('$\bar{U}_\mathrm{im}$ [m/s]','Interpreter','Latex','fontsize',12);
subplot(1,2,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    plot(omega,thetaimmean(j,:),'diamond--');
    %errorbar(omega,thetaimmean(j,:),thetaimstd(j,:),'^--');
hold on;
end
box on;
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\bar{\theta}_\mathrm{im}$ [$\circ$]','Interpreter','Latex','fontsize',12);
legend('Eroding phase','Depositing phase','Steady state'); 

figure
subplot(1,2,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    %plot(omega,ezmean(j,:),'^--');
errorbar(omega,ezmean(j,:),ezstd(j,:),'^--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
end
box on;
ylim([0 1]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\bar{e}_\mathrm{z}$ [-]','Interpreter','Latex','fontsize',12);
legend('Erosion-dominated phase','Deposition-dominated phase','Steady state'); 
subplot(1,2,2)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    %plot(omega,exmean(j,:),'^--');
errorbar(omega,exmean(j,:),exstd(j,:),'^--','LineWidth', 1.5, 'MarkerSize', 8);
hold on;
end
box on;
ylim([0 1]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\bar{e}_\mathrm{xz}$ [-]','Interpreter','Latex','fontsize',12);


figure
%subplot(1,2,1)
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    %plot(omega,thetaremean(j,:),'o--');
errorbar(omega,Eimmean(j,:),Eimstd(j,:),'^--');
hold on;
end
box on;
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);
ylabel('$\bar{E}_\mathrm{im}$ [kgm$^2$/s$^2$]','Interpreter','Latex','fontsize',12);
legend('Eroding phase','Depositing phase','Steady state'); 

figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=1:3
    %plot(omega,thetaremean(j,:),'o--');
errorbar(omega,Elossmean(j,:),Elossstd(j,:),'^--');
hold on;
end
box on;
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);
ylabel('$\bar{E}_\mathrm{loss}$ [kgm$^2$/s$^2$]','Interpreter','Latex','fontsize',12);
legend('Eroding phase','Depositing phase','Steady state'); 

%values versus time
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
t_cg=linspace(0.01,5,501);
n_dis=500/N_inter;
t_inter = linspace((t_cg(1)+t_cg(1+n_dis))*0.5,(t_cg(end)+t_cg(end-n_dis))*0.5,(length(t_cg)-1)/n_dis);
figure
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i=0:4
eval(['plot(t_inter, ex',num2str(i),');']);
hold on;
end
xlabel('t [s]');
ylabel('$\bar{\Omega}_\mathrm{p,E}$ [$\%$]','Interpreter','Latex','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);


