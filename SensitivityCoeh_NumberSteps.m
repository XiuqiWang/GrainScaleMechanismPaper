M=[0.0940 0.0952 0.0949 0.0898 0.0908 0.0905;
   0.0674 0.0704 0.0692 0.0634 0.0662 0.0650;
   0.0640 0.0652 0.0675 0.0607 0.0616 0.0640;
   0.0606 0.0620 0.0633 0.0573 0.0587 0.0600;
   0.0572 0.0597 0.0617 0.0539 0.0564 0.0584;
   0.0585 0.0603 0.0611 0.0559 0.0573 0.0580];

Ns = [1 2 3];

figure
subplot(3,1,1)
plot(Ns, M(:,1:3),'o--');
xticks([1 2 3]);
xticklabels({'500','5000','50000'});
ylabel('RE [kg/m^2/s]');
subplot(3,1,2)
plot(Ns, M(:,4:6),'o--');
xticks([1 2 3]);
xticklabels({'500','5000','50000'});
ylabel('RD [kg/m^2/s]');
subplot(3,1,3)
plot(Ns, M(:,1:3)-M(:,4:6),'o--');
hold on
plot(Ns, 0.021/5*ones(1,3),'k-');
xticks([1 2 3]);
xticklabels({'500','5000','50000'});
ylim([0.0032 0.0044]);
ylabel('RE-RD [kg/m^2/s]');
xlabel('Number of outputs');
legend('crit = 15D','crit = 20D','crit = 25D','crit = 30D','crit = 40D','Net erosion rate CG');