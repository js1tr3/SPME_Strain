
%%
set(0, 'DefaultLineLineWidth', 1.5);
figure(1)
clf
ax=[];
ax(1) = subplot(3,1,[2 3]);

plot(t,Vt)
ylabel('Voltage (V)')
xlabel('time (s)')
title('Voltage')

ax(2) = subplot(311);
plot(t,u)

ylabel('Current (A)')
title('Control')
linkaxes(ax,'x')
% print('fig1','-djpeg')

figure(2)
clf
plot(x_bat,ce(end,:),'x')

xlabel('Normalized Thickness')
ylabel('Conc (mol/m^3)')
title('Electrolyte Concentration')
% print('fig6','-djpeg')

%%
figure(3)
clf
plot(x_bat,phie(1,:),'x')

xlabel('Normalized Thickness')
ylabel('Potential(V)')
% ylim([-0.44 0])
title('\phi_e Start of Sim')
% print('fig7','-djpeg')


figure(4)
clf
plot(x_bat,phie(end,:),'x')
xlabel('Normalized Thickness')
ylabel('Potential(V)')
title('\phi_e End of Sim')
xlim([0 0.44])
% ylim([-0.44 0])
% print('fig8','-djpeg')

figure(5)
clf
plot(x_bat,eta_pl(end,:),'x')
hold on
plot(x_bat,0*x_bat,'--k')
hold off
xlabel('Normalized Thickness')
ylabel('Potential(V)')
title('\eta_{pl} End of Sim')
xlim([0 0.44])

figure(6)
clf
plot(x_bat,j_pl(end,:),'x')
hold on
plot(x_bat,0*x_bat,'--k')
hold off
xlabel('Normalized Thickness')
ylabel('Current Density (mol/m^2s)')
title('j_{pl} End of Sim')
xlim([0 0.44])

return
figure(7)
clf
plot(t,Li_loss_rate)
xlabel('Time')
ylabel('Lithium Loss Rate in mol/s')
title('Lithium Loss Rate')

figure(8)
clf
plot(t,Li_loss)
xlabel('Time')
ylabel('Lithium Loss in mol')
title('Lithium Loss')

figure(9)
clf
plot(t,Q_loss)
xlabel('Time')
ylabel('Q Loss in Ah')
title('Q Loss')



