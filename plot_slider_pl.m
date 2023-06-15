set(0, 'DefaultLineLineWidth', 1.5);
f3 = figure(1003);
k = 1;
h3 = plot(x_bat,eta_pl(k,:),'bx'); 
% h1 = plot( num2cell(x_bat), num2cell([phie_int(k,:);phie_py_int(k,:)]),'x'); 
hold on
plot(x_bat,0*x_bat,'--k')
hold off
xlabel('Normalized Thickness')
ylabel('Potential(V)')
title('\eta_{pl} ')
% legend('No Film Res','W/Film Res','Location','SouthWest')
axis([0 0.44 -0.2 0.2]) 
% [81,10,419,20]
b3 = uicontrol('Parent',f3,'Style','slider','Position',[520,50,20,340],...
              'value',k, 'min',1, 'max',length(eta_pl));

addlistener(b3,'ContinuousValueChange',@(hObject, event) makeplot_ve(hObject, event,h3,eta_pl,x_bat));

f4 = figure(1004);
k = 1;
h4 = plot(x_bat,j_pl(k,:),'bx'); 
% h1 = plot( num2cell(x_bat), num2cell([phie_int(k,:);phie_py_int(k,:)]),'x'); 
hold on
plot(x_bat,0*x_bat,'--k')
hold off
xlabel('Normalized Thickness')
ylabel('Current Density (A/m^2)')
title('j_{pl} ')
% legend('No Film Res','W/Film Res','Location','SouthWest')
axis([0 0.44 -1e-7 0]) 
% [81,10,419,20]
b4 = uicontrol('Parent',f4,'Style','slider','Position',[520,50,20,340],...
              'value',k, 'min',1, 'max',length(j_pl));

addlistener(b4,'ContinuousValueChange',@(hObject, event) makeplot_ve(hObject, event,h4,j_pl,x_bat));


function makeplot_ve(hObject,event,h,Ve,x_bat)
n = get(hObject,'Value');
k = floor(n);
% set(h,'ydata',[Ve(k,:);Ve_int(k,:)]);
set(h,{'xdata'},{x_bat},{'YData'}, {Ve(k,:)});
drawnow;
end