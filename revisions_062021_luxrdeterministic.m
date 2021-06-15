%% semi-analytical
Ka = 5; % molecules
f = 460; % fold change
K = 100; %molecules
Rtmin = 5; % molecules/cell: 24 for weakest RBS
alp = Rtmin*3.1e-4;
n = 2;
kdil = 3.1e-4;
A = 0;
kr = 0.001; kf1 = (kr)/Ka;

pars = struct('K',K','f',f,'alp',alp, 'kf1',kf1, 'kr',kr,'At',A);
ahlrange = logspace(-2,4,100);
for i = 1:length(ahlrange)
    pars.At = ahlrange(i)*40;
    [t,R] = ode15s(@luxsimpleode,[0 50]*3600, [alp/kdil 0],{},pars);
    ss(i,:) = R(end,:);
end
figure(8); 
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
loglog(ahlrange, ss(:,1)/ss(1,1),'color',rand(1,3)); hold on;%col(3,:)) ;hold on; 
xlabel('AHL (a.u.)'); ylabel('GFP (fold)')
xlim([0.02 20])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
'LineWidth', 2,'layer','top');
run luxrdata
loglog(B34(:,1)*1000,B34(:,2)/B34(1,2),'o','markersize',6,'linestyle','none','color',col(1,:));hold on
loglog(B64(:,1)*1000,B64(:,2)/B64(1,2),'o','markersize',6,'linestyle','none','color',col(2,:)); 
loglog(B32(:,1)*1000,B32(:,2)/B32(1,2),'o','markersize',6,'linestyle','none','color',col(3,:))
loglog(B31(:,1)*1000,B31(:,2)/B31(1,2),'o','markersize',6,'linestyle','none','color',col(4,:))
