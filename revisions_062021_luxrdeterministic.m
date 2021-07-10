%% semi-analytical
f = 460; % fold change
Rtmin = 2.5; % molecules/cell
Ka = 5; % molecules
K = 100; %molecules
n = 2;
kdil = 3.1e-4;
A = 0;
kr = 0.001; kf1 = (kr)/Ka;

phi = 1;
ahlrange = logspace(-4,4,100); % nM % convert to molecules
% for i = 1:length(ahlrange)
%     A = ahlrange(i);
% err = @(Ra,Rt) (Ra-0.5*((phi*Ka + phi*A + Rt) - sqrt((phi*Ka + phi*A + Rt)^2-4*phi*A*Rt)))^2 + (Rt- Rtmin*(1+f*(Ra/K)^n)/(1+(Ra/K)^n))^2; % error function: perhaps superfluous, SSE for Ra and Rt both.
% Gfp(i,:) = fmincon(@(x) err(x(1),x(2)),[A f*Rtmin],[],[],[],[],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f]); % second column of fmincon output is Rt, normalized Rt transfer function is proportional to GFP transfer function
% % alternate optimization/simultaneous nonlinear equation solver
% % err2 = @(Ra,Rt) [Ra-0.5*((Ka + A + Rt) - sqrt((Ka + A + Rt)^2-4*A*Rt)) (Rt- Rtmin*(1+f*(Ra/K)^2)/(1+(Ra/K)^2))];
% % Gfp(i,:) = lsqnonlin(@(x) err2(x(1),x(2)),[A f*Rtmin],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f])
% end

% plot solutions (normalized) against data
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
run luxrdata % data file
figure(8); %loglog(ahlrange, Gfp(:,2)/Gfp(1,2),'linewidth',1.5,'color',[0 0 1]); hold on; xlabel('AHL'); ylabel('normalized R_T'); title('dose response')
loglog(B34(:,1)*1000,B34(:,2)/B34(1,2),'o','markersize',6,'linestyle','none','color',col(1,:));hold on
loglog(B64(:,1)*1000,B64(:,2)/B64(1,2),'o','markersize',6,'linestyle','none','color',col(2,:)); 
loglog(B32(:,1)*1000,B32(:,2)/B32(1,2),'o','markersize',6,'linestyle','none','color',col(3,:))
loglog(B31(:,1)*1000,B31(:,2)/B31(1,2),'o','markersize',6,'linestyle','none','color',col(4,:))

pars = struct('K',K','f',f,'alp',Rtmin*kdil, 'kf1',kf1, 'kr',kr,'At',A);

for i = 1:length(ahlrange)
    pars.At = ahlrange(i)*1160/15; % ~1160 == kdif/kpd scaling factor * 0.6 nM --> molecules conversion
    [t,R] = ode15s(@luxsimpleode,[0 50]*3600, [Rtmin 0],{},pars);
    ss(i,:) = R(end,:);
end
figure(8); 
loglog(ahlrange, ss(:,1)/ss(1,1),'color',rand(1,3)); hold on;%col(3,:)) ;hold on; 
xlabel('AHL (a.u.)'); ylabel('GFP (fold)')
xlim([0.02 20])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
'LineWidth', 2,'layer','top');