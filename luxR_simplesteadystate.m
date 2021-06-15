%% semi-analytical
Ka = 408.8; % molecules
f = 550; % fold change
K = 2.93; %molecules
Rtmin = 10*0.45*0.01/3.1e-4; % K*2/(f+1)
alp = Rtmin*3.1e-4;
n = 2;

Rt1 = (10.^(-2:0.1:4))*Rtmin;
Ra2 = 10.^(-4:0.1:4); Rt2 = Rtmin*(1+f*(Ra2/K).^n)./(1+(Ra2/K).^n);
figure(7);  loglog(Rt2,Ra2,'r'); hold on; xlabel('R_T'); ylabel('R_a'); title('Modules')
r = logspace(-2,3.7,10);
% for i = 1:length(r)
         A = 5; %r(i);
         Ra1 = 0.5*((Ka + A + Rt1) - sqrt((Ka + A + Rt1).^2-4*A*Rt1));
         loglog(Rt1, Ra1,'b');
% end

tic
ahlrange = logspace(-2,4,50);
for i = 1:length(ahlrange)
    A = ahlrange(i);
err = @(Ra,Rt) (Ra-0.5*((Ka + A + Rt) - sqrt((Ka + A + Rt)^2-4*A*Rt)))^2 + (Rt- Rtmin*(1+f*(Ra/K)^n)/(1+(Ra/K)^n))^2;
Gfp(i,:) = fmincon(@(x) err(x(1),x(2)),[A f*Rtmin],[],[],[],[],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f]);
% err2 = @(Ra,Rt) [Ra-0.5*((Ka + A + Rt) - sqrt((Ka + A + Rt)^2-4*A*Rt)) (Rt- Rtmin*(1+f*(Ra/K)^2)/(1+(Ra/K)^2))];
% Gfp(i,:) = lsqnonlin(@(x) err2(x(1),x(2)),[A f*Rtmin],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f])
end
toc

col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
run luxrdata
figure(8); loglog(ahlrange, Gfp(:,2)/Gfp(1,2)); hold on; xlabel('AHL'); ylabel('normalized R_T'); title('dose response')
loglog(B34(:,1)*1000,B34(:,2)/B34(1,2),'o','markersize',1.75,'linestyle','none','color',col(1,:));hold on
loglog(B64(:,1)*1000,B64(:,2)/B64(1,2),'o','markersize',1.75,'linestyle','none','color',col(2,:)); 
loglog(B32(:,1)*1000,B32(:,2)/B32(1,2),'o','markersize',1.75,'linestyle','none','color',col(3,:))
loglog(B31(:,1)*1000,B31(:,2)/B31(1,2),'o','markersize',1.75,'linestyle','none','color',col(4,:))


loglog(B33(:,1)*1000,B33(:,2)/B33(1,2),'o','markersize',1.75)


Rt1 = 1*Rtmin;
xT = logspace(-2,4,25); n = 2;
alpha= (Ka+xT+Rt1);
beta = (alpha - sqrt(alpha.^2-4*xT.*Rt1));
rart = Rt1.*(1-(alpha-2*xT)./(sqrt(alpha.^2-4*xT.*Rt1)))./beta;
Ra1 = 0.5*((Ka + xT + Rt1) - sqrt((Ka + xT + Rt1).^2-4*xT*Rt1));
y = Ra1/K;
fra = (2*(f-1)*y.^2)./((1+y.^2).*(1+f*y.^2));

figure(4); loglog(xT/Ka,rart,'b'); hold on; loglog(xT/Ka, fra,'r'); loglog(xT/Ka, ones(1,length(xT)),'k')
legend('L(R_a,R_T)','L(F,R_a)'); xlabel('a_T/K_a'); ylabel('Log Gain')
figure(5); loglog(xT/Ka, rart.*fra); hold on; loglog(xT/Ka, ones(1,length(xT)),'k');
xlabel('a_T/K_a'); ylabel('L(R_a,R_T).L(F,R_a)'); title('product of log gains')

%% simple ODE model optimization
          % Ka   f   K    Rmin
x0 = log10([100  420 5  25]);
lb = log10([10   415 1   1]);
ub = log10([1000 425 100 100]);
% options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter');
% Y = fmincon(@luxerr_simplesteadystate,x0,[],[],[],[],lb,ub,[],options2)
options = optimoptions('particleswarm','SwarmSize',10,'PlotFcn',@pswplotbestf,'Display','iter');
options.MaxIterations = 150; options.FunctionTolerance = 1e-3;
for num = 1:5
   [sol,fval]= particleswarm(@luxerr_simplesteadystate,4,lb,ub,options);
   X(num,:)  = [sol,fval]
end
X = [0.900864407874771,2.61584667887231,-0.878643904047058+1,0.282558582817605,-1.59079155635523+1]% B31
X = [0.8780    2.6232    0.0857    0.2996   -0.5818];
% optimizing b32,64,34
lb = ([0    0   0]);
ub = ([10  10 10]);
options2 = optimoptions('particleswarm','SwarmSize',10,'PlotFcn',@pswplotbestf,'Display','iter');
options2.MaxIterations = 150; options2.FunctionTolerance = 1e-3;
X2 = particleswarm(@(g) luxerr_2(g,X),3,lb,ub,options2);
%% simple ODE
% Ka = 100; f = 420; K = 5; alp = 1*kdil;
% X = [0.900864407874771,2.61584667887231,-0.878643904047058+1,0.282558582817605,-1.59079155635523+1]; scl = 1;% B31
% X = [2.83174017598842,2.61922673582138,0.687894778856033,1.74137911914448]; %B31#2; scaled up, nM units
nav = 6.023e23*1e-15; scl = nav*1e-9; kdil = 3.1e-4;
% Ka = 10^X(1); f = 10^X(2); 
% K =10^X(3); alp = 1*10^X(4)*kdil; 
A = 0;
kr = 0.001; kf1 = (kr)/Ka;

pars = struct('K',K','f',f,'alp',alp, 'kf1',kf1, 'kr',kr,'At',A);
ahlrange = logspace(-2,4,100);
for i = 1:length(ahlrange)
    pars.At = ahlrange(i);
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
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
'LineWidth', 1,'layer','top');
%% simbio model
% X = [0.900864407874771,2.61584667887231,-0.878643904047058+1,-1.59079155635523+1]; scl = 1;% B31
X = [2.83174017598842,2.61922673582138,0.687894778856033,1.74137911914448]; %B31#2; scaled up, nM units
nav = 6.023e23*1e-15; scl = nav*1e-9;
Ka = 10^X(1)*scl; f = 10^X(2); K =10^X(3)*scl; alp = 15*10^X(4)*3.1e-4*scl; 

ktlnG = 0.08; kpd = 3.1e-4; kmd = 2e-3; 
ktr = 0.44*3.1e-4; ktlnR = alp*kmd/ktr; kprb = 1e-2; kprf = kprb/(K)^2; kdis = 0.001; kra = kdis/(Ka);
m1 = sbiomodel('luxr');
Robj3 = addreaction(m1,'m -> m + G');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'ktlnG','Value',ktlnG,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'ktlnG'};
Robj3 = addreaction(m1,'m -> m + R');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'ktlnR','Value',ktlnR,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'ktlnR'};
Robj3 = addreaction(m1,'p0 + 2 Ra <-> p1');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kprf','Value',kprf,'ValueUnits','1/(molecule^2*second)');
    Pobj3r = addparameter(Kobj3,'kprb','Value',kprb,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kprf','kprb'};
Robj3 = addreaction(m1,'p0 -> p0 + m');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'ktr','Value',ktr,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'ktr'};
Robj3 = addreaction(m1,'p1 -> p1 + m');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'ktract','Value',ktr*f,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'ktract'};
Robj3 = addreaction(m1,'R + a <-> Ra');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kra','Value',kra,'ValueUnits','1/(molecule*second)');
    Pobj3r = addparameter(Kobj3,'kdis','Value',kdis,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kra','kdis'};
Robj3 = addreaction(m1,'m -> null');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kmd','Value',kmd,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kmd'};    
Robj3 = addreaction(m1,'R -> null');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kpd','Value',kpd,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kpd'};    
Robj3 = addreaction(m1,'G -> null');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kpd','Value',kpd,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kpd'};    
Robj3 = addreaction(m1,'Ra -> a');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kpd','Value',kpd,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kpd'};
Robj3 = addreaction(m1,'p1 -> p0 + 2 a');
    Kobj3 = addkineticlaw(Robj3,'MassAction');
    Pobj3 = addparameter(Kobj3,'kpd','Value',kpd,'ValueUnits','1/second');
    Kobj3.ParameterVariableNames = {'kpd'};
    
for i = 1:size(m1.Species,1)
    m1.Species(i).InitialAmountUnits = 'molecule';
end

m1.species(4).InitialAmount = 1; % promoter copy #
m1.species(7).InitialAmount = 0; % ahl copies
m2 = m1;
%% simbiomodel simulations: stochastic
% % % % % stochastic setup % % % % % 
% % % % % % % % % % % % % % % % % % % % 
configsetObj = getconfigset(m1,'active');
configsetObj.SolverType = 'ssa';
configsetObj.StopTime = 3600*55;
configsetObj.SolverOptions.LogDecimation = 1;
configsetObj.RunTimeOptions.StatesToLog = 'all';
% % % pre-ahl % % %
% simdata_pre = sbiosimulate(m1,configsetObj);
simdata_pre = sbioensemblerun(m1,2,configsetObj);
for num = 1:2; dj(num) = simdata_pre(num,1).Data(end,2); end; figure; hj = histogram(dj); hj.BinWidth = 1; hj.Normalization = 'pdf';
hold on; x=0:1:100; plot(x,gampdf(x,0.5, ktlnR/kmd),'r')
std(dj)/mean(dj)
figure; plot(simdata_pre(2,1).Time/3600, simdata_pre(2,1).Data(:,2))
% initial condition for post-ahl
for i = 1:size(m1.Species,1)
    m1.Species(i).InitialAmount = simdata_pre(1,1).Data(end,i);
end

configsetObj.SolverOptions.LogDecimation = 10;
configsetObj.StopTime = 3600*15;
% figure;
nR = 10;
ahlrange = logspace(-1,1,20);
for pt = 1:length(ahlrange)
m1.species(7).InitialAmount = ahlrange(pt);
% simdata_post = sbiosimulate(m1,configsetObj);
simdata_post = sbioensemblerun(m1,nR,configsetObj);
    figure
    stochss = [];
    for i = 1:nR
%             tind = find(simdata_post(i,1).Time >= 12*3600,1,'first');
            stochss = cat(1,stochss, simdata_post(i,1).Data(end,:));
            plot(simdata_post(i,1).Time/3600, simdata_post(i,1).Data(:,2)); hold on
    end
%     histogram(stochss(:,2),'BinWidth',1);
    meanval(pt) = mean(stochss(:,2));
    cvval(pt) = std(stochss(:,2))/mean(stochss(:,2));
end
    col = [1 0.75 0
            0 0 0 
            0 0.5 1
            0 0.5 0];
    figure(2); subplot(1,2,1); loglog(ahlrange,meanval,'s--','color',col(1,:),'linewidth',2); hold on;
    xlabel('AHL'); ylabel('GFP (mean)'); title('dose response');
    subplot(1,2,2); semilogx(ahlrange,cvval,'s--','color',col(1,:),'linewidth',2'); hold on;
    xlabel('AHL'); ylabel('CV'); title('noise');
    suptitle('stochastic simulation; simple model; 100 runs')
% % % % deterministic sim % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
ahlrange = logspace(-2,2,50);
m2config = getconfigset(m2,'active');
m2config.SolverType = 'ode15s';
m2config.StopTime = 3600*15;
kbas = m2.reactions(2).KineticLaw(1).Parameters.Value;

m2.reactions(2).KineticLaw(1).Parameters.Value = kbas*1.0;
for i = 1:size(m2.Species,1)
    m2.Species(i).InitialAmount = 0;
end
m2.species(4).initialamount = 1;
preahl = sbiosimulate(m2);
for i = 1:size(m2.Species,1)
    m2.Species(i).InitialAmount = preahl.Data(end,i);
end

for k = 1:length(ahlrange)
m2.species(7).InitialAmount = ahlrange(k);
odesim = sbiosimulate(m2);
gfp(k) = odesim.Data(end,2);
end
figure(3); loglog(ahlrange,gfp/gfp(1)); hold on;