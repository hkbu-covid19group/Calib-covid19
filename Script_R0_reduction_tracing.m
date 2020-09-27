clear;

%% Parameters
% Incubation period distribution
mu=1.51263165945187;
sigma=0.578861784111941;
gamma=0.305502668928836;
t_e=6;

% Infectiousness
t_p=-0.677265707789004;
alpha_A=0.434131422318435;
alpha_B=0.540626680919779;

% percentage of traced contacts
pct=1;

% R0
R0=3.87;

%% R0_reduction
tt=0:0.1:10;

for i=1:length(tt)
    R0_reduction_r(i)=R0R_trace(tt(i),pct,[alpha_A,alpha_B,t_p],...
    [mu,sigma,gamma],t_e);
end

R0_reduction=(1-R0_reduction_r)*R0;

%% Figure
figure('position',[138 156 590 480]);hold on;

plot(tt,R0_reduction,'-','linewidth',1.5,'color','r');

plot([0 10], [R0 R0],'--','linewidth',0.5,'color','k');
plot([0 10], [1 1],'--','linewidth',0.5,'color','k');

box on
axis square
set(gca,'xlim',[0,10]);
set(gca,'ylim',[0,4]);
set(gca,'fontsize',18)
ax=gca;
set(gca,'TickLength',[0.03,0.03]);
xlabel('Contact tracing date after infection $T_{\mbox{trace}}$','fontsize',22,'Interpreter','latex')
ylabel('$R_{\mbox{E}}$ under testing and quarantine','fontsize',22,'Interpreter','latex')
title('Intervention - Contact Tracing','fontsize',22,'Interpreter','latex')

legend('100\% contacts are traced','Position',[.3 .6 .2 .2],'fontsize',16,'Interpreter','latex');
legend('boxoff');

annotation('textbox',[.25 .59 .3 .3],'String','$R_0=3.87$ at $\lambda=0.3$/day',...
    'FitBoxToText','on','Interpreter','latex','fontsize',12);
annotation('textbox',[.25 .005 .3 .3],'String','$R_0=1.0$ at $\lambda=0.0$/day',...
    'FitBoxToText','on','Interpreter','latex','fontsize',12);

%% Functions
function R0R = R0R_trace(trace_t,pct,para_infec,para_IPD,t_e)

if trace_t==0
    R0R=pct;
    return;
end

alpha_A=para_infec(1);
alpha_B=para_infec(2);
theta_p=-para_infec(3);
theta_s=theta_p-alpha_A/alpha_B/(alpha_A+alpha_B);

fun=@(T)qs_T(T,para_IPD,t_e);
temp1=integral(fun,0,trace_t+theta_s);
temp=temp1/trace_t;

R0R=pct*(1-temp);

end


function Po=P_onset_t(t, para, t_e)
mu=para(1);
sigma=para(2);
gamma=para(3);

temp1=erf((log(t_e)-mu)/(sigma*2^0.5));
temp2=2*lognpdf(t_e,mu,sigma)/gamma;

A=2/(1+temp1+temp2);

if t<=t_e
    Po=A*lognpdf(t,mu,sigma);    
else
    Po=A*lognpdf(t_e,mu,sigma)*exp(-gamma*(t-t_e));
end

end

function qs=qs_T(T, para, t_e)

fun=@(t)P_onset_t(t,para,t_e);
qs=zeros(size(T));

parfor i=1:length(T)
    if T(i)<=t_e
        qs(i)=integral(fun,0,T(i));
    else
        temp1=integral(fun,0,t_e);
        temp2=integral(fun,t_e,T(i));
        qs(i)=temp1+temp2;
    end    
end


end









