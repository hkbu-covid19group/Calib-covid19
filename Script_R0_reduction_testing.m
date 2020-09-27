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

% testing report delay
tau_d=0;

% R0
R0=3.87;

%% R0 reduction
tt=0:0.1:10;
for i=1:length(tt)
    R0_reduction_r(i)=R0R_test(tt(i),tau_d,[alpha_A,alpha_B,t_p],...
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

xlabel('Testing date after infection $T_{\mbox{test}}$','fontsize',22,'Interpreter','latex')
ylabel('$R_{\mbox{E}}$ under testing and quarantine','fontsize',22,'Interpreter','latex')
title('Intervention - Testing','fontsize',22,'Interpreter','latex')

legend('No testing report delay','Position',[.5 .35 .2 .2],'fontsize',16);
legend('boxoff');

annotation('textbox',[.3 .59 .3 .3],'String','$R_{\mbox{E}}=3.87$ at $\lambda=0.3$/day',...
    'FitBoxToText','on','Interpreter','latex','fontsize',12);
annotation('textbox',[.3 .005 .3 .3],'String','$R_{\mbox{E}}=1.0$ at $\lambda=0.0$/day',...
    'FitBoxToText','on','Interpreter','latex','fontsize',12);

%% Functions
function R0R = R0R_test (test_t,tau_d,para_infec,para_IPD,t_e)

alpha_A=para_infec(1);
alpha_B=para_infec(2);
theta_p=-para_infec(3);

qa1=P_onset_t(test_t+theta_p,para_IPD,t_e)/alpha_A;
temp1=exp(-alpha_A*tau_d)*alpha_B/(alpha_A+alpha_B);
temp2=exp(-alpha_B*tau_d)*alpha_A/(alpha_A+alpha_B);
temp3=(1-exp(-(alpha_A+alpha_B)*tau_d))*alpha_A*alpha_B/(alpha_A+alpha_B)^2;
temp4=exp(-alpha_B*tau_d)*alpha_A*alpha_A/(alpha_A+alpha_B);

fun=@(t)P_onset_t(test_t-t+theta_p,para_IPD,t_e)/alpha_A.*exp(-alpha_B*t);
if test_t+theta_p<=t_e
    temp5=integral(fun,0,test_t);
else
    temp51=integral(fun,0,test_t+theta_p-t_e);
    temp52=integral(fun,test_t+theta_p-t_e,test_t);
    temp5=temp51+temp52;
end

R0R=qa1*(temp1+temp2+temp3)+temp4*temp5;

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










