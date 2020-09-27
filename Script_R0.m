clear;

%% Parameters
% Incubation period distribution
mu=1.51263165945187;
t_e=6;
sigma=0.578861784111941;
gamma=0.305502668928836;

% Infectiousness
t_p=-0.677265707789004;
alpha_A=0.434131422318435;
alpha_B=0.540626680919779;


%% R0
lambda=flip(0.6:-0.01:-gamma);
R0=fun_R0_lambda(lambda,[mu,sigma,gamma],t_e,[alpha_A,alpha_B,t_p]);

%R0=0 at lambda=-gamma
lambda=[-gamma lambda];
R0=[0 R0];


%% R0 at lambda=0.3
R0_lambda03=R0(lambda==0.3);
disp(['R_0 at lambda=0.3 is ' num2str(R0_lambda03)]);

%% Figure
figure('position',[138 156 590 480]);hold on;

plot([0 5], [0.3 0.3],'--','linewidth',0.5,'color','k');
plot([R0_lambda03 R0_lambda03], [-0.5 0.5],'--','linewidth',0.5,'color','r');
plot([0 5], [-gamma -gamma],'--','linewidth',0.5,'color','b');

plot(R0,lambda,'linewidth',1.5,'color','r');

box on
axis square

set(gca,'xlim',[0,5]);
set(gca,'ylim',[-0.5,0.5]);
set(gca,'fontsize',18)
ax=gca;
set(gca,'TickLength',[0.03,0.03]);

xlabel('Effective reproduction number $R_{\mbox{E}}$','fontsize',22,'Interpreter','latex')
ylabel('Daily growth rate $\lambda$ (1/day)','fontsize',22,'Interpreter','latex')

%% Functions
function R0 = fun_R0_lambda(lambda,para_IP,t_e,para_infec)

alpha_A=para_infec(1);
alpha_B=para_infec(2);
theta_p=-para_infec(3);
theta_s=theta_p-alpha_A/alpha_B/(alpha_A+alpha_B);

R0=[];
for j=1:length(lambda)
    fun=@(t)P_onset_ln_R0(t,para_IP,t_e,lambda(j));
    temp1=integral(fun,0,t_e);
    temp2=integral(fun,t_e,Inf);
    R0(j)=exp(-1*lambda(j)*theta_s)/(temp1+temp2);
end

end

function Po=P_onset_ln_R0(t, para, t_e, lambda)
mu=para(1);
sigma=para(2);
gamma=para(3);

temp1=erf((log(t_e)-mu)/(sigma*2^0.5));
temp2=2*lognpdf(t_e,mu,sigma)/gamma;

A=2/(1+temp1+temp2);

if t<=t_e
    Po=A*lognpdf(t,mu,sigma).*exp(-lambda.*t);    
else
    Po=A*lognpdf(t_e,mu,sigma).*exp(gamma*t_e).*exp(-(gamma+lambda)*t);
end

end







