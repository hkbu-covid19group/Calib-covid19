clear;

%% load data
load './Sampledata/IncubationPeriod.mat';

%% Maximum likelihood estimation
t_e=6;
fun = @(para)IP_func_ln_c(para, t_e, all_s, all_l);
[para, Lval]= fminsearch(fun,[0.1, 0.1]);

%% Display estimated parameters
disp(['Natural logarithm of likelihood=' num2str(-Lval)]);
mu=para(1);
disp(['mu=' num2str(mu)]);
sigma=para(2);
disp(['sigma=' num2str(sigma)]);
disp(['t_e=' num2str(t_e)]);
gamma=-1*lognpdf_d(t_e,mu,sigma)/lognpdf(t_e,mu,sigma);
disp(['gamma=' num2str(gamma)]);

%% Figure
figure('position',[138 156 590 480]);hold on;

tt=0:0.1:30;
plot(tt,P_onset(tt,[para, gamma],t_e),'linewidth',1.5,'color','r');

box on
axis square
set(gca,'xlim',[0,20]);
set(gca,'ylim',[0.001,1]);
set(gca,'fontsize',18)
ax=gca;
ax.YScale='log';
set(gca,'TickLength',[0.03,0.03]);

xlabel('Days after infection $t$','fontsize',22,'Interpreter','latex')
ylabel('Symptom onset time distribution $P_{\mbox{O}}(t)$','fontsize',22,'Interpreter','latex')

%% Functions
function L = IP_func_ln_c( para, t_e, IPl, IPu )

L=0;
mu=para(1);
sigma=para(2);
gamma=-1*lognpdf_d(t_e,mu,sigma)/lognpdf(t_e,mu,sigma);

fun=@(t)P_onset(t,[mu,sigma,gamma],t_e);

for i=1:length(IPl)
    if (IPl(i)-0.5<t_e) && (IPu(i)+0.5>t_e)
        temp1=integral(fun,IPl(i)-0.5,t_e);
        temp2=integral(fun,t_e,IPu(i)+0.5);
        temp=temp1+temp2;
    else
        temp=integral(fun,IPl(i)-0.5,IPu(i)+0.5);
    end

    L=L+log(temp);
end

L=-L;

end

function Po=P_onset(t, para, t_e)
mu=para(1);
sigma=para(2);
gamma=para(3);

temp1=erf((log(t_e)-mu)/(sigma*2^0.5));
temp2=2*lognpdf(t_e,mu,sigma)/gamma;
A=2/(1+temp1+temp2);

Po(t<=t_e)=A*lognpdf(t(t<=t_e),mu,sigma); 
Po(t>t_e)=A*lognpdf(t_e,mu,sigma)*exp(-gamma*(t(t>t_e)-t_e));

end

function y = lognpdf_d(x,mu,sigma)
%Derivative of lognorm pdf

temp1=-1*(log(x)-mu+sigma.^2)./(x.^2.*sigma.^3.*sqrt(2*pi));
temp2=exp(-0.5 * ((log(x) - mu)./sigma).^2);
y=temp1.*temp2;

end






