clear;

%% Parameters
% Incubation period distribution
mu=1.51263165945187;
t_e=6;
sigma=0.578861784111941;
gamma=0.305502668928836;

% Infectiousness
t_p=0.677265707789004;
alpha_A=0.434131422318435;
alpha_B=0.540626680919779;

%% Convolution

tt=-12:0.05:24;
SI_dis=SI_convolution(tt,[mu sigma gamma],t_e,[alpha_A,alpha_B,t_p]);

%% Figure
figure('position',[138 156 590 480]);hold on;

plot(tt,SI_dis,'linewidth',1.5,'color','r');

box on
axis square

set(gca,'ylim',[0.001,0.3]);
set(gca,'xlim',[-12,24]);
set(gca,'fontsize',18)
ax=gca;
ax.YScale='log';
set(gca,'xtick',-12:6:24)
set(gca,'TickLength',[0.03,0.03]);

xlabel('Days after symptom onset of infector $t$','fontsize',22,'Interpreter','latex')
ylabel('Serial interval distribution $P_{SI}(t)$','fontsize',22,'Interpreter','latex')

%% Functions

function dis_si = SI_convolution(SI,para_IP,t_e,para_infec)

mu=para_IP(1);
sigma=para_IP(2);
gamma=para_IP(3);

alpha_A=para_infec(1);
alpha_B=para_infec(2);
t_p=para_infec(3);

dis_si=zeros(size(SI));

for i=1:length(SI)
    fun=@(t)(P_infec1(t,alpha_A,alpha_B,t_p).*P_onset_ln(SI(i)-t, [mu,sigma,gamma], t_e));
    if t_p>=SI(i)
        temp1=integral(fun,-Inf,SI(i)-t_e);
        temp2=integral(fun,SI(i)-t_e,SI(i));
        temp=temp1+temp2;
    elseif t_p>=SI(i)-t_e
        temp1=integral(fun,-Inf,SI(i)-t_e);
        temp2=integral(fun,SI(i)-t_e,t_p);
        temp3=integral(fun,t_p,SI(i));
        temp=temp1+temp2+temp3;
    else
        temp1=integral(fun,-Inf,t_p);
        temp2=integral(fun,t_p,SI(i)-t_e);
        temp3=integral(fun,SI(i)-t_e,SI(i));
        temp=temp1+temp2+temp3;
    end
    dis_si(i)=temp;
end
end

function Ps=P_infec1(t,alpha_A,alpha_B,t_p)
    
    A=alpha_A.*alpha_B./(alpha_A+alpha_B);
    
    Ps(t<=t_p)=A*exp(alpha_A*(t(t<=t_p)-t_p));
    Ps(t>t_p)=A*exp(-1*alpha_B*(t(t>t_p)-t_p));

end

function Po=P_onset_ln(t, para, t_e)
mu=para(1);
sigma=para(2);
gamma=para(3);

temp1=erf((log(t_e)-mu)/(sigma*2^0.5));
temp2=2*lognpdf(t_e,mu,sigma)/gamma;

A=2/(1+temp1+temp2);

Po(t<=t_e)=A*lognpdf(t(t<=t_e),mu,sigma); 
Po(t>t_e)=A*lognpdf(t_e,mu,sigma)*exp(-gamma*(t(t>t_e)-t_e));

end






















