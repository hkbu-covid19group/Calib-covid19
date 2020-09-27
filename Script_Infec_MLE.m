clear;

%% load data
load './Sampledata/ExposureWindow.mat';

%% Maximum likelihood estimation
fun = @(para)ML_infec_ew(para,W_l,W_u);
[para, Lval]= fminsearch(fun,[0.1, 0.1, -0.1]);

%% Display estimated parameters
disp(['Natural logarithm of likelihood=' num2str(-Lval)]);
alpha_A=para(1);
disp(['alpha_A=' num2str(alpha_A)]);
alpha_B=para(2);
disp(['alpha_B=' num2str(alpha_B)]);
t_P=-para(3);
disp(['theta_P=' num2str(-t_P)]);
theta_s=-t_P-alpha_A/alpha_B/(alpha_A+alpha_B);
disp(['theta_S=' num2str(theta_s)]);

%% Figure
%figure('position',[138 156 200 160]);hold on;
figure('position',[138 156 590 480]);hold on;

tt=-8:0.01:6;
plot(tt,P_infec(tt,alpha_A,alpha_B,t_P),'linewidth',1.5,'color','r');

box on
axis square
set(gca,'xlim',[-8,6]);
set(gca,'ylim',[0,0.3]);
set(gca,'fontsize',18)
ax=gca;
set(gca,'TickLength',[0.03,0.03]);

xlabel('Days after symptom onset $t$','fontsize',22,'Interpreter','latex')
ylabel('Normalized infectiousness $P_{\mbox{I}}(t)$','fontsize',22,'Interpreter','latex')

%% Functions
function L = ML_infec_ew ( para, W_l, W_u )

alpha_A=para(1);
alpha_B=para(2);
t_p=para(3);

L=0;
fun=@(t)P_infec(t,alpha_A,alpha_B,t_p);

for i=1:length(W_l)
    if (W_l(i)-0.5<t_p) && (W_u(i)+0.5>t_p)
        temp1=integral(fun,W_l(i)-0.5,t_p);
        temp2=integral(fun,t_p,W_u(i)+0.5);
        temp=temp1+temp2;
    else
        temp=integral(fun,W_l(i)-0.5,W_u(i)+0.5);
    end

    L=L+log(temp); 
end

L=-L;

end

function Ps=P_infec(t,alpha_A,alpha_B,t_p)
    
    A=alpha_A.*alpha_B./(alpha_A+alpha_B);
    Ps(t<=t_p)=A*exp(alpha_A*(t(t<=t_p)-t_p));
    Ps(t>t_p)=A*exp(-1*alpha_B*(t(t>t_p)-t_p));
    
end








