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

% R0
R0=3.87;

% Contact tracing on day 4 after infection
day_ct=4;

%% R0 matrix

R0_ct_100=R0R_trace(day_ct,1,[alpha_A,alpha_B,t_p],[mu,sigma,gamma],t_e);
R0_mat_day4=[];
for i=0:1:1000
    for j=0:1:1000
        p_ct=i*0.001;
        p_mask=j*0.001;
        temp=R0*(1-R0_ct_100*p_ct);
        R0_mat_day4(i+1,j+1)=temp*(1-0.5*p_mask)^2;
    end
end

% Contour for R0=1
cl_r01_day4=[];
for i=1:1:1001
    cl_r01_day4(i)=find(R0_mat_day4(i,:)<=1,1);
end

%% Figure
figure('position',[138 156 590 480]);hold on;
cl_map=cl_map_gen(min(min(R0_mat_day4)),max(max(R0_mat_day4)),500);
imagesc(R0_mat_day4);colormap(cl_map);hold on;
plot(cl_r01_day4-1,0:1:1000,'--','linewidth',1,'color','k');

box on
axis square
set(gca,'xlim',[0,1000]);
set(gca,'ylim',[0,1000]);
set(gca,'xtick',0:200:1000)
set(gca,'ytick',0:200:1000)
set(gca,'xticklabels',{'0%','20%','40%','60%','80%','100%'});
set(gca,'yticklabels',{'0%','20%','40%','60%','80%','100%'});
set(gca,'fontsize',18)
ax=gca;
set(gca,'TickLength',[0.03,0.03]);
xlabel('Percentage of mask-wearing population','fontsize',22,'Interpreter','latex')
ylabel('Percentage of traced contacts','fontsize',22,'Interpreter','latex')
title('Basic reproduction number $R_e$','fontsize',22,'Interpreter','latex')

c = colorbar;
c.Label.FontSize=18;

legend('$R_0=1.0$','Location','northwest','fontsize',22,'Interpreter','latex');
legend('boxoff');

%% Functions
function R0R = R0R_trace(trace_t,pct,para_infec,para_IPD,t_e)

alpha_A=para_infec(1);
alpha_B=para_infec(2);
theta_p=-para_infec(3);
theta_s=theta_p-alpha_A/alpha_B/(alpha_A+alpha_B);

fun=@(t)P_onset_t(t,para_IPD,t_e);
if trace_t-theta_s<=t_e
    qs=integral(fun,0,trace_t-theta_s);
else
    temp1=integral(fun,0,t_e);
    temp2=integral(fun,t_e,trace_t-theta_s);
    qs=temp1+temp2;
end

R0R=pct*(1-qs);

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

function cl_map = cl_map_gen (R_min,R_max,N_cl)
% R=R_min at (0 0 0.5)
% R=1 at (0 1 1)
% R=2 at (1 1 0)
% R=3 at (1 0 0)
% R=R_max at (0 0 0)

index1=round((1-R_min)/(R_max-R_min)*N_cl);
index2=round(1/(R_max-R_min)*N_cl)+index1;
index3=round(1/(R_max-R_min)*N_cl)+index2;

cl_map=zeros(N_cl,3);
%R_0=R_min
cl_map(1:index1,1)=0;
cl_map(1:index1,2)=linspace(0,1,index1);
cl_map(1:index1,3)=linspace(0.5,1,index1);
%R_0=1
cl_map(index1:index2,1)=linspace(0,1,index2-index1+1);
cl_map(index1:index2,2)=1;
cl_map(index1:index2,3)=linspace(1,0,index2-index1+1);
%R_0=2;
cl_map(index2:index3,1)=1;
cl_map(index2:index3,2)=linspace(1,0,index3-index2+1);
cl_map(index2:index3,3)=0;
%R_0=3;
cl_map(index3:N_cl,1)=linspace(1,0,N_cl-index3+1);
cl_map(index3:N_cl,2)=0;
cl_map(index3:N_cl,3)=0;

end












