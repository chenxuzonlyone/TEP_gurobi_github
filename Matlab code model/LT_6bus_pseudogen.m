clear
clc
cd('/Users/zhangcaihua/Desktop/TEP_gurobi/data')
% all use pu convert theta(voltage angle) to pu(rad)
%% basic information 
%from/to/cost/reactance/capability  
line_info = xlsread('mj_6bus_cost.xlsx',1);% all line
%from/to/cost/reactance/capability  
Eline_info = xlsread('mj_6bus_cost.xlsx',2);%existing line
%from/to/cost/reactance/capability  
Cline_info = xlsread('mj_6bus_cost.xlsx',3);%candidate line
%bus_num/gen_max/gen_actual
gen_info = xlsread('mj_6bus_cost.xlsx',4);
%bus_num/load
load_info = xlsread('mj_6bus_cost.xlsx',5);
%bus number info
bus_info = xlsread('mj_6bus_cost.xlsx',6);

%% time record start
tic;% measure time for execution time
%% relax indices information                                    (control)
% relax variable collection. 0= relax; 1=don't relax
relax_slack_variable1=0;%0 regard as no slack variable%%%%%%%%%%%%%%%%%%%%%%%
relax_slack_variable2=0;%0 regard as no slack variable%%%%%%%%%%%%%%%%%%%%%%%
relax_gen_variable=1; %0 regard gen as known variables; 1 regard as unknown
relax_candidateline_variable=1;   %empty matrix means no candidate line
                                  %1 means candidate lines exist
relax_existingline_variable=[];   %empty matrix means no existing line
                                  %1 means existing lines exist
load_change=1;%                                               (control)

%% given parameter 
% indices information
line_Efrom=1*relax_existingline_variable;%existing line info indices
line_Eto=2*relax_existingline_variable;
line_Ecost=3*relax_existingline_variable;
line_Ereactance=4*relax_existingline_variable;
line_Elimit=5*relax_existingline_variable;
line_Cfrom=1*relax_candidateline_variable;% Candidate line info indices
line_Cto=2*relax_candidateline_variable;
line_Ccost=3*relax_candidateline_variable;
line_Creactance=4*relax_candidateline_variable;
line_Climit=5*relax_candidateline_variable;
gen_busnum=1;% gen info indices (generation at bus #)
gen_min=2;
gen_max=3;
gen_fixed=4;
gen_cost=5;
load_busnum=1;% load info indices
load_fixed=2;

% basic number info
nbus=bus_info;%                                                 (control)
nEline=size(Eline_info,1);
nCline=size(Cline_info,1);
nGen=size(gen_info,1);
nload=size(load_info,1);
nMLC=1;%nMaxLineConnect                                         (control)
nCline_total=nEline*(nMLC-1)+nCline*nMLC;
nplan_year=5;% numbers of planning years                       (control)
%Fkl_0+Fkl_1               Xkl                          I  Gk   S1 S2 thata
%Fkl_0,Fkl_1,              Xkl,                            Gk,  Rk1Rk2theta,     Gk_ts, I]
%Aeq_NB_t= [flow_NB_t,      x_NB_t, -Kp_NB_t, s1_NB_t,-s2_NB_t, theta_NB_t, -Kp_ts_NB_t, I_NB_t];%generation reagrded as unknown variables
nAV=(nEline*nMLC+nCline*nMLC)+(nEline*(nMLC-1)+nCline*nMLC)*2+(nGen)+(nbus*3)+nCline_total; %numbers of  annual variable 

%basic gen & load info 
gen_P=gen_info(:,gen_fixed); %%%%%%%%%%%%%%%fixed generation output //not include in the C code
load_P=load_info(:,load_fixed)*load_change;
load_increase_factor=0.01;%                                       (control)
discount_rate=0.1; % discount rate
M=7*10^4;%                                                         (control)
M_kl=7*10^4;% candidate branch flow different M meanning same value
loss_penalty=10^4;
ref_position=1; % reference bus position
duration_time=8760; %year duration time (hour) 365*24*60
capacity_factor=0.6;% capacity factor of generator 50%
million_transfer=10^6;% transfer unit from dollar to million dollar
candidate_line_pool=[1:8]';
max_open_num=0;%nCline: all lines are controlable. 1: 1 line is controable (control) // not included in the c file

%% parameter calculation
%bus related incidence matrix
% bus branch incidence matrix
Kl_E=zeros(nbus,nEline);
Kl_C=zeros(nbus,nCline);
%existing 
for i=1:nEline
    Kl_E(Eline_info(i,line_Efrom),i)= 1;
    Kl_E(Eline_info(i,line_Eto  ),i)=-1;
end
%candidate
for i=1:nCline
    Kl_C(Cline_info(i,line_Cfrom),i)= 1;
    Kl_C(Cline_info(i,line_Cto  ),i)=-1;
end
 % bus generation incidence matrix
 Kp=zeros(nbus,nGen);
 for i=1:nGen
    Kp(gen_info(i,gen_busnum),i)= 1;
 end
% bus pseudo generator incidence matrix (ts: transmission switching)
Kp_ts=zeros(nbus,nCline);
for i=1:nCline
    Kp_ts(Cline_info(i,line_Cfrom),i)= 1;
    Kp_ts(Cline_info(i,line_Cto  ),i)=-1;
end
% bus load incidence matrix
 Kd=zeros(nbus,nload);
 for i=1:nload
    Kd(load_info(i,load_busnum),i)= 1;
 end
% bus slack variable incidence matrix
Kr1=eye(nbus);
Kr2=eye(nbus);
% bus angle(theta) incidence matrix
Ka=eye(nbus); 
% branch reactance matrix
X_E=diag(Eline_info(:,line_Ereactance));
X_C=diag(Cline_info(:,line_Creactance));

%% nodal balance constraints
Aeq_NB=zeros(nplan_year*nbus,nplan_year*nAV);
beq_NB=zeros(nplan_year*nbus,1);
for i=1:nplan_year

load_Pt=load_P*(1+load_increase_factor)^(i-1);% load will increase annually according to increase factor

flow_NB_t=[Kl_E,repmat(Kl_E,1,nMLC-1),repmat(Kl_C,1,nMLC)];% include exist & candidate flow (non_zero)(checked)
x_NB_t=[repmat(zeros(nbus,nEline),1,nMLC-1),repmat(zeros(nbus,nCline),1,nMLC)];%decision variable (zero)(checked)
gen_NB_t=Kp*gen_P;%%%%%%%%% regarded as known variables (zero)
Kp_NB_t=Kp*(relax_gen_variable-0);%(non_zero)(checked)
load_NB_t=Kd*load_Pt;%%(non_zero)
s1_NB_t=Kr1*(relax_slack_variable1-0);%slack variable 1 %%%%%%%% this constraint can be relaxed (zero) (checked)
s2_NB_t=Kr2*(relax_slack_variable2-0);%slack variable 2 %%%%%%%% this constraint can be relaxed (zero) (checked)
theta_NB_t=Ka*0;%(zero) (checked)
Kp_ts_NB_t=Kp_ts;% pseudo generator(non_zero) (checked)
I_NB_t=zeros(nbus,size(x_NB_t,2));% Installation status change variable. When installation status change, it equals to 1; otherwise it equals to 0; (zero)

%variab [   Fkl_0,Fkl_1,    Xkl,    Gk,         Rk1,   Rk2,      theta,     Gk_ts,       I]
Aeq_NB_t= [flow_NB_t,      x_NB_t, -Kp_NB_t, s1_NB_t,-s2_NB_t, theta_NB_t, -Kp_ts_NB_t, I_NB_t];%generation reagrded as unknown variables
beq_NB_t=-load_NB_t+gen_NB_t*(1-relax_gen_variable);

% annual update
Aeq_NB(1+(i-1)*nbus:nbus+(i-1)*nbus, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aeq_NB_t; % update annually 
beq_NB(1+(i-1)*nbus:nbus+(i-1)*nbus,1)=beq_NB_t; % update annually
end

%% status change constraints
Aeq_SC=zeros(nplan_year*nCline_total,nplan_year*nAV);
beq_SC=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
flowE_SC_t=zeros(nCline_total,nEline);%(zero) (existing line setting is Not included in the C code)
flowC_SC_t=zeros(nCline_total,nCline_total);%(zero) (checked)
x_SC_C_t= eye(nCline_total);% current year status (non_zero) (checked)   
x_SC_P_t= eye(nCline_total);% previous year status (non_zero) (checked)
gen_SC_t=zeros(nCline_total,nGen);%(zero) (checked)
s1_SC_t=zeros(nCline_total,nbus);%(zero) (checked)
s2_SC_t=zeros(nCline_total,nbus);%(zero) (checked)
theta_SC_t=zeros(nCline_total,nbus);%(zero) (checked)
Kp_ts_SC_t=zeros(nCline_total,nCline_total);% pseudo generator (zero)
I_SC_C_t=eye(nCline_total);% current year status change (non_zero) (checked)
I_SC_P_t=zeros(nCline_total,nCline_total);% previous year status change (zero)
%variab     [Fkl_0,     Fkl_1,              Xkl,   Gk,      Rk1,        Rk2,   theta,     Gk_ts         I]
Aeq_SC_C_t= [flowE_SC_t,flowC_SC_t,      x_SC_C_t, gen_SC_t,s1_SC_t, s2_SC_t, theta_SC_t, Kp_ts_SC_t -I_SC_C_t];
beq_SC_C_t= zeros(nCline_total,1);

% annual update
Aeq_SC(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aeq_SC_C_t; % update annually 
beq_SC(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=beq_SC_C_t; % update annually

% previous year update(Assuming initial status of every line is 0, first year status change equal to first year status)
if i>1
    %variab [    Fkl_0,     Fkl_1,       Xkl(Prev),Gk,      Rk1,     Rk2,     theta,      Gk_ts     I]
    Aeq_SC_P_t= [flowE_SC_t,flowC_SC_t, -x_SC_P_t, gen_SC_t,s1_SC_t, s2_SC_t, theta_SC_t, Kp_ts_SC_t I_SC_P_t];
    Aeq_SC(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-2)*nAV:nAV+(i-2)*nAV)=Aeq_SC_P_t; % previous year update 
end
end
% Aeq_SC=[];
% beq_SC=[];
%% Installation status maintains constraints. (Inequality)
Aineq_SM=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_SM=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
flowE_SM_t=zeros(nCline_total,nEline); %(zero) (not included in the C model)
flowC_SM_t=zeros(nCline_total,nCline_total); %(zero) (checked)
x_SM_C_t= eye(nCline_total);% current year status (nonzero) (checked)
x_SM_P_t= eye(nCline_total);% previous year status (nonzero) (checked)
gen_SM_t=zeros(nCline_total,nGen); %(zero) (checked)
s1_SM_t=zeros(nCline_total,nbus);%(zero) (checked)
s2_SM_t=zeros(nCline_total,nbus);%(zero) (checked)
theta_SM_t=zeros(nCline_total,nbus);%(zero) (checked) 
Kp_ts_SM_t=zeros(nCline_total,nCline_total);% pseudo generator (zero) (checked) 
I_SM_t=zeros(nCline_total,nCline_total); %(zero) 

%variab [Fkl_0,         Fkl_1,            Xkl,      Gk,      Rk1,     Rk2,     theta,      Gk_ts,       I]
Aineq_SM_C_t= [flowE_SM_t,flowC_SM_t,    -x_SM_C_t, gen_SM_t,s1_SM_t, s2_SM_t, theta_SM_t, Kp_ts_SM_t I_SM_t]; 
bineq_SM_C_t= zeros(nCline_total,1);

% annual update
Aineq_SM(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_SM_C_t; % update annually 
bineq_SM(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_SM_C_t; % update annually

% previous year update(Assuming initial status of every line is 0, first year status change equal to first year status)
if i>1
    %variab [        Fkl_0,     Fkl_1,     Xkl(Prev),Gk,      Rk1,     Rk2,     theta,      Gk_ts,     I]
    Aineq_SM_P_t= [flowE_SM_t,flowC_SM_t,  x_SM_P_t, gen_SM_t,s1_SM_t, s2_SM_t, theta_SM_t, Kp_ts_SM_t I_SM_t];
    Aineq_SM(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-2)*nAV:nAV+(i-2)*nAV)=Aineq_SM_P_t; % previous year update 
end
end
% Aineq_SM=[];
% bineq_SM=[];
%% power flow

% existing power flow
Aeq_EPF=zeros(nplan_year*nEline,nplan_year*nAV);
beq_EPF=zeros(nplan_year*nEline,1);
for i=1:nplan_year
flowE_EPF_t=eye(nEline);
flowC_EPF_t=zeros(nEline,nEline*(nMLC-1)+nCline*nMLC);
x_EPF_t=zeros(nEline,nEline*(nMLC-1)+nCline*nMLC);
gen_EPF_t=zeros(nEline,nGen);
s1_EPF_t=zeros(nEline,nbus);
s2_EPF_t=zeros(nEline,nbus);
theta_EPF_t=(X_E^-1)*Kl_E';
Kp_ts_EPF_t=zeros(nEline,nCline_total);% pseudo generator
I_EPF_t=zeros(nEline,nCline_total);

%variab [Fkl_0,         Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,         Rk1,   Rk2,     theta,       Gk_ts,      I ]
Aeq_EPF_t=[flowE_EPF_t, flowC_EPF_t,           x_EPF_t,gen_EPF_t,s1_EPF_t,s2_EPF_t,-theta_EPF_t, Kp_ts_EPF_t,I_EPF_t];
beq_EPF_t=zeros(nEline,1);

% annual update
Aeq_EPF(1+(i-1)*nEline:nEline+(i-1)*nEline, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aeq_EPF_t;
beq_EPF(1+(i-1)*nEline:nEline+(i-1)*nEline,1)=beq_EPF_t;
end

% candidate power flow (Positive)
Aineq_P_CPF=zeros(nplan_year*nCline_total,nplan_year*nAV);% positive 
bineq_P_CPF=zeros(nplan_year*nCline_total,1);
Aineq_N_CPF=zeros(nplan_year*nCline_total,nplan_year*nAV);% negative
bineq_N_CPF=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
flowE_CPF_t=zeros(nCline_total,nEline); %(zero) (not include in the C model)
flowC_CPF_t=eye(nCline_total);%(nonzero) (checked)
x_CPF_t=eye(nCline_total)*M;%(nonzero) (checked)
gen_CPF_t=zeros(nCline_total,nGen);%(zero) (checked)
s1_CPF_t=zeros(nCline_total,nbus);%(zero) (checked)
s2_CPF_t=zeros(nCline_total,nbus);%(zero) (checked)
theta_CPF_t=zeros(nCline_total,nbus);%(nonzero) (checked)
Kp_ts_CPF_t=zeros(nCline_total,nCline_total);% pseudo generator %(zero) (checked)
I_CPF_t=zeros(nCline_total,nCline_total);%(zero) (checked)

for j=1:nMLC-1 %existing line's candidate (not include in the C model)
    theta_CPF_t(1+(j-1)*nEline:nEline+(j-1)*nEline,:)=(X_E^-1)*Kl_E';% X^-1*Kl=F
end
for j=1:nMLC %pure candidate lines
    theta_CPF_t(1+nEline*(nMLC-1)+(j-1)*nCline:nCline+nEline*(nMLC-1)+(j-1)*nCline,:)=(X_C^-1)*Kl_C';
end

% candidate power flow (Positive)
%    variab [   Fkl_0,      Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,         Rk1,   Rk2,     theta,      Gk_ts,          I]
Aineq_P_CPF_t=[flowE_CPF_t, flowC_CPF_t,           x_CPF_t,gen_CPF_t,s1_CPF_t,s2_CPF_t,-theta_CPF_t,Kp_ts_CPF_t, I_CPF_t];% Positive
bineq_P_CPF_t=ones(nCline_total,1)*M;

% candidate power flow (Negative)
%    variab [   Fkl_0,      Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,         Rk1,   Rk2,      theta,       Gk_ts,       I]
Aineq_N_CPF_t=[ flowE_CPF_t,-flowC_CPF_t,           x_CPF_t,gen_CPF_t,s1_CPF_t,s2_CPF_t, theta_CPF_t,Kp_ts_CPF_t, I_CPF_t];% Negative
bineq_N_CPF_t=ones(nCline_total,1)*M;

% annual update
Aineq_P_CPF(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_P_CPF_t;
bineq_P_CPF(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_P_CPF_t;
Aineq_N_CPF(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_N_CPF_t;
bineq_N_CPF(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_N_CPF_t;
end
% Aineq_P_CPF=[];
% bineq_P_CPF=[];
% Aineq_N_CPF=[];
% bineq_N_CPF=[];
%% branch limit
% existing (not included in the C model)
Aineq_P_EBL=zeros(nplan_year*nEline,nplan_year*nAV);
bineq_P_EBL=zeros(nplan_year*nEline,1);
Aineq_N_EBL=zeros(nplan_year*nEline,nplan_year*nAV);
bineq_N_EBL=zeros(nplan_year*nEline,1);

for i=1:nplan_year
% existing branch limit (Positive)
flowE_EBL_t=eye(nEline);
flowC_EBL_t=zeros(nEline,nCline_total);
x_EBL_t=zeros(nEline,nCline_total);
gen_EBL_t=zeros(nEline,nGen);
s1_EBL_t=zeros(nEline,nbus);
s2_EBL_t=zeros(nEline,nbus);
theta_EBL_t=zeros(nEline,nbus);
Kp_ts_EBL_t=zeros(nEline,nCline_total);% pseudo generator
I_EBL_t=zeros(nEline,nCline_total);

%    variab [ Fkl_0,        Fkl_1(Fkl_1E,Fkl_1C),   Xkl,  Gk,       Rk1,    Rk2,        theta,        Gk_ts,        I]
Aineq_P_EBL_t=[ flowE_EBL_t, flowC_EBL_t,           x_EBL_t,gen_EBL_t,s1_EBL_t,s2_EBL_t, theta_EBL_t, Kp_ts_EBL_t, I_EBL_t];% Positive
bineq_P_EBL_t=Eline_info(:,line_Elimit);

% existing branch limit (Negative)
%    variab [ Fkl_0,        Fkl_1(Fkl_1E,Fkl_1C),   Xkl,  Gk,       Rk1,    Rk2,        theta,      Gk_ts,      I]
Aineq_N_EBL_t=[-flowE_EBL_t, flowC_EBL_t,           x_EBL_t,gen_EBL_t,s1_EBL_t,s2_EBL_t, theta_EBL_t,Kp_ts_EBL_t,I_EBL_t];% Negative
bineq_N_EBL_t=Eline_info(:,line_Elimit);

% annual update
Aineq_P_EBL(1+(i-1)*nEline:nEline+(i-1)*nEline, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_P_EBL_t;
bineq_P_EBL(1+(i-1)*nEline:nEline+(i-1)*nEline,1)=bineq_P_EBL_t;
Aineq_N_EBL(1+(i-1)*nEline:nEline+(i-1)*nEline, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_N_EBL_t;
bineq_N_EBL(1+(i-1)*nEline:nEline+(i-1)*nEline,1)=bineq_N_EBL_t;
end

% candidate
Aineq_P_CBL=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_P_CBL=zeros(nplan_year*nCline_total,1);
Aineq_N_CBL=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_N_CBL=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
% candidate branch limit (Positive)
flowE_CBL_t=zeros(nCline_total,nEline); %(zero) (not include in the C model)
flowC_CBL_t=eye(nCline_total);%(nonzero) (checked)
x_CBL_t=diag([repmat(Eline_info(:,line_Elimit)',1,nMLC-1),repmat(Cline_info(:,line_Climit)',1,nMLC)])-eye(nCline_total)*M; %(nonzero) (checked)
gen_CBL_t=zeros(nCline_total,nGen);%(zero)(checked)
s1_CBL_t=zeros(nCline_total,nbus);%(zero)(checked)
s2_CBL_t=zeros(nCline_total,nbus);%(zero)(checked)
theta_CBL_t=zeros(nCline_total,nbus);%(zero)(checked)
Kp_ts_CBL_t=zeros(nCline_total,nCline_total);% pseudo generator (zero)(checked)
I_CBL_t=zeros(nCline_total,nCline_total);%(zero)(checked)

%    variab [Fkl_0,         Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,        Rk1,    Rk2,        theta,    Gk_ts,      I]
Aineq_P_CBL_t=[flowE_CBL_t, flowC_CBL_t,          -x_CBL_t,gen_CBL_t,s1_CBL_t,s2_CBL_t,theta_CBL_t, Kp_ts_CBL_t, I_CBL_t]; %(Positive)
bineq_P_CBL_t=ones(nCline_total,1)*M;

% candidate branch limit (Negative)
%    variab [Fkl_0,         Fkl_1(Fkl_1E,Fkl_1C),Xkl,    Gk,        Rk1,    Rk2,        theta,      Gk_ts,       I]
Aineq_N_CBL_t=[flowE_CBL_t,-flowC_CBL_t,          -x_CBL_t,gen_CBL_t,s1_CBL_t,s2_CBL_t,theta_CBL_t, Kp_ts_CBL_t, I_CBL_t]; %(Negative)
bineq_N_CBL_t=ones(nCline_total,1)*M;

% annual update
Aineq_P_CBL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_P_CBL_t;
bineq_P_CBL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_P_CBL_t;
Aineq_N_CBL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_N_CBL_t;
bineq_N_CBL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_N_CBL_t;
end
% Aineq_P_CBL=[];
% bineq_P_CBL=[];
% Aineq_N_CBL=[];
% bineq_N_CBL=[];
%% generation limit
Aineq_R_GL=zeros(nplan_year*nGen,nplan_year*nAV);
bineq_R_GL=zeros(nplan_year*nGen,1);
Aineq_L_GL=zeros(nplan_year*nGen,nplan_year*nAV);
bineq_L_GL=zeros(nplan_year*nGen,1);

for i=1:nplan_year
flowE_GL_t=zeros(nGen,nEline); %(zero) (not included in the C model)
flowC_GL_t=zeros(nGen,nCline_total);%(zero)(checked)
x_GL_t=zeros(nGen,nCline_total);%(zero)(checked)
gen_GL_t=eye(nGen);%(nonzero)(checked)
s1_GL_t=zeros(nGen,nbus);%(zero)(checked)
s2_GL_t=zeros(nGen,nbus);%(zero)(checked)
theta_GL_t=zeros(nGen,nbus);%(zero)(checked)
Kp_ts_GL_t=zeros(nGen,nCline_total);% pseudo generator (zero)(checked)
I_GL_t=zeros(nGen,nCline_total);%(zero)(checked)

%right
% variab   [    Fkl_0,     Fkl_1(Fkl_1E,Fkl_1C),    Xkl,   Gk,      Rk1,        Rk2,    theta,   Gk_ts,     I]
Aineq_R_GL_t=[flowE_GL_t,  flowC_GL_t,            x_GL_t, gen_GL_t, s1_GL_t, s2_GL_t, theta_GL_t,Kp_ts_GL_t,I_GL_t];
bineq_R_GL_t= gen_info(:,gen_max);

%left
% variab   [  Fkl_0,       Fkl_1(Fkl_1E,Fkl_1C),    Xkl,   Gk,        Rk1,   Rk2,       theta,     Gk_ts,   I]
Aineq_L_GL_t=[flowE_GL_t,  flowC_GL_t,            x_GL_t,-gen_GL_t, s1_GL_t, s2_GL_t, theta_GL_t,Kp_ts_GL_t,I_GL_t];
bineq_L_GL_t=-gen_info(:,gen_min);

% annual update
Aineq_R_GL(1+(i-1)*nGen:nGen+(i-1)*nGen, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_R_GL_t;
bineq_R_GL(1+(i-1)*nGen:nGen+(i-1)*nGen,1)=bineq_R_GL_t;
Aineq_L_GL(1+(i-1)*nGen:nGen+(i-1)*nGen, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_L_GL_t;
bineq_L_GL(1+(i-1)*nGen:nGen+(i-1)*nGen,1)=bineq_L_GL_t;
end
% Aineq_R_GL=[];
% bineq_R_GL=[];
% Aineq_L_GL=[];
% bineq_L_GL=[];
%% slack limit (already included in the upper & lower bound. So not show in the C program)
% slack 1
Aineq_S1L=zeros(nplan_year*nbus,nplan_year*nAV);
bineq_S1L=zeros(nplan_year*nbus,1);

for i=1:nplan_year
flowE_S1L_t=zeros(nbus,nEline);
flowC_S1L_t=zeros(nbus,nCline_total);
x_S1L_t=zeros(nbus,nCline_total);
gen_S1L_t=zeros(nbus,nGen);
s1_S1L_t=eye(nbus);
s2_S1L_t=zeros(nbus,nbus);
theta_S1L_t=zeros(nbus,nbus);
Kp_ts_S1L_t=zeros(nbus,nCline_total);% pseudo generator
I_S1L_t=zeros(nbus,nCline_total);

% variab  [Fkl_0,           Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,         Rk1,    Rk2,      theta,       Gk_ts,       I]
Aineq_S1L_t=[flowE_S1L_t, flowC_S1L_t,            x_S1L_t, gen_S1L_t, s1_S1L_t, s2_S1L_t, theta_S1L_t, Kp_ts_S1L_t I_S1L_t];
bineq_S1L_t=inf*ones(nbus,1);

% annual update
Aineq_S1L(1+(i-1)*nbus:nbus+(i-1)*nbus, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_S1L_t;
bineq_S1L(1+(i-1)*nbus:nbus+(i-1)*nbus,1)=bineq_S1L_t;
end

% slack 2
Aineq_S2L=zeros(nplan_year*nbus,nplan_year*nAV);
bineq_S2L=zeros(nplan_year*nbus,1);

for i=1:nplan_year
flowE_S2L_t=zeros(nbus,nEline);
flowC_S2L_t=zeros(nbus,nCline_total);
x_S2L_t=zeros(nbus,nCline_total);
gen_S2L_t=zeros(nbus,nGen);
s1_S2L_t=zeros(nbus,nbus);
s2_S2L_t=eye(nbus);
theta_S2L_t=zeros(nbus,nbus);
Kp_ts_S2L_t=zeros(nbus,nCline_total);% pseudo generator
I_S2L_t=zeros(nbus,nCline_total);

% variab  [Fkl_0,           Fkl_1(Fkl_1E,Fkl_1C),   Xkl,    Gk,         Rk1,    Rk2,      theta,       Gk_ts,       I]
Aineq_S2L_t=[flowE_S2L_t, flowC_S2L_t,            x_S2L_t, gen_S2L_t, s1_S2L_t, s2_S2L_t, theta_S2L_t, Kp_ts_S2L_t, I_S2L_t];
bineq_S2L_t=inf*ones(nbus,1);

% annual update
Aineq_S2L(1+(i-1)*nbus:nbus+(i-1)*nbus, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_S2L_t;
bineq_S2L(1+(i-1)*nbus:nbus+(i-1)*nbus,1)=bineq_S2L_t;
end
Aineq_S1L=[];
bineq_S1L=[];
Aineq_S2L=[];
bineq_S2L=[];
%% switch line number limit
% % flowE_SLN=zeros(1,nEline);
% % flowC_SLN=zeros(1,nCline_total);
% % x_SLN=ones(1,nCline_total);
% % gen_SLN=zeros(1,nGen);
% % s1_SLN=zeros(1,nbus);
% % s2_SLN=zeros(1,nbus);
% % theta_SLN=zeros(1,nbus);
% % 
% % %variab [     Fkl_0,     Fkl_1,     Xkl,   Gk,      Rk1,    Rk2,    theta    ]
% % Aineq_P_SLN= [flowE_SLN, flowC_SLN,-x_SLN, gen_SLN, s1_SLN, s2_SLN, theta_SLN];%generation reagrded as unknown variables
% % bineq_P_SLN=-nCline+max_open_num;

Aineq_P_SLN=[];% get rid of max OFF line limit
bineq_P_SLN=[];% get rid of max OFF line limit

%% pseudo generator line relationship
Aineq_P_PS=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_P_PS=zeros(nplan_year*nCline_total,1);
Aineq_N_PS=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_N_PS=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
flowE_PS_t=zeros(nCline_total,nEline); %(zero) (not included in the C model)
flowC_PS_t=eye(nCline_total);%(nonzero) (checked)
x_PS_t=eye(nCline_total)*M;%(nonzero) (checked)
gen_PS_t=zeros(nCline_total,nGen);%(zero) (checked) 
s1_PS_t=zeros(nCline_total,nbus);%(zero) (checked)
s2_PS_t=zeros(nCline_total,nbus);%(zero) (checked)
theta_PS_t=zeros(nCline_total,nbus);%(zero) (checked)
Kp_ts_PS_t=eye(nCline_total);% pseudo generator (nonzero) (checked)
I_PS_t=zeros(nCline_total,nCline_total);%(zero) (checked)

%positive
%variab     [Fkl_0,      Fkl_1,      Xkl,       Gk,     Rk1,    Rk2,        theta,     Gk_ts,       I]
Aineq_P_PS_t= [flowE_PS_t, flowC_PS_t,-x_PS_t, gen_PS_t, s1_PS_t, s2_PS_t, theta_PS_t, -Kp_ts_PS_t, I_PS_t];%generation reagrded as unknown variables
bineq_P_PS_t=zeros(nCline_total,1);

%negative 
%variab     [Fkl_0,      Fkl_1,      Xkl,       Gk,     Rk1,    Rk2,        theta,     Gk_ts,       I]
Aineq_N_PS_t= [flowE_PS_t,-flowC_PS_t,-x_PS_t, gen_PS_t, s1_PS_t, s2_PS_t, theta_PS_t,  Kp_ts_PS_t, I_PS_t];%generation reagrded as unknown variables
bineq_N_PS_t=zeros(nCline_total,1);

% annual update
Aineq_P_PS(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_P_PS_t;
bineq_P_PS(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_P_PS_t;
Aineq_N_PS(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_N_PS_t;
bineq_N_PS(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_N_PS_t;
end
% Aineq_P_PS=[];
% bineq_P_PS=[];
% Aineq_N_PS=[];
% bineq_N_PS=[];
%% Pesudo generator limit
Aineq_P_PSL=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_P_PSL=zeros(nplan_year*nCline_total,1);
Aineq_N_PSL=zeros(nplan_year*nCline_total,nplan_year*nAV);
bineq_N_PSL=zeros(nplan_year*nCline_total,1);

for i=1:nplan_year
flowE_PSL_t=zeros(nCline_total,nEline); %(zero)(not include in C model)
flowC_PSL_t=zeros(nCline_total,nCline_total);%(zero) (checked)
x_PSL_t=eye(nCline_total)*M;%(nonzero) (checked)
gen_PSL_t=zeros(nCline_total,nGen);%(zero) (checked)
s1_PSL_t=zeros(nCline_total,nbus);%(zero) (checked)
s2_PSL_t=zeros(nCline_total,nbus);%(zero) (checked)
theta_PSL_t=zeros(nCline_total,nbus);%(zero) (checked)
Kp_ts_PSL_t=eye(nCline_total);% pseudo generator (nonzero) (checked)
I_PSL_t=zeros(nCline_total,nCline_total);%(zero) (checked)

%positive
%variab     [   Fkl_0,      Fkl_1,      Xkl,       Gk,      Rk1,        Rk2,        theta,     Gk_ts,       I]
Aineq_P_PSL_t= [flowE_PSL_t, flowC_PSL_t, x_PSL_t, gen_PSL_t, s1_PSL_t, s2_PSL_t, theta_PSL_t, Kp_ts_PSL_t, I_PSL_t];%generation reagrded as unknown variables
bineq_P_PSL_t=ones(nCline_total,1)*M;

%negative
%variab     [   Fkl_0,      Fkl_1,      Xkl,       Gk,      Rk1,        Rk2,        theta,     Gk_ts,       I]
Aineq_N_PSL_t= [flowE_PSL_t, flowC_PSL_t, x_PSL_t, gen_PSL_t, s1_PSL_t, s2_PSL_t, theta_PSL_t,-Kp_ts_PSL_t, I_PSL_t];%generation reagrded as unknown variables
bineq_N_PSL_t=ones(nCline_total,1)*M;

% annual update
Aineq_P_PSL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_P_PSL_t;
bineq_P_PSL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_P_PSL_t;
Aineq_N_PSL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total, 1+(i-1)*nAV:nAV+(i-1)*nAV)=Aineq_N_PSL_t;
bineq_N_PSL(1+(i-1)*nCline_total:nCline_total+(i-1)*nCline_total,1)=bineq_N_PSL_t;
end
% Aineq_P_PSL=[];
% bineq_P_PSL=[];
% Aineq_N_PSL=[];
% bineq_N_PSL=[];

%% Initial status & upper bound & lower bound & ctype (C checked)
lb=zeros(nplan_year*nAV,1);
ub=zeros(nplan_year*nAV,1);
ctype=[];

for i=1:nplan_year
%  set the candidate line status 
line_initial=ones(nCline_total,1);% In fact, this constraint is used to set the limit of binary variable
line_initial(candidate_line_pool)=0;% the meaning is to let the candidate line binary variable limit from 0 to 1

% reference bus position
theta_ref_limit_P= inf*ones(nbus,1);% Positive bus angle limit
theta_ref_limit_P(ref_position,1)=0;% reference position bus angle = 0
theta_ref_limit_N=-inf*ones(nbus,1);% Negative bus angle limit
theta_ref_limit_N(ref_position,1)=0;% reference position bus angle = 0

% variab  [Fkl_0,        Fkl_1(Fkl_1E,Fkl_1C),    Xkl,                  Gk,              Rk1,             Rk2,             theta,               Gk_ts,                      I]
lb_t=[-inf*ones(nEline,1);-inf*ones(nCline_total,1);line_initial;        zeros(nGen,1);   zeros(nbus,1);   zeros(nbus,1);   theta_ref_limit_N; -inf*ones(nCline_total,1); zeros(nCline_total,1)];
ub_t=[ inf*ones(nEline,1); inf*ones(nCline_total,1);ones(nCline_total,1);inf*ones(nGen,1);inf*ones(nbus,1);inf*ones(nbus,1);theta_ref_limit_P;  inf*ones(nCline_total,1); ones(nCline_total,1)];
% variab  [Fkl_0,                Fkl_1(Fkl_1E,Fkl_1C),    Xkl,                      Gk,                 Rk1,                Rk2,             theta,               Gk_ts,                      I]
ctype_t=[repmat('C',1,nEline),repmat('C',1,nCline_total),repmat('B',1,nCline_total),repmat('C',1,nGen),repmat('C',1,nbus),repmat('C',1,nbus),repmat('C',1,nbus), repmat('C',1,nCline_total),  repmat('B',1,nCline_total)];

% annual update
lb(1+(i-1)*nAV:nAV+(i-1)*nAV, 1)=lb_t;
ub(1+(i-1)*nAV:nAV+(i-1)*nAV, 1)=ub_t;
ctype=[ctype,ctype_t];
end

%% all
Aineq=[Aineq_SM;Aineq_P_CPF;Aineq_N_CPF;Aineq_P_EBL;Aineq_N_EBL;Aineq_P_CBL;Aineq_N_CBL;Aineq_R_GL;Aineq_L_GL;Aineq_S1L;Aineq_S2L;Aineq_P_SLN;Aineq_P_PS;Aineq_N_PS;Aineq_P_PSL;Aineq_N_PSL];
bineq=[bineq_SM;bineq_P_CPF;bineq_N_CPF;bineq_P_EBL;bineq_N_EBL;bineq_P_CBL;bineq_N_CBL;bineq_R_GL;bineq_L_GL;bineq_S1L;bineq_S2L;bineq_P_SLN;bineq_P_PS;bineq_N_PS;bineq_P_PSL;bineq_N_PSL];
Aeq=[Aeq_NB;Aeq_SC;Aeq_EPF];
beq=[beq_NB;beq_SC;beq_EPF];

% objective function (C checked)
f=zeros(nplan_year*nAV,1);

for i=1:nplan_year
flowE_F_t=zeros(nEline,1);
flowC_F_t=zeros(nCline_total,1);
x_F_t= zeros(nCline_total,1);
%coeff=(capacity_factor*duration_time)/((1+discount_rate)^(i-1)*million_transfer);
gen_F_t=(gen_info(:,gen_cost)*capacity_factor*duration_time)/((1+discount_rate)^(i-1)*million_transfer);
s1_F_t=(loss_penalty*ones(nbus,1)*duration_time/((1+discount_rate)^(i-1)*million_transfer))*relax_slack_variable1;
s2_F_t=(loss_penalty*ones(nbus,1)*duration_time/((1+discount_rate)^(i-1)*million_transfer))*relax_slack_variable2;
theta_F_t=zeros(nbus,1);
Kp_ts_F_t=zeros(nCline_total,1);% pseudo generator
I_F_t=[repmat(Eline_info(:,line_Ecost),(nMLC-1),1);repmat(Cline_info(:,line_Ccost),nMLC,1)]/((1+discount_rate)^(i-1));

%     [Fkl_0,   Fkl_1(Fkl_1E,Fkl_1C), Xkl,      Gk,         Rk1,    Rk2,   theta,       Gk_ts,      I]
f_t=[flowE_F_t; flowC_F_t;            x_F_t;    gen_F_t;    s1_F_t; s2_F_t;theta_F_t; Kp_ts_F_t; I_F_t];

% annual update
f(1+(i-1)*nAV:nAV+(i-1)*nAV,1)=f_t; % update annually
end

% options = cplexoptimset('cplex');% set the tolerance equal to zero
% options.mip.tolerances.absmipgap = 0;% absolute mip gap = 0
% options.mip.tolerances.mipgap = 0;% mip gap = 0
% options.mip.tolerances.integrality = 0;% mip integrality = 0
% [x_P,fval,exitflag,output]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
[x_P,fval,exitflag,output]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],[]);
timeEnd=toc;% execution time of calculation

fval
exitflag
% % data recording
% x_output=[x_output,x];% data recording
% time_cost=[time_cost;[timeEnd,fval]];
% 
% save('x_output','x_output');
% save('time_cost','time_cost');


