%% Program to test the 9! groups prioritization

clc
clear all
close all
profile on
tic
%% Parameters:
global b_si b_sy b_ri b_ry v rate_i
global mu_I mu_Y r condition Vacc_limit age_limit test_comb current_group real_order

dI=13; % Infectious period + latent period.
b_si=1/dI;                                              %=beta1   |                          % beta1>=beta2
% b_sy=0.1*b_si;                                        %=beta2   | TO BE DETERMINED         % beta1>=beta3
b_ri=0.01*b_si;  %Denmark estimation 1% reinfection     %=beta3   |                          % beta3>beta4
alpha=0.5;                                              %Half likely for a recovered individual to be infected by a reinfected individual than by a first infected individual.
b_ry=alpha*b_ri;                                        %beta4= beta2*beta3/beta1;           % beta2>beta4
                                                        % beta1/beta2 = beta3/beta4
b_sy=b_ry*b_si/b_ri;
load IFR_9Groups
load Hosp_COVID_Spain_9G                                % % of hospitalizations respect the infections within a given group.
load UCI_COVID_Spain_9G                                 % % of patients requiring UCI respect the infections within a given group.

% IFR_9Groups(1): IFR for individuals 0-9
% IFR_9Groups(2): IFR for individuals 10-19
% ...
% IFR_9Groups(9): IFR for individuals >80


for i=1:length(IFR_9Groups)
    r(i)=(1/dI)*(1-(IFR_9Groups(i)/100));
    mu_I(i)=(1/dI)*(IFR_9Groups(i)/100);
    mu_Y(i)=0;
end

% load CM_2020_9Groups
% load N_2020_9Groups
% M=CM_2020_9Groups; %Matrix 9x9
% N=N_2020_9Groups;

% % To load original contact matrix
% load SIYRD_9G_contactMatrix
% M = SIYRD_9G_contactMatrix;

% To load reduced contact matrix
M = readmatrix('SIYRD_9G_contactMatrix_reduced.csv');

% To load population densities
load N_9x9_2020_ESP
N = N_9x9_2020_ESP;


% load('NAM_SPAIN_2020.mat')

% Para calcular YLL o PYLL (Years of Potential Life Lost)
% Once loaded, the variable's name is NAM. Column 1: 0-85 (ages)           --> We must select the age limit delimiting the G1 and G2
%                                          Column 2: Population with age i.
% Row 1: 0 years..
% Row 80: 79 years
% Row 84: 83 years
% Row 85: >=84 years
c1=1;
c2=10;
for i=1:9
    Av_age(i) = sum(NAM(c1:c2,1).*NAM(c1:c2,2))/sum(NAM(c1:c2,2));                 % Average age G1.
    c1=c2+1;
    c2=c1+9;
    if i==8
        c2=length(NAM);
    end
end
SLE = 83.6;

%% Initial conditions %
I_0=[ 1 1 1 1 1 1 1 1 1];             %Vector containing the initial number of first infected individuals for each group (1 x m)
count=1;
m=length(M);
for i=1:m
%     data_init=[N(i)-I_0 I_0(i) 0 0 0 0 0];
    yd0(count)=N(i)-I_0(i);
    yd0(count+1)=I_0(i);
    yd0(count+2)=0;
    yd0(count+3)=0;
    yd0(count+4)=0;
    yd0(count+5)=0;
    yd0(count+6)=0;
    count=(count+6)+1;
end


%% Solution of the model

% Define vaccination coverage. When the population of susceptible drops below (1-Vacc_limit) it automatically jumps to the other group
Vacc_limit=0.7;

% Define vaccination rates
%rates=[0.05*10^-2 0.1*10^-2 0.5*10^-2 1*10^-2 1.5*10^-2 2*10^-2];
rates=0.05*10^-2;

% Define age groups: G1: 0-9; G2: 10-19; G3: 20-29; G4: 30-39; G5: 40-49; G6: 50-59; G7: 60-69; G8: 70-79
age_limit=8;

for j=1:length(rates)+1
    if j==1 % no vaccination scenario is simulated
        condition=1;
        rate_i=0;
        v=zeros(1,9);
        tmax=365;
        [t,yd] = ode15s(@(t,y)All_combinations_SIYRD_dydn(t,y,M,N),[0,tmax],yd0);
        t_noVacc=t;
        yd_noVacc=yd;
        S_noVacc = 0;
        I_noVacc = 0;
        Y_noVacc = 0;
        R_noVacc = 0;
        D_noVacc = 0;
        CI_noVacc= 0;
        YLL_noVacc=0;
        Hosp_noVacc=0;
        UCI_noVacc=0;

        for t=1:9
            % Now we are no longer interested into group data. Since we aim
            % to evaluate 9! (362,880) vaccination strategies, we will just
            % look at the absolute numbers of deaths and infections. We may
            % analyze them against time or vaccination rates, but we will
            % not continue with the division into the two groups. We will
            % work now with 9 groups and all strategies starting from
            % G1-...-G9 prioritization will be tested.

            S_noVacc = S_noVacc + yd(:,-6+t*7);
            I_noVacc = I_noVacc + yd(:,-5+t*7);
            Y_noVacc = Y_noVacc + yd(:,-4+t*7);
            R_noVacc = R_noVacc + yd(:,-3+t*7);
            D_noVacc = D_noVacc + yd(:,-2+t*7);
            CI_noVacc= CI_noVacc + yd(:,t*7); %Cumulative new infections

            Hosp_noVacc = Hosp_noVacc + max(yd(:,t*7))*Hosp_COVID_Spain_9G(t);  % Total number of patients requiring hospitalization
            UCI_noVacc  = UCI_noVacc + max(yd(:,t*7))*UCI_COVID_Spain_9G(t);    % Total number of patients requiring UCI-hospitalization
            YLL_noVacc= YLL_noVacc + (SLE-Av_age(t))*max(yd(:,-2+t*7));

        end
        Deaths_noVacc=max(D_noVacc);
        Inf_tot_noVacc=max(CI_noVacc);


    else
        condition = 2; % to simulate 9! combinations for each vaccination rate
        rate_i = rates(j-1);

        num_comb = factorial(9);
        ele_comb = 1:9;
        per_comb = perms(ele_comb);
        tmax=365; % 1 year simulation

%        check_comb=1:40320:length(per_comb);
%        check_comb=80641:factorial(8)/8:120960;  % Focus just on those starting by the 7th group
%        check_comb=40321:factorial(8)/8:80640;    % Focus just on those starting by the 8th group



%        for i=1:length(check_comb)
        for i = 1:num_comb
            v=zeros(1,9);
            test_comb = per_comb(i,:); % It gives me the vaccination order we will follow
            current_group = test_comb(1);
            real_order=zeros(1,9);
            real_order(1)=test_comb(1);
            v(test_comb(1)) = rate_i*sum(N);

            options = odeset('MaxStep',0.1);
            [t,yd] = ode15s(@(t,y)All_combinations_SIYRD_dydn(t,y,M,N),[0,tmax],yd0,options);
            t_Vacc{i,j-1}=t;
            yd_Vacc{i,j-1}=yd;

            S_Vacc = 0;
            I_Vacc = 0;
            Y_Vacc = 0;
            R_Vacc = 0;
            D_Vacc = 0;
            CI_Vacc= 0;
            V_Vacc = 0; % Vaccinated individuals
            YLL_Vacc_num=0;
            Hosp_Vacc_num=0;
            UCI_Vacc_num=0;

            for t=1:9

                S_Vacc = S_Vacc + yd(:,-6+t*7);
                I_Vacc = I_Vacc + yd(:,-5+t*7);
                Y_Vacc = Y_Vacc + yd(:,-4+t*7);
                R_Vacc = R_Vacc + yd(:,-3+t*7);
                D_Vacc = D_Vacc + yd(:,-2+t*7);
                V_Vacc = V_Vacc + yd(:,-1+t*7);
                CI_Vacc= CI_Vacc + yd(:,t*7); %Cumulative new infections

                % We need to divide the final results by 100 because these
                % were percentages.
                Hosp_Vacc_num = Hosp_Vacc_num + max(yd(:,t*7))*Hosp_COVID_Spain_9G(t);       % Total number of patients requiring hospitalization
                UCI_Vacc_num = UCI_Vacc_num + max(yd(:,t*7))*UCI_COVID_Spain_9G(t);          % Total number of patients requiring UCI-hospitalization

                YLL_Vacc_num= YLL_Vacc_num + (SLE-Av_age(t))*max(yd(:,-2+t*7));

            end
            Death_tot(i,j-1)=max(D_Vacc);
            Inf_tot(i,j-1) = max(CI_Vacc);
            Total_vacc(i,j-1) = max(V_Vacc);
            yd_Vacc_byGroups{i,j-1}=[S_Vacc I_Vacc Y_Vacc R_Vacc D_Vacc CI_Vacc];
            vacc_order{i,j-1}=real_order;
            YLL_Vacc(i,j-1) = YLL_Vacc_num;
            Hosp_Vacc(i,j-1) = Hosp_Vacc_num;
            UCI_Vacc(i,j-1) = UCI_Vacc_num;

        end


      end
end

RD_vacc=100-(Death_tot*100/Deaths_noVacc);
RI_vacc=100-(Inf_tot*100/Inf_tot_noVacc);
RYLL_vacc=100-(YLL_Vacc*100/YLL_noVacc);
RH_vacc=100-(Hosp_Vacc*100/Hosp_noVacc);
RUCI_vacc=100-(UCI_Vacc*100/UCI_noVacc);

% Save results
save('v_1_results.mat','Death_tot','Inf_tot','Total_vacc', 'vacc_order','YLL_Vacc','Hosp_Vacc','UCI_Vacc')

toc
