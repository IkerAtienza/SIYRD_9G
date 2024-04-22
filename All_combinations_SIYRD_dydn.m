function dydt = All_combinations_SIYRD_dydn(t,y,M,N)


global b_si b_sy b_ri b_ry v rate_i test_comb
global mu_I mu_Y r condition Vacc_limit current_group real_order

    m=length(M);  %Number of age population groups

    % y(1+7n): S_m      --> Susceptible individuals from group-m
    % y(2+7n): I_m  	--> First-infected individuals from group-m
    % y(3+7n): Y_m      --> Reinfected individuals from group-m
    % y(4+7n): R_m      --> Recovered individuals from group-m
    % y(5+7n): D_m      --> Dead individuals from group-m.
    % y(6+7n): ## of vaccines in subgroup i
    % y(7+7n): ## of new infections in subgroup i

    %[M,N1,N2]=Matrix_conversion;
    % Closure condition: M_ij * N_i = M_ji * N_j --> (M_ij/M_ji)=Nj/Ni.
    
    % Now, we do not know a priori which group will be vaccinated first
    % neither which order vaccination will follow. In each iteration, such
    % order will be given by the elements in the vector test_comb.
    
    c=1;
    c1=0;
    c2=2;
    if condition == 2
        
        if y(-6+7*current_group) <= (1-Vacc_limit) * N(current_group) && v(current_group)~=0 
            v(current_group) = 0;
            
            while c1==0 && c2<=length(test_comb)
                if y(-6+7*test_comb(c2)) > (1-Vacc_limit) * N(test_comb(c2)) && v(test_comb(c2))==0 
                    current_group = test_comb(c2);
                    real_order(c2)=test_comb(c2);
                    v(current_group) = rate_i*sum(N);
                    c1=1;
                end 
                c2=c2+1;
            end
        end
    end
            
            
    
    
%     %% G1-first vaccination:
%     %   -  If one subgroup G1(0-9), G1(10-19), G1(20-29), G1(30-39) and
%     %   G1(40-49) runs out of susceptible individual (y(1+6n)<...), the
%     %   remaining vaccination rates must be readapted accordingly.
%     if condition==2
%         c1=1;
%         c2=0;
%         G1_stop(1)=0;
%         N_ok=0;
%         for i=1:age_limit
%             if y(-6+7*i)<=(1-Vacc_limit)*N(i) && v(i)~=0 && lim1(i)==0
%                G1_stop(c1)=i;
%                c1=c1+1;
%                v(i)=0;
%                lim1(i)=y(-6+7*i);
%             elseif y(-5+6*i)>(1-Vacc_limit)*N(i) && v(i)~=0 && lim1(i)==0
%                 N_ok=N_ok+N(i);
%                 c2=c2+1;
%             end
%         end
%         % First, with this for loop and if-statement I check whether any
%         % subgroup of G1 has already reached its limit.
%         if (G1_stop(1)~=0 && c2>=1) %% Only if there are more subgroups in G1 available to receive more vaccines I go inside this for loop, otherwise, I assign the vaccines
%             % the subgroups within G2. 
%              for i=1:age_limit
%                  if v(i)~=0
%                     v(i)=rate_i*sum(N)*N(i)/N_ok;
%                  end
%              end
%         elseif (G1_stop(1)~=0  && c2==0)
%             %I start vaccination of G2-subgroups  
%              for i=(age_limit+1):9
%                   v(i)=rate_i*sum(N)*N(i)/sum(N((age_limit+1):9));  % First older group, composed in this case by 50-59, 60-69, 70-79 and >80
%              end
%         end
% 
%         c3=1;
%         c4=0;
%         G2_stop(1)=0;
%         for i=age_limit+1:9
%             if y(-6+7*i)<=(1-Vacc_limit)*N(i) && v(i)~=0 && lim2(i)==0
%                G2_stop(c3)=i;
%                c3=c3+1;
%                v(i)=0;
%                lim2(i)=y(-6+7*i);
%             else
%                 N_ok=N_ok+N(i);
%                 c4=c4+1;
%             end
%         end
% 
%         if (G2_stop(1)~=0 && c4>=1)
%              for i=age_limit+1:9
%                  if v(i)~=0
%                     v(i)=rate_i*sum(N)*N(i)/N_ok;
%                  end
%              end
%         end 
%     end
 
 
 
%% Simultaneous vaccination 
% if condition==4
%         c1=1;
%         c2=0;
%         GT_stop(1)=0;
%         N_ok=0;
%         for i=1:9
%             if y(-6+7*i)<=(1-Vacc_limit)*N(i) && v(i)~=0 && lim2(i)==0
%                GT_stop(c1)=i;
%                c1=c1+1;
%                v(i)=0;
%                lim2(i)=y(-6+7*i);
%             elseif y(-5+6*i)>(1-Vacc_limit)*N(i) && v(i)~=0 && lim2(i)==0
%                 N_ok=N_ok+N(i);
%                 c2=c2+1;
%             end
%         end
%         
%         if (GT_stop(1)~=0 && c2>=1)
%              for i=1:9
%                  if v(i)~=0
%                     v(i)=rate_i*sum(N)*N(i)/N_ok;
%                  end
%              end
%         end 
% end


%%    
    count=1;
    for i=1:m
        A=0;
        B=0;
        for j=1:m
            A=A+M(i,j)*(y(-5+7*j)/N(j));       %We are taking just I_m (Parentheses in (56) first part without betas)
            B=B+M(i,j)*(y(-4+7*j)/N(j));     % We are taking just Y_m (Parentheses in (56) second part without betas)
        end
            EQ1 = -(b_si*A + b_sy*B)*y(count)- v(i);

            EQ2 =  (b_si*A + b_sy*B)*y(count) - r(i)*y(count+1) - mu_I(i)*y(count+1);

            EQ3 =  (b_ri*A + b_ry*B)*y(count+3) - r(i)*y(count+2) - mu_Y(i)*y(count+2); 

            EQ4 = -(b_ri*A + b_ry*B)*y(count+3) + r(i)*y(count+1) + r(i)*y(count+2) + v(i); 

            EQ5 =  mu_I(i)*y(count+1) + mu_Y(i)*y(count+2);
            
            EQ6 = v(i);

            EQ7= (b_si*A + b_sy*B)*y(count) +(b_ri*A + b_ry*B)*y(count+3); % Equation with cumulative infections.

        dydt(count,1)=EQ1;
        dydt(count+1,1)=EQ2;
        dydt(count+2,1)=EQ3;
        dydt(count+3,1)=EQ4;
        dydt(count+4,1)=EQ5;
        dydt(count+5,1)=EQ6;
        dydt(count+6,1)=EQ7;
        count=(count+6)+1;
    end
    
end