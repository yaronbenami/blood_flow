clear all
clc

N_seg = 6;  %number of segments
N_node = 6; %number of nodes

a=100;   %segment lengths [microns]
d=20;   %segment diameter [microns]


mu = 1.85e-6; %reference viscosity  [g/micron/s]
gamma=0.5/300;  %5.8660e-04;  %0.0001;  %pressure gradient   [g/micron/s^2]

dt = 2*a*mu/(gamma*d^2);
t=0:dt:30000*dt;

H_0 = 0.45;
alpha = 0.5;
L_14=a; L_25=a; L_45=a; L_46=a; L_56=a; L_63=a;
D_14=alpha*d; D_25=d; D_45=d; D_46=alpha*d; D_56=d; D_63=d;

seg_length = [L_14,L_25,L_45,L_46,L_56,L_63];
seg_diameter = [D_14,D_25,D_45,D_46,D_56,D_63];

% array of node's indices connected to a specific node
conections_node = zeros (3,N_node);
conections_node(:,1) = [4;0;0]; conections_node(:,2) = [5;0;0]; conections_node(:,3) = [6;0;0]; conections_node(:,4) = [1;5;6];
conections_node(:,5) = [2;4;6]; conections_node(:,6) = [3;4;5];

% array of segment's indices connected to a specific node
conections_seg = zeros (3,N_node);
conections_seg(:,1) = [1;0;0]; conections_seg(:,2) = [2;0;0]; conections_seg(:,3) = [6;0;0]; conections_seg(:,4) = [1;3;4];
conections_seg(:,5) = [2;3;5]; conections_seg(:,6) = [6;4;5];

conections_seg_direction = zeros (3,N_node);  % the direction in the conected segment is either 1 for outgoing or -1 for ingoing flow

% array of node's indices connected to a specific segment
seg_nodes = zeros (2,N_seg);
seg_nodes(:,1) = [1;4]; seg_nodes(:,2) = [2;5]; seg_nodes(:,3) = [4;5]; seg_nodes(:,4) = [4;6]; seg_nodes(:,5) = [5;6]; seg_nodes(:,6) = [3;6];

% parameters for calculating viscosity
seg_f = 1 + 10*(seg_diameter/10).^12;
seg_C = (0.8+exp(-0.075*seg_diameter)).*(-1+1./seg_f) + 1./seg_f;
seg_mu45 = 6*exp(-0.085*seg_diameter) + 3.2 -2.44*exp(-0.06*seg_diameter.^0.645);
seg_beta = (seg_diameter./(seg_diameter-1.1)).^2;


seg_perm = zeros(1,N_seg);

seg_ind_t = ones(1,N_seg);
seg_ind_min_change = ones(1,N_seg);
seg_ind_t_tot = 2*ones(1,N_seg);
ini_flag = zeros(1,N_seg);
flag_tot = zeros(1,N_seg);

H_0_t = zeros(length(t),N_seg);
H_L_t = zeros(length(t),N_seg);
H_avg_t = zeros(length(t),N_seg);
mu_avg_t = zeros(length(t),N_seg);

P_t = zeros(length(t),N_node);
Q_t = zeros(length(t),N_seg);
U_t = zeros(length(t),N_seg);

% Initial state
epsilon=1e-6;

P_t(1,:)=[0.5 0.5 0 0 0 0];
H_avg_t(1,:)=[H_0 H_0 0 H_0 H_0 H_0];  % initial condition for average hematocrit
H_0_t(1,:)= H_avg_t(1,:);  % initial condition for H0
H_L_t(1,:)= H_avg_t(1,:);  % initial condition for HL

% boundary conditions
P_t(:,1:2)=0.5;
P_t(:,3)=0;
H_0_t(:,1)=H_0;
H_0_t(:,2)=H_0;
H_0_t(200,2)=H_0-epsilon;

% initial condition for average viscosity
for jj=1:N_seg
    
    mu_avg_t(1,jj) = mu*seg_beta(jj)*(1 + seg_beta(jj)*(seg_mu45(jj)-1)/((1-0.45)^seg_C(jj)-1)*(-1 + (1-H_avg_t(1,jj))^seg_C(jj)));
    
end

P_node = P_t(1,:);

%Pries 1989 with smoothing function
% Q_0 = 0.4/(alpha*d);
% p = 1+6.98*(1-H_0)/(alpha*d);
% r = -6.96/(alpha*d)*log(1/alpha);
% 
% x_0 = Q_0 + 1e-6;
% x_sol_0 = fsolve(@(x)fun_smooth_spline(x,H_0,alpha,d,0),x_0);
% a_3_0 = exp(r)*(x_sol_0-Q_0)^p/(exp(r)*(x_sol_0-Q_0)^p + (1-x_sol_0-Q_0)^p)/x_sol_0^3;
% 
% x_1 = 1 - Q_0 - 1e-6;
% x_sol_1 = fsolve(@(x)fun_smooth_spline(x,H_0,alpha,d,1),x_1);
% a_3_1 = (exp(r)*(x_sol_1-Q_0)^p/(exp(r)*(x_sol_1-Q_0)^p + (1-x_sol_1-Q_0)^p)-1)/(x_sol_1-1)^3;


for ii=1:length(t)  % begin time loop
    
    
    seg_perm = seg_diameter.^4./(16*seg_length.*mu_avg_t(ii,:));  % segment permeabilities
    
    % calculate pressure nodes using Gauss-Seidel iterative method
    eps=1;
    while (eps > 1e-20)
        
        P_node_p = P_node; %saving previous pressure iteration
        
        for jj=4:6  % running over bulk nodes
            
            P_node(jj)= (P_node(conections_node(1,jj))*seg_perm(conections_seg(1,jj)) + P_node(conections_node(2,jj))*seg_perm(conections_seg(2,jj)) + P_node(conections_node(3,jj))*seg_perm(conections_seg(3,jj)))...
                /(seg_perm(conections_seg(1,jj)) + seg_perm(conections_seg(2,jj)) + seg_perm(conections_seg(3,jj)));
        end
        
        eps=max(abs((P_node(4:6)-P_node_p(4:6))./P_node(4:6)));
        
    end
    
    P_t(ii,:) = P_node;
    
    for jj=1:N_seg
        U_t(ii,jj) = seg_diameter(jj)^2/(32*mu_avg_t(ii,jj)*seg_length(jj))*(P_node(seg_nodes(1,jj))-P_node(seg_nodes(2,jj)));
        % if U>0 the flow is from node with smaller index to larger index
    end
    
    Q_t(ii,:) = pi*U_t(ii,:).*seg_diameter.^2/4;
    
    
    % calcualte the retarded time in each segment
    if (ii>1)
        for jj=1:N_seg
            
            if(any(sign(U_t(seg_ind_min_change(jj)+1:ii,jj))~=sign(U_t(seg_ind_min_change(jj):ii-1,jj))))  %if velocity changed its sign
                
                ind_change_last_two = find(sign(U_t(seg_ind_min_change(jj)+1:ii,jj))~=sign(U_t(seg_ind_min_change(jj):ii-1,jj)),2,'last') + seg_ind_min_change(jj);
                
                if (ind_change_last_two(end)==ii)
                    
                    seg_ind_t(jj) = ind_change_last_two(end)-2;
                    
                    if(length(ind_change_last_two)>1)    % if there were at least two changes in velocity sign
                        seg_ind_min_change(jj) = ind_change_last_two(1);
                    end
                    
                end
                ind_change = ind_change_last_two(end);
                
                if (ind_change<ii)  %otherwise can't calculate integral
                    L_change = abs(trapz(t(ind_change:ii),U_t(ind_change:ii,jj))); %propgation length since last sign change
                    
                    if (L_change<seg_length(jj))   % if the hematocrit comes for the same end of the segment
                        
                        I_max = abs(trapz(t(seg_ind_min_change(jj):ind_change),U_t(seg_ind_min_change(jj):ind_change,jj)));
                        
                        if (I_max >= L_change)
                            
                            for kk = seg_ind_t(jj):-1:seg_ind_min_change(jj)
                                I = abs(trapz(t(kk:ind_change),U_t(kk:ind_change,jj)));
                                if (I >= L_change)
                                    tau_ind = kk;
                                    seg_ind_t(jj)=kk;
                                    break;
                                end
                            end
                            if (U_t(ii,jj)>0)
                                H_L_t(ii,jj)=H_L_t(tau_ind,jj);
                            else
                                H_0_t(ii,jj)=H_0_t(tau_ind,jj);
                            end
                            
                        end
                        
                        I_tot = trapz(t(seg_ind_t_tot(jj)-1:ii),U_t(seg_ind_t_tot(jj)-1:ii,jj));
                        
                        if(abs(I_tot) > seg_length(jj))   % this is the case where velocity change sign but there is an accumulative effect
                            
                            H_max=0;
                            
                            for kk=seg_ind_t_tot(jj)-1:ii-1
                                I = abs(trapz(t(kk:ii),U_t(kk:ii,jj)));
                                if (I <= seg_length(jj))
                                    seg_ind_t_tot(jj) = kk;
                                    if (I_tot>0)
                                        if (H_0_t(kk,jj) > H_max)
                                            H_max = H_0_t(kk,jj);
                                        end
                                    else
                                        if (H_L_t(kk,jj) > H_max)
                                            H_max = H_L_t(kk,jj);
                                        end
                                    end
                                    break;
                                else
                                    if (I_tot>0)
                                        if (H_0_t(kk,jj) > H_max)
                                            H_max = H_0_t(kk,jj);
                                        end
                                    else
                                        if (H_L_t(kk,jj) > H_max)
                                            H_max = H_L_t(kk,jj);
                                        end
                                    end
                                end
                            end
                            
                            if (I_tot>0)
                                H_L_t(ii,jj) = H_max ;
                            else
                                H_0_t(ii,jj) = H_max ;
                            end
                            
                            
                        else
                            if (I_tot>0)
                                seg_ind_t_tot(jj) =  seg_ind_t_tot(jj)-1 + find(U_t(seg_ind_t_tot(jj):ii,jj)>0,1);
                            else
                                seg_ind_t_tot(jj) =  seg_ind_t_tot(jj)-1 + find(U_t(seg_ind_t_tot(jj):ii,jj)<0,1);
                            end
                        end
                        
                        
                        
                    else  % if the hematocrit comes for the other end of the segment
                        seg_ind_min_change(jj) = ind_change;
                        seg_ind_t(jj) = seg_ind_min_change(jj);
                        tau_ind = seg_ind_min_change(jj);
                        if (U_t(ii,jj)>0)
                            H_L_t(ii,jj)=H_0_t(tau_ind,jj);
                        else
                            H_0_t(ii,jj)=H_L_t(tau_ind,jj);
                        end
                        
                    end
                    
                end
                
            else  %if the velocity at the same direction since last time the propagation length reached the segment length
                
                tau_ind=1;
                
                for kk=seg_ind_t(jj):ii
                    I = abs(trapz(t(kk:ii),U_t(kk:ii,jj)));
                    if (I <= seg_length(jj))    % I either changed from > to < or if it is straight < then tau=t(1) (initial condition)
                        if (ini_flag(jj)==1)
                            tau_ind = kk;
                            seg_ind_t(jj)=kk;
                            break;
                        else
                            tau_ind = 1;
                            seg_ind_t(jj)=1;
                            break;
                        end
                    else
                        ini_flag(jj)=1;
                    end
                end
                
                if (U_t(ii,jj)>0)
                    H_L_t(ii,jj)=H_0_t(tau_ind,jj);
                else
                    H_0_t(ii,jj)=H_L_t(tau_ind,jj);
                end
                
            end
            
            
        end
    end
    
    % calculating hematocrit at the beginning of each segment
    for jj=4:6  % running over bulk nodes
        
        % conections_seg_direction is 1 for outgoing flow and -1 for ingoing flow
        conections_seg_direction(1,jj) = (P_node(jj) - P_node(conections_node(1,jj)))/abs(P_node(jj) - P_node(conections_node(1,jj)));
        if(P_node(jj) - P_node(conections_node(1,jj))==0)
            conections_seg_direction(1,jj) = 0;
        end
        conections_seg_direction(2,jj) = (P_node(jj) - P_node(conections_node(2,jj)))/abs(P_node(jj) - P_node(conections_node(2,jj)));
        if(P_node(jj) - P_node(conections_node(2,jj))==0)
            conections_seg_direction(2,jj) = 0;
        end
        conections_seg_direction(3,jj) = (P_node(jj) - P_node(conections_node(3,jj)))/abs(P_node(jj) - P_node(conections_node(3,jj)));
        if(P_node(jj) - P_node(conections_node(3,jj))==0)
            conections_seg_direction(3,jj) = 0;
        end
        
        ind_conv = conections_seg(conections_seg_direction(:,jj)<1,jj);  % the segments flowing into the junction
        ind_div = conections_seg(conections_seg_direction(:,jj)==1,jj);    % the segments flowing out of the junction
        
        Q_conv = Q_t(ii,ind_conv);
        Q_div = Q_t(ii,ind_div);
            
        if (conections_seg_direction(1,jj)+conections_seg_direction(2,jj)+conections_seg_direction(3,jj) < 1)
            % two converging segments  -- mass balance
            
            % calculating H at the end of the converging segments
            H_conv=zeros(1,2);
            for kk=1:2
                if (Q_conv(kk)>0)
                    H_conv(kk) = H_L_t(ii,ind_conv(kk));
                else
                    H_conv(kk) = H_0_t(ii,ind_conv(kk));
                end
            end
            if (Q_div>0)
                H_0_t(ii,ind_div) = (H_conv(1)*abs(Q_conv(1)) + H_conv(2)*abs(Q_conv(2)))/abs(Q_div);
            else
                H_L_t(ii,ind_div) = (H_conv(1)*abs(Q_conv(1)) + H_conv(2)*abs(Q_conv(2)))/abs(Q_div);
            end
            
        else    % two diverging segments -- haematocrit splitting rule
            
            % calculating H at the end of the converging segment
            if (Q_conv>0)
                H_conv = H_L_t(ii,ind_conv);
            else
                H_conv = H_0_t(ii,ind_conv);
            end
            
            psi_1 = abs(Q_div(1)/Q_conv);
            psi_2 = 1 - psi_1;
            
            %Pries (Gardner et. al's (2010) paper)
            Q_0=0.4/seg_diameter(ind_conv);

            if (psi_1 >= Q_0 && psi_1 <= 1-Q_0)
                p=1 + 6.98*(1-H_conv)/seg_diameter(ind_conv);
                exp_r = (seg_diameter(ind_div(1))/seg_diameter(ind_div(2)))^(-6.96/seg_diameter(ind_conv));
                phi_1 = exp_r*(psi_1-Q_0)^p/(exp_r*(psi_1-Q_0)^p + (1-psi_1-Q_0)^p);
            else
                if(psi_1 < Q_0)
                    phi_1=0;
                else
                    phi_1=1;
                end
            end

            %Pries with smoothing
%             Q_0=0.4/seg_diameter(ind_conv);
% 
%             if (psi_1 >= x_sol_0 && psi_1 <= x_sol_1)
%                 p=1 + 6.98*(1-H_conv)/seg_diameter(ind_conv);
%                 exp_r = (seg_diameter(ind_div(1))/seg_diameter(ind_div(2)))^(-6.96/seg_diameter(ind_conv));
%                 phi_1 = exp_r*(psi_1-Q_0)^p/(exp_r*(psi_1-Q_0)^p + (1-psi_1-Q_0)^p);
%             else
%                 if(psi_1 < x_sol_0)
%                     phi_1 = a_3_0*psi_1^3;
%                 else
%                     phi_1 = 1 + a_3_1*(psi_1-1)^3;
%                 end
%             end

            %Klitzman & Johnson (1982)
%             p = 2.001;
%             phi_1 = psi_1^p/(psi_1^p + (1-psi_1)^p);
            
            phi_2 = 1 - phi_1;
            
            if (Q_div(1)>0)
                H_0_t(ii,ind_div(1)) = phi_1/psi_1*H_conv;
            else
                H_L_t(ii,ind_div(1)) = phi_1/psi_1*H_conv;
            end
            if (Q_div(2)>0)
                H_0_t(ii,ind_div(2)) = phi_2/psi_2*H_conv;
            else
                H_L_t(ii,ind_div(2)) = phi_2/psi_2*H_conv;
            end
            
        end
        
    end
    
    % advancing hematocrit to next time step
    if (ii<length(t))
        for jj=1:N_seg
            H_avg_t(ii+1,jj) = H_avg_t(ii,jj) + dt*U_t(ii,jj)/seg_length(jj)*(H_0_t(ii,jj) - H_L_t(ii,jj));
            D = mu*seg_beta(jj)*seg_beta(jj)*(seg_mu45(jj)-1)/((1-0.45)^seg_C(jj)-1);
            mu_avg_t(ii+1,jj) = mu_avg_t(ii,jj) + dt*D*U_t(ii,jj)/seg_length(jj)*((1-H_0_t(ii,jj))^seg_C(jj) - (1-H_L_t(ii,jj))^seg_C(jj));
        end
    end
    
    ii   % show time step on screen
    
end  % end time loop

Q_m = mean(Q_t(:,2));
U_m = mean(U_t(:,2));
t_norm = t*U_m/a;
figure(2)
hold on
plot(t_norm,H_avg_t,'k-',t_norm,repmat(H_avg_t(1,:),length(t),1),'k--')
figure(1)
hold on
plot(t_norm,Q_t(:,3)./Q_m,'k-')





