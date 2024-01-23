function [sh,T] = Experiment(par0,par1,par2,par3,par4,exp)
% par0=N=8,16,24;
% par1=p_sh=0.05;
% par2=sev=1,2,3,4,5,6;
% par3=T_onset=1,2,3,4,5;
% par4=conf=0:0.1:2;
% Overall, 408 jobs

sh=0;
Runs=100;                         % number of runs
T=1000;                           % time steps
Output=NaN(4*Runs,T);            % Output

for run=1:Runs
    sh=sh+1
    seedStat=0;
    %Seed
    if seedStat==0
        rng('shuffle','twister')
        seed=rng;
        %  save seed;
    else
        %  load seed;
        if exist('seed')
            rng(seed);
        else
            rng('shuffle','twister');
            seed=rng;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Main parameters
    r=0.35;                          % synergy effect
    s1=0.25;                         % sharing rate 1
    s2=1;                            % sharing rate 2
    
    p_sh=par1;                       % probability that shock occurs
    T_sh=10;                         % How long the schock lasts
    T_onset=par3;                    % time for onset
    sev=par2;                        % severity: 6 - max, 0 - min
    
    conf=par4;                       % the strength of conformity in collective actions
    nonmat_z=0.5;                    % non-material payoff associated with z=1
    nonmat_x=0.25;                   % non-material payoff associated with x=1
    
    % Other parameters
    N=par0;                          % number of individuals
    c_x=1;                           % labor costs
    c_y=4;                           % labor costs
    c_z=1;                           % labor costs
    gamma=0.5;                       % Cobb_Douglas production function parameter
    R=16*ones(N,1);
    B=56;
    Z0=0.5*N;                        % the corresponding half effort
    mu1=0.25;                        % updating probability
    
    % Environmental effects
    X0=0.25*N*ones(1,T);             % half effort
    a=zeros(1,T);                    % unmitigated shock on B
    b=0.4*ones(1,T);                 % mitigated shock on R
    
    X0_0=0.25*N;                            % initial HE in CA to mitigate schock
    X0_max=0.25*N+0.0625*N*sev;
    a_0=0;                                  % unmitig. regular costs in regular CA
    b_0=0.4;                                % mitig. regular costs in ind. resource
    a_max=0.0+0.125*sev;
    b_max=0.4+0.1*sev;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create shocks:
    Create_shocks_v1s
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model variables
    x=NaN(N,T);                       % contribution to communal work
    y=NaN(N,T);                       % labor effort
    z=NaN(N,T);                       % contribution to a regular CA
    pi=NaN(N,T);                      % material payoff
    
    % Initial conditions
    y0=1;                             % mean initial labor effort
    sigma_y=0.2;                      % std in initial values of y
    y(:,1)=y0+sigma_y*randn(N,1);     % initial values of y
    y(y(:,1)<0,1)=0;                  % corrections
    
    px0=0.2;                          % initial probability for x to be 1
    pz0=0.2;                          % initial probability for z to be 1
    for i=1:N
        if rand<px0
            x(i,1)=1;                 % initial values of x
        else
            x(i,1)=0;
        end
        if rand<pz0
            z(i,1)=1;                 % initial values of z
        else
            z(i,1)=0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The dynamics
    for t=2:T
        for i=1:N
            if rand<mu1
                CH_y=cell(2,2);
                Ui=NaN(2,2);
                for ch_x=0:1
                    for ch_z=0:1
                        X_aux=sum(x(:,t-1))-x(i,t-1);
                        X_bar=sum(x(:,t-1))/N;
                        Z_aux=sum(z(:,t-1))-z(i,t-1);
                        Z_bar=sum(z(:,t-1))/N;
                        A1_aux=1-b(t)*(1-(X_aux+ch_x)/(X_aux+ch_x+X0(t)/((1+(Z_aux+ch_z)/N)^r)));
                        if Z_aux+ch_z~=0
                            A2_aux=B*(1-a(t))*((Z_aux+ch_z)/(Z_aux+ch_z+Z0/((1+(X_aux+ch_x)/N)^r)))*((1-s2)*ch_z/(Z_aux+ch_z)+s2/N);
                        else
                            A2_aux=0;
                        end
                        A3_aux=-c_x*ch_x-c_z*ch_z;
                        A4_aux=sum(R.*(y(:,t-1).^gamma))-R(i)*(y(i,t-1)^gamma);
                        fun=@(ch_y)-utility_v1s(ch_y,R,c_y,A1_aux,A2_aux,A3_aux,i,gamma,conf,ch_x,ch_z,X_bar,Z_bar,A4_aux,s1,N,nonmat_z,nonmat_x);
                        options = optimoptions('fmincon','Display','off');
                        [CH_y{ch_x+1,ch_z+1},Ui(ch_x+1,ch_z+1)] = fmincon(fun,1,[],[],[],[],0,inf,[],options);
                    end
                end
                [~,max_idx] = max(-Ui,[],'all','linear');                                        % find optimal x and z
                [ch_x_opt,ch_z_opt]=ind2sub([2 2],max_idx);
                x(i,t)=ch_x_opt-1;                                                               % save optimal x
                z(i,t)=ch_z_opt-1;                                                               % save optimal x
                y(i,t)=CH_y{ch_x_opt,ch_z_opt};                                                  % save optimal y                                                                      % save optimal r
            else
                x(i,t)=x(i,t-1);                                                                 % save optimal e
                z(i,t)=z(i,t-1);                                                                 % save optimal s
                y(i,t)=y(i,t-1);                                                                 % save optimal r
            end
        end
        X_tot=sum(x(:,t));
        Z_tot=sum(z(:,t));
        A1=1-b(t)*(1-X_tot/(X_tot+X0(t)/((1+Z_tot/N)^r)));
        A4=sum(R.*(y(:,t).^gamma));
        for i=1:N
            if Z_tot~=0
                A2=B*(1-a(t))*(Z_tot/(Z_tot+Z0/((1+X_tot/N)^r)))*((1-s2)*z(i,t)/Z_tot+s2/N);
            else
                A2=0;
            end
            A3=-c_x*x(i,t)-c_z*z(i,t);
            pi(i,t)=A1*((1-s1)*R(i)*(y(i,t)^gamma)+s1*A4/N)+A2+A3-c_y*y(i,t);
        end
    end
    
    Output(run,:)=sum(x,1);
    Output(Runs+run,:)=mean(y,1);
    Output(2*Runs+run,:)=sum(z,1);
    Output(3*Runs+run,:)=mean(pi,1);
    
end

writematrix(Output,['Experiment_' num2str(exp) num2str(par0) num2str(par1) num2str(par2) num2str(par3) num2str(par4) '.txt'])

end