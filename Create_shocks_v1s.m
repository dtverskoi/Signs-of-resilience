X0_ind=NaN(1,T_sh);                                % one shock
a_ind=NaN(1,T_sh);  
b_ind=NaN(1,T_sh);

for t1=1:T_onset                                        % constructing one shock: onset up
    X0_ind(t1)=X0_0+(X0_max-X0_0)*t1/T_onset;
    a_ind(t1)=a_0+(a_max-a_0)*t1/T_onset;
    b_ind(t1)=b_0+(b_max-b_0)*t1/T_onset;
end

for t1=T_onset+1:T_sh-T_onset                            % constructing one shock: peak
    X0_ind(t1)=X0_max;
    a_ind(t1)=a_max;
    b_ind(t1)=b_max;
end

for t1=T_sh-T_onset+1:T_sh
    X0_ind(t1)=X0_max-(X0_max-X0_0)*(t1-T_sh+T_onset)/T_onset;
    a_ind(t1)=a_max-(a_max-a_0)*(t1-T_sh+T_onset)/T_onset;
    b_ind(t1)=b_max-(b_max-b_0)*(t1-T_sh+T_onset)/T_onset;
end

t1=100;
while t1<=T
    if rand<p_sh
        X0(t1:t1+T_sh-1)=X0_ind;
        a(t1:t1+T_sh-1)=a_ind;
        b(t1:t1+T_sh-1)=b_ind;
        t1=t1+T_sh+1;
    else
        t1=t1+1;
    end
end