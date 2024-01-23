function [u] = utility_v1s(ch_y,R,c_y,A1_aux,A2_aux,A3_aux,i,gamma,conf,ch_x,ch_z,X_bar,Z_bar,A4_aux,s1,N,nonmat_z,nonmat_x)

pi_aux=A1_aux*((1-s1)*R(i)*(ch_y^gamma)+s1*(A4_aux+R(i)*(ch_y^gamma))/N)+A2_aux+A3_aux-c_y*ch_y;
u=pi_aux-conf*(ch_x-X_bar)^2-conf*(ch_z-Z_bar)^2+nonmat_z*ch_z+nonmat_x*ch_x;

end
