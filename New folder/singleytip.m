%Venna Region info region when radial displacement almost equal sqrt(2)

Zeta_age_vena_s = (180/pi)*-log((0.71-0.702)/(1-0.702))/(A_tip_s)
Y_tip_venna_s= 0.702+(1-0.702)*exp(-A_tip_s*(Zeta_age_vena_s*pi/180))
z_venna_s = Z_Tip_passage_s + K_2_tip_s*((Zeta_age_vena_s*pi/180)-Zeta_passage_s)

%Ytip for Single Rotot (convergance divargance where the throat at venna region) 

if Zeta_age_s< Zeta_age_vena_s
    Y_tip_s= 0.702+(1-0.702)*exp(-A_tip_s*(Zeta_age_s*pi/180))
else
    Y_tip_c_s= Y_tip_venna_s +((0.2^(1/A_tip_s))*exp(A_tip_s*(Zeta_age_s*pi/180)))
    if Y_tip_c_s<1
        Y_tip_s=Y_tip_c_s
    else
        Y_tip_s=1
    end
end
         %%% Velocity Profile for single rotor at plane downstream = H %%%
v_ind_up_s=(Lamda_up(N))*omga*R
Tip_Thrust_r_s= sqrt(CT_up/2)/(Lamda_up(N))
v_center_s= v_ind_up_s*((1-(Y_tip_s^3))/(Y_tip_s^3))
Lamda_Center_s= v_center_s/(omga*R)