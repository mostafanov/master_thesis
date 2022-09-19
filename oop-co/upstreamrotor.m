function [Lamda_ind_up_s, v_ind_up_s_0, v_center_s, Lamda_Center_s, v_ind_up_s, CT_up, y_s, Lamda_D_old, v_center_old, V_D_old, v_ind_up_old, Lamda_2_s_old, Lamda_2_s, v_center_2_co_old, cQ_up] = upstreamrotor(sigma, seta_o_up, seta_twist, omga, Cl_alpha, R_cut, b, delta_r, taper_slope, K_int_upstream, Nb, roh, Mioh, N, NN, K, k_br, H, C_root, taper_ratio, R, R_cut_factor)
    function [sigma, seta_o_up, seta_twist, omga, Cl_alpha, R_cut, b, delta_r, taper_slope, K_int_upstream, Nb, roh, Mioh, N, NN, K, k_br, H, C_root, taper_ratio, R, R_cut_facto] = rotorparameters(sigma, seta_o_up, seta_twist, omga, Cl_alpha, R_cut, b, delta_r, taper_slope, K_int_upstream, Nb, roh, Mioh, N, NN, K, k_br, H, C_root, taper_ratio, R, R_cut_factor)
    end

        %-------------------------Inital guess upstream rotor----------------------%
 for i=1:N
     r_up(i)=((delta_r*i))/b;                                              % radius from center of each element to the end of span of a single blade
     seta_up(i)= seta_o_up+ r_up(i)*seta_twist;                            % pitch angle for each element
     c_up(i)=(C_root)+(taper_slope*r_up(i)*b);                             % chord of each element
     Re_up_T(i) = (roh*(omga* (((r_up(i))*b)+R_cut))*c_up(i))/Mioh;        % Reynolds number at each element

%---------------------------Lift slope for upstream rotor------------------%
  %%calculating lift slope for each element%%
    if seta_up(i) <= 2*(pi/180)
         if Re_up_T(i) < 70000
             Cl_alpha_up(i)= ((0.2325/70000)* Re_up_T(i))*180/pi
         else Cl_alpha_up(i)=((0.21+(20000*0.7/Re_up_T(i)))/2)*180/pi
         end

     elseif  2*(pi/180)< seta_up(i) <= 4*(pi/180)
         if Re_up_T(i) < 70000
             Cl_alpha_up(i)= ((0.152/70000)* Re_up_T(i))*180/pi
         elseif 70000 <= Re_up_T(i) <= 150000 
             Cl_alpha_up(i)=0.152*180/pi
         else Cl_alpha_up(i)=((0.21+(20000*0.7/Re_up_T(i)))/2)*180/pi
         end
         
     else Cl_alpha_up(i)= Cl_alpha

     end
     
 %-------------------upstream rotor performance parameter calculation------

     Lamda_up(i)=((Cl_alpha_up(i)*sigma)/16)*((sqrt (1+((32*seta_up(i)*r_up(i))/(Cl_alpha_up(i)*sigma))))-1);
     B_up(i) = 1-((Lamda_up(i)*sqrt(2)))
     delta_Ct_up(i)=((Cl_alpha_up(i)*sigma)/2)*((seta_up(i)*(r_up(i)*B_up(i))^2)-(Lamda_up(i)*r_up(i)/B_up(i)))*(delta_r/b) ;
     Cl_up(i)=Cl_alpha_up(i)*(seta_up(i)-(Lamda_up(i)/(r_up(i)*B_up(i))));
     V_up(i)= (omga* (((r_up(i))*b)+R_cut))
     Re_up(i)= (roh*V_up(i)*c_up(i))/Mioh
     Cd_0_up(i)=2/((Re_up(i)*0.75/1000)+35)
     
     cdup_1_o (i) =(Lamda_up(i)/(r_up(i)))*pi/180
     cdup_2_o (i) =(seta_up(i))
     Cd_o_up(i) = (Cd_0_up(i)+(0.02*cdup_2_o (i))+(0.4*(cdup_2_o (i))^2));
     CQ_o_up(i)=(sigma/8)*Cd_o_up(i)
     seta_o_up = seta_up(i);
 end
 total_delta_Ct_up =(sum (delta_Ct_up))
 Ct_up=total_delta_Ct_up*Nb
 CT_up= Ct_up
 delta_CQ_i_up=delta_Ct_up.*Lamda_up
 CQ_i_up = sum (delta_CQ_i_up)
 total_CQ_o_up= sum (CQ_o_up)*(delta_r/b)
 B_loss_up= 1-(sqrt(Ct_up/2))
 cQ_up=(CQ_i_up*(K)*K_int_upstream*Nb/B_loss_up)+(total_CQ_o_up*Nb)
 
               %%%Single Rotor tip vortex radial displacement%%%

K_1_tip_s = 0.25*((Ct_up/(sigma*Nb))+(0.001*seta_twist))
Zeta_passage_s = 2*pi/Nb
Z_Tip_passage_s = K_1_tip_s*Zeta_passage_s
K_2_tip_s = (1+(0.01*seta_twist)) * sqrt(Ct_up)
A_tip_s = 0.145+(27*Ct_up)
if Z_Tip_passage_s > H
    Zeta_age_s= (H/K_1_tip_s)*180/pi   
else 
    Zeta_age_s= (((H-Z_Tip_passage_s)/K_2_tip_s)+(Zeta_passage_s))*180/pi    
end

%Venna Region info region when radial displacement almost equal sqrt(2)

Zeta_age_vena_s = (180/pi)*-log((0.781-0.78)/(1-0.78))/(A_tip_s)
Y_tip_venna_s= 0.78+(1-0.78)*exp(-A_tip_s*(Zeta_age_vena_s*pi/180))
z_venna_s = Z_Tip_passage_s + K_2_tip_s*((Zeta_age_vena_s*pi/180)-Zeta_passage_s)

%Ytip for Single Rotot (convergance divargance where the throat at venna region) 

if Zeta_age_s< Zeta_age_vena_s
    Y_tip_s= 0.78+(1-0.78)*exp(-A_tip_s*(Zeta_age_s*pi/180))
else
    Y_tip_c_s= Y_tip_venna_s +((0.2^(1/A_tip_s))*exp(A_tip_s*(Zeta_age_s*pi/180)))
    if Y_tip_c_s<1
        Y_tip_s=Y_tip_c_s
    else
        Y_tip_s=1
    end
end
         %%% Velocity Profile for single rotor at plane downstream = H %%%
Lamda_ind_up_s=Lamda_up(N)
v_ind_up_s_0=Lamda_ind_up_s*omga*R
v_center_s= 0.5*v_ind_up_s_0*((1-(Y_tip_s^2))/(Y_tip_s^2))
Lamda_Center_s= v_center_s/(omga*R)
v_ind_up_s=v_ind_up_s_0+v_center_s

% Intial Value for Co-axial Calculation

CT_up_old = CT_up
x_s  =(b-(1-Y_tip_s)*R)/b
X_s = round(x_s*N)
y_s =X_s/N
Lamda_D_old = Lamda_up(N*y_s)
v_center_old= v_center_s
V_D_old= Lamda_D_old*(omga*R)
v_ind_up_old=v_ind_up_s_0
Lamda_1_s =Nb*sum(Lamda_up)*delta_r/b
Lamda_2_s= Lamda_1_s/(Y_tip_s^2)
Lamda_2_s_old=Lamda_2_s
v_center_2_co_old=0
end
