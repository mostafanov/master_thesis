%input
sigma=0.042/2;                                                             % solidaty ratio0.042 0.027
seta_o_up=12*(pi/180)
H=0.4
delta_seta =6*(pi/180)
omga= 130
C_root = 0.062                                                             % 0.062 0.294
taper_ratio = 0.637                                                       % 0.637 0.388 
R= 0.762 %3.81 % ;                                                            %rotor redius
Nb =2
R_cut_factor = 0.167                                                       % 0.167 0.13
seta_twist=0*(pi/180);                                                     %twist angle
roh = 1.225
Mioh = 0.0000181206

                        
                                      %%%%parameter%%%%%

Cl_alpha=2*pi;                                                             % airfoil coefficient of lift slope
seta_o_down=seta_o_up+delta_seta;                                          %pitch angle downstream rotor
R_cut=R_cut_factor*R;                                                      %rotor cut off
N =400                                                                     %descritization elements
NN=20
b = R-R_cut                                                                %blade span
delta_r =(b/N)                               
taper_slope = ((C_root*taper_ratio)-C_root)/b

B_loss=(1-((C_root*(1+(0.7*taper_ratio)))/(1.5*R)))

K =1.15


if H < 0.35
K_int_upstream =1
else
    K_int_upstream =1
end
k_br =0.5


                                   %%%%%%%upstream rotor %%%%%%%
 for i=1:N
     r_up(i)=((delta_r*i))/b;
     seta_up(i)= seta_o_up+ r_up(i)*seta_twist;
     c_up(i)=(C_root)+(taper_slope*r_up(i)*b)
     Re_up_T(i) = (roh*(omga* (((r_up(i))*b)+R_cut))*c_up(i))/Mioh

                                        %%Lift slope%%
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

                    %%upstream rotor performance parameter

     Lamda_up(i)=((Cl_alpha_up(i)*sigma)/16)*((sqrt (1+((32*seta_up(i)*r_up(i))/(Cl_alpha_up(i)*sigma))))-1);
     B_up(i) = 1-((Lamda_up(i)*sqrt(2)))
     delta_Ct_up(i)=((Cl_alpha_up(i)*sigma)/2)*((seta_up(i)*(r_up(i)*B_up(i))^2)-(Lamda_up(i)*r_up(i)/B_up(i)))*(delta_r/b) ;
     Cl_up(i)=Cl_alpha_up(i)*(seta_up(i)-(Lamda_up(i)/r_up(i)));
     V_up(i)= sqrt (((omga* (((r_up(i))*b)+R_cut))^2)+((Lamda_up(i)*omga* (((r_up(i))*b)+R_cut))^2))
     Re_up(i)= (roh*V_up(i)*c_up(i))/Mioh
     Cd_0_up(i)=2/((Re_up(i)*0.75/1000)+35)
     
     cdup_1_o (i) =(Lamda_up(i)/(r_up(i)))*pi/180
     cdup_2_o (i) =(seta_up(i)-cdup_1_o (i))
     Cd_o_up(i) = (Cd_0_up(i)+(0.02*cdup_2_o (i))+(0.4*(cdup_2_o (i))^2))*(delta_r/b);
     CQ_o_up(i)=(sigma/8)*Cd_o_up(i)
     seta_o_up = seta_up(i);
 end
 total_delta_Ct_up =(sum (delta_Ct_up))
 Ct_up=total_delta_Ct_up*Nb
 CT_up= Ct_up*B_loss
 delta_CQ_i_up=delta_Ct_up.*Lamda_up
 CQ_i_up = sum (delta_CQ_i_up)
 total_CQ_o_up= sum (CQ_o_up)
 cQ_up=(CQ_i_up*(K)*K_int_upstream)+(total_CQ_o_up*Nb)
 CT_up_old = CT_up
               %%%Single Rotor tip vortex radial displacement%%%

K_1_tip_s = 0.25*((CT_up/(sigma*Nb))+(0.001*seta_twist))
Zeta_passage_s = 2*pi/Nb
Z_Tip_passage_s = K_1_tip_s*Zeta_passage_s
K_2_tip_s = (1+(0.01*seta_twist)) * sqrt(Ct_up)
A_tip_s = 0.145+(27*CT_up)
if Z_Tip_passage_s > H
    Zeta_age_s= (H/K_1_tip_s)*180/pi   
else 
    Zeta_age_s= (((H-Z_Tip_passage_s)/K_2_tip_s)+(Zeta_passage_s))*180/pi    
end

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
v_ind_up_s=sqrt(Ct_up/2)*omga*R
Tip_Thrust_r_s= sqrt(CT_up/2)/(Lamda_up(N))
v_center_s= v_ind_up_s*((1-(Y_tip_s^3))/(Y_tip_s^3))
Lamda_Center_s= v_center_s/(omga*R)

