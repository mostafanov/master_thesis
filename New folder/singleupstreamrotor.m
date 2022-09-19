
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
 Ct_up=total_delta_Ct_up
 CT_up= Ct_up
 delta_CQ_i_up=delta_Ct_up.*Lamda_up
 CQ_i_up = sum (delta_CQ_i_up)
 total_CQ_o_up= sum (CQ_o_up)
 cQ_up=(CQ_i_up*(K)*K_int_upstream)+(total_CQ_o_up*Nb)

 CT_up_old = CT_up
               %%%Single Rotor tip vortex radial displacement%%%

K_1_tip_s = 0.25*((CT_up/(sigma*2))+(0.001*seta_twist))
Zeta_passage_s = 2*pi/Nb
Z_Tip_passage_s = K_1_tip_s*Zeta_passage_s
K_2_tip_s = (1+(0.01*seta_twist)) * sqrt(Ct_up)
A_tip_s = 0.145+(27*CT_up)
if Z_Tip_passage_s > H
    Zeta_age_s= (H/K_1_tip_s)*180/pi   
else 
    Zeta_age_s= (((H-Z_Tip_passage_s)/K_2_tip_s)+(Zeta_passage_s))*180/pi    
end