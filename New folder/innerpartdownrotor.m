                                  %%% Inner Part %%%
r_down_in = zeros(1,X);           
seta_down_in = zeros(1,X)
for  i=1:X
     r_down_in(i)=((delta_r*i))/b;
     seta_down_in(i)= seta_o_down+ r_down_in(i)*seta_twist;
     v_c_in(i) = (v_root)+((v_ind_up-v_root)*r_down_in(i)/(y))
     lamda_climb_in(i)=v_c_in(i)/(omga*R);
     c_down_in(i) =(C_root)+(taper_slope*r_down_in(i)*b)

    %Lift Slope Condition
                         
     V_down_in_L(i) =sqrt (((omga* (((r_down_in(i))*b)+R_cut))^2)+((v_c_in(i))^2))
     Re_down_in_L(i) = (roh*V_down_in_L(i)*c_down_in(i))/Mioh
    if seta_down_in(i) <= 2*(pi/180)
         if Re_down_in_L(i) < 70000
              Cl_alpha_down_in(i) = ((0.2325/70000)* Re_down_in_L(i))*180/pi
         else  Cl_alpha_down_in(i) =((0.21+(20000*0.7/Re_down_in_L(i)))/2)*180/pi
         end

     elseif  2*(pi/180) < seta_down_in(i) <= 4*(pi/180)
         if Re_down_in_L(i) < 70000
              Cl_alpha_down_in(i)= ((0.152/70000)* Re_down_in_L(i))*180/pi
         elseif 70000 <= Re_down_in_L(i) <= 150000 
             Cl_alpha_down_in(i)=0.152*180/pi
         else Cl_alpha_down_in(i)=((0.21+(20000*0.7/Re_down_in_L(i)))/2)*180/pi
         end
         
     else Cl_alpha_down_in(i)= Cl_alpha

     end  
    
    %Lamda_down_in(i)=sqrt ((((sigma*Cl_alpha_down_in(i))/16)-(lamda_climb_in(i)/2))^2+(sigma*Cl_alpha_down_in(i)*seta_down_in(i)*r_down_in(i)/8))-((sigma*Cl_alpha_down_in(i)/16)-(lamda_climb_in(i)/2)); 
    Lamda_down_in(i)=((Cl_alpha_down_in(i)*sigma)/16)*((sqrt (1+((32*seta_down_in(i)*r_down_in(i))/(Cl_alpha_down_in(i)*sigma))))-1);
    B_down_in(i) = 1-(Lamda_down_in(i)*sqrt(2))
    delta_Ct_down_in(i)=((Cl_alpha_down_in(i)*sigma)/2)*((seta_down_in(i)*(r_down_in(i)*B_down_in(i))^2)-((Lamda_down_in(i))*r_down_in(i)/B_down_in(i)))*(delta_r/b);
    %Lamda_D(i)=((Cl_alpha_down_in(i)*sigma)/16)*((sqrt (1+((32*seta_down_in(i)*r_down_in(i))/(Cl_alpha_down_in(i)*sigma))))-1);
    Lamda_D(i)= Lamda_down_in(i)-lamda_climb_in(i)
    Cl_down_in(i)=Cl_alpha_down_in(i)*(seta_down_in(i)-( Lamda_down_in(i)/r_down_in(i))); 
    V_down_in(i)=sqrt (((omga* (((r_down_in(i))*b)+R_cut))^2)+((Lamda_down_in(i)*omga* (((r_down_in(i))*b)+R_cut)))^2)
    Re_down_in(i)= (roh*V_down_in(i)*c_down_in(i))/Mioh
    Cd_0_down_in(i)=2/((Re_down_in(i)*0.75/1000)+35)
    cd_downin_1_o(i) = (pi/180)*(Lamda_down_in(i))/r_down_in(i)
    cd_downin_2_o(i) = seta_down_in(i)-cd_downin_1_o(i)
    Cd_o_down_in(i) = (Cd_0_down_in(i)+(0.02*cd_downin_2_o(i))+(0.4*cd_downin_2_o(i)^2))*(delta_r/b);
    CQ_o_down_in(i)=(k_br*sigma/4)*Cd_o_down_in(i)
    seta_o_down=seta_down_in(i);
end