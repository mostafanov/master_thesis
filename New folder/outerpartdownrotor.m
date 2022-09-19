                                           %%% Outer Part %%%

r_down_out = zeros(X,N);           
seta_down_out = zeros(X,N)
for i= (X:N)
    r_down_out(i)=((delta_r*i))/b;
    seta_down_out(i)= seta_o_down+ r_down_out(i)*seta_twist;
    v_c_out(i)= v_root+((v_ind_up-v_root)*(r_down_out(i))/((y)))
    lamda_climb_out(i)=v_c_out(i)/(omga*R);
    c_down_out(i) =(C_root)+(taper_slope*r_down_out(i)*b)
                              
                           %%%lift slope condition
                         
    V_down_out_L(i)=sqrt (((omga* (((r_down_out(i))*b)+R_cut))^2)+(v_c_out(i))^2)
    Re_down_out_L(i)= (roh*V_down_out_L(i)*c_down_out(i))/Mioh
    if seta_down_out(i) <= 2*(pi/180)
         if Re_down_out_L(i) < 70000
             Cl_alpha_down_out(i)= ((0.2325/70000)* Re_down_out_L(i))*180/pi
         else Cl_alpha_down_out(i)=((0.21+(20000*0.7/Re_down_out_L(i)))/2)*180/pi
         end

     elseif  2*(pi/180)< seta_down_out(i) <= 4*(pi/180)
         if Re_down_out_L(i) < 70000
             Cl_alpha_down_out(i)= ((0.152/70000)* Re_down_out_L(i))*180/pi
         elseif 70000 <= Re_down_out_L(i) <= 150000 
             Cl_alpha_down_out(i)=0.152*180/pi
         else Cl_alpha_down_out(i)=((0.21+(20000*0.7/Re_down_out_L(i)))/2)*180/pi
         end
         
     else Cl_alpha_down_out(i)= Cl_alpha

     end  

    Lamda_down_out(i)=sqrt ((((sigma*Cl_alpha_down_out(i))/16)-(lamda_climb_out(i)/2))^2+(sigma*Cl_alpha_down_out(i)*seta_down_out(i)*r_down_out(i)/8))-((sigma*Cl_alpha_down_out(i)/16)-(lamda_climb_out(i)/2)); 
    %Lamda_down_out(i)=((Cl_alpha_down_out(i)*sigma)/16)*((sqrt (1+((32*seta_down_out(i)*r_down_out(i))/(Cl_alpha_down_out(i)*sigma))))-1);
    B_down_out(i)=1-(Lamda_down_out(i)*sqrt(2))
    delta_Ct_down_out(i)=((Cl_alpha_down_out(i)*sigma)/2)*((seta_down_out(i)*(r_down_out(i)*B_down_out(i))^2)-(Lamda_down_out(i)*r_down_out(i)/B_down_out(i)))*(delta_r/b); 
    Cl_down_out(i)=Cl_alpha_down_out(i)*(seta_down_out(i)-(Lamda_down_out(i)/r_down_out(i)));
    V_down_out(i)=sqrt (((omga* (((r_down_out(i))*b)+R_cut))^2)+((Lamda_down_out(i)*omga* (((r_down_out(i))*b)+R_cut)))^2)
    Re_down_out(i)= (roh*V_down_out(i)*c_down_out(i))/Mioh
    Cd_0_down_out(i)=2/((Re_down_out(i)*0.75/1000)+35)
    cd_downout_1_o(i) = (Lamda_down_out(i)/r_down_out(i))*pi/180
    cd_downout_2_o(i)= seta_down_out(i)-cd_downout_1_o(i)
    Cd_o_down_out(i) = (Cd_0_down_out(i)+(0.02*cd_downout_2_o(i))+(0.4*(cd_downout_2_o(i))^2))*(delta_r/b)
    CQ_o_down_out(i)=(k_br*sigma/4)*Cd_o_down_out(i)
    seta_o_down=seta_down_out(i);
 end