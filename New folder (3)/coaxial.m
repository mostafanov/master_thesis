%input
sigma=0.042/2;                                                             % solidaty ratio0.042 0.027
seta_o_up=2*(pi/180)
H=0.2 
delta_seta =2*(pi/180)
omga= 130
C_root = 0.062                                                             % 0.062 0.294
taper_ratio = 0.637                                                      % 0.637 0.388 
R=  0.762  %3.81 ;                                                            %rotor redius
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

%B_loss=(1-((C_root*(1+(0.7*taper_ratio)))/(1.5*R)))
B_loss=1

K =1.15


if H < 0.35
K_int_upstream =1.4
else
    K_int_upstream =1.21
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
     Cl_up(i)=Cl_alpha_up(i)*(seta_up(i)-(Lamda_up(i)/(r_up(i)*B_up(i))));
     %V_up(i)= sqrt (((omga* (((r_up(i))*b)+R_cut))^2)+((Lamda_up(i)*omga* (((r_up(i))*b)+R_cut))^2))
     V_up(i)= (omga* (((r_up(i))*b)+R_cut))
     Re_up(i)= (roh*V_up(i)*c_up(i))/Mioh
     Cd_0_up(i)=2/((Re_up(i)*0.75/1000)+35)
     
     cdup_1_o (i) =(Lamda_up(i)/(r_up(i)))*pi/180
     %cdup_2_o (i) =(seta_up(i)-cdup_1_o (i))
     cdup_2_o (i) =(seta_up(i))
     Cd_o_up(i) = (Cd_0_up(i)+(0.02*cdup_2_o (i))+(0.4*(cdup_2_o (i))^2));
     CQ_o_up(i)=(sigma/8)*Cd_o_up(i)
     seta_o_up = seta_up(i);
 end
 total_delta_Ct_up =(sum (delta_Ct_up))
 Ct_up=total_delta_Ct_up*Nb*B_loss
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

 
                             %%%%%%%%%%% Downstream rotor %%%%%%%%%%%

for k=1:NN

    CT_up(k)=CT_up_old
    Lamda_D=Lamda_D_old
    V_D= Lamda_D*(omga*R)
    v_center=v_center_old
    v_ind_up_co= v_ind_up_old
   
    v_ind = v_ind_up_co
    Lamda_2_co =Lamda_2_s_old
    
   

%%%Y_tip=Y_tip_s
V_S_Ratio= sqrt(Lamda_2_s)
V_CO_Ratio= sqrt((sqrt (Lamda_2_co^2)))
Y_tip_co= sqrt ((Y_tip_s*V_S_Ratio/V_CO_Ratio)^2)
if  0.702 < Y_tip_co & Y_tip_co < 1
    Y_tip=Y_tip_co
else
Y_tip=Y_tip_s
end

%Elements Inside the Slip Stream Region
x  =(b-(1-Y_tip)*R)/b
X = round(x*N)
y =X/N

%Velocity Profile at the downstream rotor plane

                                       
if y < 1
v_center_1_co= v_ind_up_s_0-v_ind_up_co
Lamda_Center_1= v_center_1_co/(omga*R)
v_center_2_co= (0.5*v_ind_up_s_0*((1-(Y_tip^3))/(Y_tip^3)))-(V_D*0.5)
Lamda_Center_2= v_center_2_co/(omga*R)
v_center = v_center_1_co+v_center_2_co
v_root = v_center+((v_ind-v_center)*R_cut_factor/Y_tip)
Lamda_Center = v_center/(omga*R)

%Interferance Factor 1.41 when H=0 and 1.16 at venna 
K_int_downstream_slope = (0.99-0.72)/(1.41-1.16)
K_int_downstream  =(1.41-(K_int_downstream_slope*0.99))+(K_int_downstream_slope*Y_tip)

                                  %%% Inner Part %%%
r_down_in = zeros(1,X);           
seta_down_in = zeros(1,X)
for  i=1:X
     r_down_in(i)=((delta_r*i))/b;
     seta_down_in(i)= seta_o_down+ r_down_in(i)*seta_twist;
     v_c_in(i) = (v_root)+((v_ind-v_root)*r_down_in(i)/(y))
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
    
    Lamda_down_in(i)=sqrt ((((sigma*Cl_alpha_down_in(i))/16)-(lamda_climb_in(i)/2))^2+(sigma*Cl_alpha_down_in(i)*seta_down_in(i)*r_down_in(i)/8))-((sigma*Cl_alpha_down_in(i)/16)-(lamda_climb_in(i)/2)); 
    
    B_down_in(i) = 1-(Lamda_down_in(i)*sqrt(2))
    delta_Ct_down_in(i)=((Cl_alpha_down_in(i)*sigma)/2)*((seta_down_in(i)*(r_down_in(i)*B_down_in(i))^2)-((Lamda_down_in(i))*r_down_in(i)/B_down_in(i)))*(delta_r/b);
    Lamda_D(i)=((Cl_alpha_down_in(i)*sigma)/16)*((sqrt (1+((32*seta_down_in(i)*r_down_in(i))/(Cl_alpha_down_in(i)*sigma))))-1);
    %Lamda_D(i)= Lamda_down_in(i)-lamda_climb_in(i)
    Cl_down_in(i)=Cl_alpha_down_in(i)*(seta_down_in(i)-( Lamda_down_in(i)/r_down_in(i))); 
    %V_down_in(i)=sqrt (((omga* (((r_down_in(i))*b)+R_cut))^2)+((Lamda_down_in(i)*omga* (((r_down_in(i))*b)+R_cut)))^2)
    V_down_in(i)=(omga* (((r_down_in(i))*b)+R_cut))
    Re_down_in(i)= (roh*V_down_in(i)*c_down_in(i))/Mioh
    Cd_0_down_in(i)=2/((Re_down_in(i)*0.75/1000)+35)
    cd_downin_1_o(i) = (pi/180)*(Lamda_down_in(i))/r_down_in(i)
    %cd_downin_2_o(i) = seta_down_in(i)-cd_downin_1_o(i)
    cd_downin_2_o(i) = seta_down_in(i)
    Cd_o_down_in(i) = (Cd_0_down_in(i)+(0.02*cd_downin_2_o(i))+(0.4*cd_downin_2_o(i)^2))*(delta_r/b);
    CQ_o_down_in(i)=(k_br*sigma/4)*Cd_o_down_in(i)
    seta_o_down=seta_down_in(i);
end
                                           %%% Outer Part %%%

r_down_out = zeros(X,N);           
seta_down_out = zeros(X,N)
for i= (X:N)
    r_down_out(i)=((delta_r*i))/b;
    seta_down_out(i)= seta_o_down+ r_down_out(i)*seta_twist;
    v_c_out(i)= v_root+((v_ind-v_root)*(r_down_out(i))/((y)))
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
    %V_down_out(i)=sqrt (((omga* (((r_down_out(i))*b)+R_cut))^2)+((Lamda_down_out(i)*omga* (((r_down_out(i))*b)+R_cut)))^2)
    V_down_out(i)=(omga* (((r_down_out(i))*b)+R_cut))
    Re_down_out(i)= (roh*V_down_out(i)*c_down_out(i))/Mioh
    Cd_0_down_out(i)=2/((Re_down_out(i)*0.75/1000)+35)
    cd_downout_1_o(i) = (Lamda_down_out(i)/r_down_out(i))*pi/180
    %cd_downout_2_o(i)= seta_down_out(i)-cd_downout_1_o(i)
    cd_downout_2_o(i)= seta_down_out(i)
    Cd_o_down_out(i) = (Cd_0_down_out(i)+(0.02*cd_downout_2_o(i))+(0.4*(cd_downout_2_o(i))^2))*(delta_r/b)
    CQ_o_down_out(i)=(k_br*sigma/4)*Cd_o_down_out(i)
    seta_o_down=seta_down_out(i);
 end
total_delta_Ct_down_in= sum (delta_Ct_down_in)
total_delta_Ct_down_out=sum (delta_Ct_down_out)
Ct_down =(total_delta_Ct_down_in+total_delta_Ct_down_out)*Nb*B_loss
CT_down= Ct_down/K_int_downstream
delta_CQ_i_down_in= (delta_Ct_down_in.*Lamda_down_in)
CQ_i_down_in= sum (delta_CQ_i_down_in)
delta_CQ_i_down_out = (delta_Ct_down_out.*Lamda_down_out)
CQ_i_down_out= sum (delta_CQ_i_down_out)
total_CQ_o_down_in= sum (CQ_o_down_in)
total_CQ_o_down_out= sum (CQ_o_down_out)
B_loss_down= 1-(sqrt(CT_down/2))
cQ_down =((CQ_i_down_in+CQ_i_down_out)*(K)*K_int_downstream*Nb/B_loss_down)+((total_CQ_o_down_in+total_CQ_o_down_out)*Nb)
                                 
                                    %%% upstream Thrust correction%%%
Lamda_2= (sum (Lamda_down_in))*(delta_r/b)*Nb
%Lamda_2= (sum (Lamda_D)+sum (lamda_climb_in))*(delta_r/b)*Nb
Lamda_D_tip=Lamda_D(X)
lamda_Tip_relation = (Lamda_ind_up_s-Lamda_Center_1)/Lamda_ind_up_s
%Ct_downt_relation= (sqrt(Ct_up/2))*((3*Y_tip^2)-1)/(2*Y_tip^2)
CT_up_new_co =Ct_up*(lamda_Tip_relation)^2

if CT_up_new_co > Ct_up
    CT_up_new (k)=  Ct_up
       
else
    CT_up_new (k)=CT_up_new_co
end

elseif y == 1 & H < z_venna 
 v_center_1_co= v_ind_up_s_0-v_ind
v_center_2_co= (0.5*v_ind_up_s_0*((1-(Y_tip^3))/(Y_tip^3)))-(V_D*0.5)
v_center = v_center_1_co+v_center_2_co
  Lamda_Center = v_center/(omga*R)
  v_root= v_center+((v_ind-v_center)*R_cut_factor/Y_tip)
  K_int_downstream = 1.41
  for i=1:N
    r_down(i)=((delta_r*i))/b;
    seta_down(i)= seta_o_down+  r_down(i)*seta_twist;
    v_c(i)= v_root+((v_ind-v_root)*r_down(i))
    lamda_climb(i)=v_c(i)/(omga*R);
    c_down(i)=(C_root)+(taper_slope*r_down(i)*R)

                                        %%Lift slope%%
    V_down_L(i)=sqrt ((omga* (((r_down(i))*b)+R_cut))^2+(v_c(i))^2)
    Re_down_L(i)= (roh*V_down_L(i)*c_down(i))/Mioh                                    
    
     if seta_down(i) <= 2*(pi/180)
         if Re_down_T(i) < 70000
             Cl_alpha_down(i)= ((0.2325/70000)* Re_down_T(i))*180/pi
         else Cl_alpha_down(i)=((0.21+(20000*0.7/Re_down_T(i)))/2)*180/pi
         end

     elseif  2*(pi/180)< seta_down(i) <= 4*(pi/180)
         if Re_down_L(i) < 70000
             Cl_alpha_down(i)= ((0.152/70000)* Re_down_L(i))*180/pi
         elseif 70000 <= Re_down_L(i) <= 150000 
             Cl_alpha_down(i)=0.152*180/pi
         else Cl_alpha_down(i)=((0.21+(20000*0.7/Re_down_L(i)))/2)*180/pi
         end
         
     else Cl_alpha_down(i)= Cl_alpha

     end  

    Lamda_down(i)=sqrt ((((sigma*Cl_alpha_down(i))/16)-(lamda_climb(i)/2))^2+(sigma*Cl_alpha_down(i)*seta_down(i)*r_down(i)/8))-((sigma*Cl_alpha_down(i)/16)-(lamda_climb(i)/2)); 
    B_down(i)= 1-(Lamda_down(i)*sqrt(2))
    delta_Ct_down(i)=((Cl_alpha_down(i)*sigma)/2)*((seta_down(i)*(r_down(i)*B_down(i))^2)-((Lamda_down(i))*r_down(i)/B_down(i)))*(delta_r/b)
    Lamda_D(i)=((Cl_alpha_down(i)*sigma)/16)*((sqrt (1+((32*seta_down(i)*r_down(i))/(Cl_alpha_down(i)*sigma))))-1);
    %Lamda_D(i)= Lamda_down(i)-lamda_climb(i)
    Cl_down(i)=Cl_alpha_down(i)*(seta_down(i)-( Lamda_down(i)/r_down(i)));
    %V_down(i)=sqrt (((omga* (((r_down(i))*b)+R_cut))^2)+((Lamda_down(i)*omga* (((r_down(i))*b)+R_cut)))^2)
    V_down(i)=(omga* (((r_down(i))*b)+R_cut))
    Re_down(i)= (roh*V_down(i)*c_down(i))/Mioh
    Cd_0_down(i)=2/((Re_down(i)*0.75/1000)+35)
    cd_down_1_o(i) = (pi/180)*(Lamda_down(i))/r_down(i)
    %cd_down_2_o(i) = seta_down_in(i)-cd_down_1_o(i)
    cd_down_2_o(i) = seta_down_in(i)
    Cd_o_down_(i) = (Cd_0_down_(i)+(0.02*cd_down_2_o(i))+(0.4*cd_down_2_o(i)^2))*(delta_r/b);
    CQ_o_down(i)=(k_br*sigma/4)*Cd_o_down(i)
    seta_o_down=seta_down(i);
  end
 total_Ct_down =(sum (delta_Ct_down))
 Ct_down=total_Ct_down*Nb
 
 CT_down= Ct_down/K_int_downstream
 delta_CQ_i_down=delta_Ct_down.*Lamda_down/B_down
 CQ_i_down = sum (delta_CQ_i_down)
 total_CQ_o_down= sum (CQ_o_down)
 B_loss_down= 1-(sqrt(CT_down/2)/Nb)
 cQ_down=(CQ_i_down*(K)*K_int_downstream/B_loss_down)+(total_CQ_o_down*Nb)
 Lamda_D_tip=Lamda_D(N)
                                    %%% upstream Thrust correction%%%
Lamda_2= (sum (Lamda_down))*(delta_r/b)*Nb
%Lamda_2= (sum (Lamda_D)+sum (lamda_climb))*(delta_r/b)*Nb
lamda_Tip_relation = (Lamda_ind_up_s-Lamda_Center_1)/Lamda_ind_up_s
%Ct_downt_relation= (sqrt(Ct_up/2))*((3*Y_tip^2)-1)/(2*Y_tip^2)
CT_up_new_co =Ct_up*(lamda_Tip_relation)^2
if CT_up_new_co > Ct_up
    CT_up_new (k)=  Ct_up
else
    CT_up_new (k)=CT_up_new_co
end
else 
v_center= 0
v_root= 0
K_int_downstream=1
for i=1:N
     r_down(i)=((delta_r*i))/b;
     seta_down(i)= seta_o_down+ r_down(i)*seta_twist;
     Lamda_down(i)=((Cl_alpha*sigma)/16)*((sqrt (1+((32*seta_down(i)*r_down(i))/(Cl_alpha*sigma))))-1);
     B_down(i) = 1-(Lamda_down(i)*sqrt(2)/1)
     delta_Ct_down(i)=((Cl_alpha*sigma)/2)*((seta_down(i)*(r_down(i)*B_down(i))^2)-(Lamda_down(i)*r_down(i)/B_down(i)))*(delta_r/R) ;
     Cl_down(i)=Cl_alpha*(seta_down(i)-(Lamda_down(i)/r_down(i)));
     c_down(i)=(C_root)+(taper_slope*r_down(i)*b)   
     %V_down(i)= sqrt ((omga* (((r_down(i))*b)+R_cut))^2+(Lamda_down(i)*omga* (((r_down(i))*b)+R_cut))^2)
     V_down(i)= (omga* (((r_down(i))*b)+R_cut))
     Re_down(i)= (roh*V_down(i)*c_down(i))/Mioh
     Cd_0_down(i)=2/((Re_down(i)*0.75/1000)+35)
     cd_down_1_o(i)= (Lamda_down(i)/r_down(i))*pi/180
     cd_down_2_o(i) = seta_down(i)-cd_down_1_o(i)
     cd_down_2_o(i) = seta_down(i)
     Cd_o_down(i) = (Cd_0_down(i)+(0.02*cd_down_2_o(i))+(0.4*(cd_down_2_o(i))^2))*(delta_r/b);
     CQ_o_down(i)=(k_br*sigma/4)*Cd_o_down(i)
     seta_o_up = seta_up(i);
 end
 total_Ct_down =(sum (delta_Ct_down))
 Ct_down=total_Ct_down*Nb

 CT_down= Ct_down
 delta_CQ_i_down=delta_Ct_down.*Lamda_down/B_down
 CQ_i_down = sum (delta_CQ_i_down)
 total_CQ_o_down= sum (CQ_o_down)
 B_loss_down= 1-(sqrt(CT_down/2)/Nb)
 cQ_down=(CQ_i_down*(K)*K_int_downstream*Nb/B_loss_down)+(total_CQ_o_down*Nb)
CT_up_new(k) = Ct_up

end

CT_up_old= CT_up_new(k)

Lamda_D_old = Lamda_D_tip
v_center_old = v_center

v_ind_up_old=sqrt(CT_up_new(k)/2)*omga*R
Lamda_2_s_old=Lamda_2
v_center_2_co_old =v_center_2_co
end


                                          %%%%%coaxial parameter%%%%%
                                        
 CT_up_new_f= CT_up_new(NN)                                    

 CT = (CT_down+ CT_up_new_f)
 cp = cQ_down+ cQ_up

 %figure (101); 
 %plot(r_up,Lamda_up,'ro')
 %plot(r_down_in,Lamda_down_in)
 %plot(r_down_out,v_c_out)
 %grid on;