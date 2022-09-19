%input
sigma=0.042;                                                               % solidaty ratio
seta_o_up=2*(pi/180);                                                      %pitch angle upstream rotor
H=0.2;                                                                     %spaceing

delta_seta=2*(pi/180)                                                      %pitch angle difference
seta_twist=0*(pi/180);                                                     %twist angle
omga=130;
roh= 1.225
Mioh= 0.0000181206
C_root= 0.062
taper_ratio= 0.667                                                             %0.645
R=  .762;                                                                 %rotor redius
Nb=2
R_cut_factor=0.168
%parameter

Cl_alpha=2*pi;                                                             % airfoil coefficient of lift slope
seta_o_down=seta_o_up+delta_seta;                                          %pitch angle downstream rotor
R_cut=R_cut_factor*R;                                                      %rotor cut off
N=400                                                                       %descritization elements
b= R-R_cut                                                                 %blade span
delta_r=(b/N)                               
taper_slope= ((C_root*taper_ratio)-C_root)/b
B=(1-((C_root*(1+(0.7*taper_ratio)))/(1.5*R)))


K=1.15
K_int=1.16
%k_br=0.6
k_br=0.5
%k_br=0.85

                                   %%%%%%%upstream rotor %%%%%%%
 for i=1:N
     r_up(i)=((delta_r*i))/b;
     seta_up(i)= seta_o_up+ r_up(i)*seta_twist;
     Lamda_up(i)=((Cl_alpha*sigma)/16)*((sqrt (1+((32*seta_up(i)*r_up(i))/(Cl_alpha*sigma))))-1);
     B_up(i)=1-((1-B)*((r_up(i))^2)*1.5)
     delta_Ct_up(i)=((Cl_alpha*sigma)/2)*((seta_up(i)*(r_up(i)* B_up(i))^2)-(Lamda_up(i)*r_up(i)* B_up(i)))*delta_r ;
     Cl_up(i)=Cl_alpha*(seta_up(i)-(Lamda_up(i)/r_up(i)));
     
     c_up(i)=(C_root)+(taper_slope*r_up(i)*R)
     %V_up(i)=  omga* (((r_up(i))*b)+R_cut)
     V_up(i)= sqrt ((omga* (((r_up(i))*b)+R_cut))^2+(Lamda_up(i)*omga* (((r_up(i))*b)+R_cut))^2)
     Re_up(i)= (roh*V_up(i)*c_up(i))/Mioh
     Cd_0_up(i)=2/((Re_up(i)*0.75/1000)+35)
     Cd_o_up(i) = (Cd_0_up(i)+(0.02*(seta_up(i)-(Lamda_up(i)/r_up(i)))*pi/180)+(0.4*((seta_up(i)-(Lamda_up(i)/r_up(i)))*pi/180)^2));
     Cp_o_up(i)=(k_br*sigma/4)*Cd_o_up(i)
     seta_o_up = seta_up(i);
 end
 s =(sum (delta_Ct_up))
 Ct_up=s*Nb
 f=delta_Ct_up.*Lamda_up
 d = sum (f)
 C_up= sum (Cp_o_up)
 cQ_up=(d*2*K*K_int)+(C_up*delta_r*Nb) 

                                  %%%%%%% Downstream rotor %%%%%%%


                                         
K_1_tip = 0.25*((Ct_up/sigma)+(0.001*seta_twist))
Zeta_passage = 360/Nb
Z_Tip_passage = K_1_tip*Zeta_passage
K_2_tip = (1+(0.01*seta_twist)) * sqrt(Ct_up)
Zeta_age= ((H-Z_Tip_passage)/K_2_tip)+Zeta_passage
A_tip = 0.145+(27*Ct_up)
if Z_Tip_passage > H
    Y_tip= 0.78+(1-0.78)*exp(-A_tip*(H/K_1_tip))
else Y_tip= 0.78+(1-0.78)*exp(-A_tip*Zeta_age)
end
%beta=30*(pi/180);
%beta=(45-(45*4/(4+(seta_o_up/(pi/180)))))*(pi/180)
%x=1-(H*tan(beta));
X= round(Y_tip*N)
y=X/N

%v_ind_up=(Lamda_up(N))*omga*((r_up(N)*b)+R_cut)
 v_ind_up=(Lamda_up(N))*omga*R
%%v_c= v_ind_up*(1+(((H)^2)/4))
%%lamda_climb=v_c/(omga*R);                                                  %velocity ratio for climb                                         

v_center= v_ind_up*(((1-R_cut_factor)/(Y_tip^2))-1)
v_root= v_center+((v_ind_up-v_center)*R_cut_factor/Y_tip)

                                       %inner part
for i=1:y*N
    r_down_in(i)=((delta_r*i))/b;
    seta_down_in(i)= seta_o_down+  r_down_in(i)*seta_twist;

    v_c(i)= v_root+((v_ind_up-v_root)*r_down_in(i)/Y_tip)
    lamda_climb(i)=v_c(i)/(omga*R);

    Lamda_down_in(i)=sqrt ((((sigma*Cl_alpha)/16)-(lamda_climb(i)/2))^2+(sigma*Cl_alpha*seta_down_in(i)*r_down_in(i)/8))-((sigma*Cl_alpha/16)-(lamda_climb(i)/2)); 
    B_down_in(i)=1-((1-B)*((r_down_in(i))^2)*1.5)
    delta_Ct_down_in(i)=((Cl_alpha*sigma)/2)*((seta_down_in(i)*(r_down_in(i)*B_down_in(i))^2)-(Lamda_down_in(i)*r_down_in(i)*B_down_in(i)))*delta_r;
    Cl_down_in(i)=Cl_alpha*(seta_down_in(i)-( Lamda_down_in(i)/r_down_in(i)));
    c_down_in(i)=(C_root)+(taper_slope*r_down_in(i)*R)
  
    %V_down_in(i)=  omga* (((r_down_in(i))*b)+R_cut)
    V_down_in(i)=sqrt ((omga* (((r_down_in(i))*b)+R_cut))^2+(Lamda_down_in(i)*omga* (((r_down_in(i))*b)+R_cut))^2)
    Re_down_in(i)= (roh*V_down_in(i)*c_down_in(i))/Mioh
    Cd_0_down_in(i)=2/((Re_down_in(i)*0.75/1000)+35)
    Cd_o_down_in(i) = (Cd_0_down_in(i)+(0.02*pi/180*(seta_down_in(i)-(Lamda_down_in(i)/r_down_in(i))))+(0.4*((seta_down_in(i)-(Lamda_down_in(i)/r_down_in(i)))*pi/180)^2));
    Cp_o_down_in(i)=(k_br*sigma/4)*Cd_o_down_in(i)
    seta_o_down=seta_down_in(i);
end
                                           %outer part

 
for i=y*N:N
    r_down_out(i)=((delta_r*i))/b;
    seta_down_out(i)= seta_o_down+ r_down_out(i)*seta_twist;
    Lamda_down_out(i)=((Cl_alpha*sigma)/16)*((sqrt (1+((32*seta_down_out(i)*r_down_out(i))/(Cl_alpha*sigma))))-1);
    B_down_out(i)=1-((1-B)*((r_down_out(i))^2)*1.5)
    delta_Ct_down_out(i)=((Cl_alpha*sigma)/2)*((seta_down_out(i)*(r_down_out(i)*B_down_out(i))^2)-(Lamda_down_out(i)*r_down_out(i)*B_down_out(i)))*delta_r ;
    Cl_down_out(i)=Cl_alpha*(seta_down_out(i)-( Lamda_down_out(i)/r_down_out(i)));
     c_down_out(i)=(C_root)+(taper_slope*r_down_out(i)*R)

     %V_down_out(i)=  omga* (((r_down_out(i))*b)+R_cut)
     V_down_out(i)=sqrt ((omga* (((r_down_out(i))*b)+R_cut))^2+(Lamda_down_out(i)*omga* (((r_down_out(i))*b)+R_cut))^2)
    Re_down_out(i)= (roh*V_down_out(i)*c_down_out(i))/Mioh
    Cd_0_down_out(i)=2/((Re_down_out(i)*0.75/1000)+35)
    Cd_o_down_out(i) = (Cd_0_down_out(i)+(0.02*pi/180*(seta_down_out(i)-(Lamda_down_out(i)/r_down_out(i))))+(0.4*((seta_down_out(i)-(Lamda_down_out(i)/r_down_out(i)))*pi/180)^2));
    Cp_o_down_out(i)=(k_br*sigma/4)*Cd_o_down_out(i)
    seta_o_down = seta_down_out(i);
end
a= sum (delta_Ct_down_in)
b=sum (delta_Ct_down_out)
Ct_down =(a+b)*Nb
q= (delta_Ct_down_in.*Lamda_down_in)
g= sum (q)
v = (delta_Ct_down_out.*Lamda_down_out)
l= sum (v)
C_down_in= sum (Cp_o_down_in)
C_down_out= sum (Cp_o_down_out)
cQ_down =((g+l)*2*K*K_int)+((C_down_in+C_down_out)*delta_r*Nb)


                                          %%%%%coaxial parameter%%%%%
 
 CT = Ct_down+ Ct_up
 cp = cQ_down+ cQ_up
 %CT_sigma=CT/((sigma)^2)
 %figure (101);
 
 %plot(r_up,Lamda_up,'ro')
 %plot(r_up/R,Cl_up,'ro')
 %plot(r_down_in/R,Cl_down_in,'ro')
 %plot(r_down_out/R,Cl_down_out,'ro')
 %plot(r_up/R,a_up*r_up,'ro')
 %grid on;