                             %%%%%%%%%%% Downstream rotor %%%%%%%%%%%

for k=1:NN

    CT_up(k)=CT_up_old
    Lamda_D=Lamda_D_old
    V_D= Lamda_D*(omga*R)
    v_center_c=v_center_old
    v_ind= v_ind_up_s

%%%Y_tip=Y_tip_s
V_S_Ratio= (v_center_s+v_ind_up_old)^1/3
V_CO_Ratio= (v_center_c+v_ind+V_D)^1/3

Y_tip=Y_tip_s*V_S_Ratio/V_CO_Ratio
%Elements Inside the Slip Stream Region
x  =(b-(1-Y_tip)*R)/b
X = round(x*N)
y =X/N

%Velocity Profile at the downstream rotor plane

Lmada_induced_Tip = sqrt(sqrt(CT_up(k)^2)/2)/Tip_Thrust_r_s 
v_ind_up =(Lmada_induced_Tip)*omga*R                                        
if y < 1
 
v_center = (v_ind_up_old/(Y_tip^3))-V_D-v_ind_up
v_root = v_center+((v_ind_up-v_center)*R_cut_factor/Y_tip)
Lamda_Center = v_center/(omga/R)

%Interferance Factor 1.41 when H=0 and 1.16 at venna 
K_int_downstream_slope = (0.99-0.72)/(1.41-1.16)
K_int_downstream  =(1.41-(K_int_downstream_slope*0.99))+(K_int_downstream_slope*Y_tip)
