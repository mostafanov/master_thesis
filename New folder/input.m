%input
sigma=0.042;                                                             % solidaty ratio0.042 0.027
seta_o_up=8*(pi/180)
H=0.3
delta_seta =4*(pi/180)
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