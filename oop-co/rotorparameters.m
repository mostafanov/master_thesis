function [sigma, seta_o_up, delta_seta, seta_twist, omga, Cl_alpha, seta_o_down, R_cut, b, delta_r, taper_slope, K_int_upstream, Nb, roh, Mioh, N, NN, K, k_br] = rotorparameters(Rotor_Solidity, Up_root_Pitch, H, del_pitch, rpm, C_root, taper_ratio, R, R_cut_factor)
    function [Rotor_Solidity, Up_root_Pitch, H, del_pitch, rpm, C_root, taper_ratio, R, R_cut_factor] = userinput(Rotor_Solidity, Up_root_Pitch, H, del_pitch, rpm, C_root, taper_ratio, R, R_cut_factor)
    end
sigma=Rotor_Solidity/2;                                      
seta_o_up=Up_root_Pitch*(pi/180);
delta_seta = del_pitch*(pi/180);
seta_twist=blade_twist*(pi/180);
omga= rpm*2*(pi/60);                                                       % rad/s
Cl_alpha=2*pi;                                                             % airfoil coefficient of lift slope
seta_o_down=seta_o_up+delta_seta;                                          % pitch angle downstream rotor
R_cut=R_cut_factor*R;                                                      % rotor cut off
N =400                                                                     % descritization elements
NN=20                                                                      % Iteration count
b = R-R_cut                                                                % single blade span from root to tip excluding the hub
delta_r =(b/N)                                                             % span of each element                              
taper_slope = ((C_root*taper_ratio)-C_root)/b                              % slope of the chord changing with the span
K =1.15                                                                    % Correction factor for induced drag to calculate induced power
if H < 0.35
K_int_upstream =1.4                                                        % interference/interaction factor for coaxial
else
    K_int_upstream =1.21
end
k_br =0.5                                                                  % factor for profile drag

end