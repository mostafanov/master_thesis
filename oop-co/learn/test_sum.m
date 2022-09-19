function [sigma, seta_o_up, delta_seta, seta_twist, omga, Cl_alpha, seta_o_down, R_cut, b, delta_r, taper_slope, K_int_upstream, Nb, roh, Mioh, N, NN, K, k_br] = test_sum(Rotor_Solidity, Up_root_Pitch, H, del_pitch, rpm, C_root, taper_ratio, R, R_cut_factor, blade_twist, Nb, roh, Mioh)

out.sigma=Rotor_Solidity/2;                                      
out.seta_o_up=Up_root_Pitch*(pi/180);
out.delta_seta = del_pitch*(pi/180);
out.seta_twist=blade_twist*(pi/180);
out.omga= rpm*2*(pi/60);                                                       % rad/s
out.Cl_alpha=2*pi;                                                             % airfoil coefficient of lift slope
out.seta_o_down=out.seta_o_up+delta_seta;                                          % pitch angle downstream rotor
out.R_cut=R_cut_factor*R;                                                      % rotor cut off
out.N =400                                                                     % descritization elements
out.NN=20                                                                      % Iteration count
out.b = R-R_cut                                                                % single blade span from root to tip excluding the hub
out.delta_r =(b/N)                                                             % span of each element                              
out.taper_slope = ((C_root*taper_ratio)-C_root)/b                              % slope of the chord changing with the span
out.K =1.15                                                                    % Correction factor for induced drag to calculate induced power
if H < 0.35
out.K_int_upstream =1.4                                                        % interference/interaction factor for coaxial
else
   out.K_int_upstream =1.21
end
out.k_br =0.5                                                                  % factor for profile drag

end