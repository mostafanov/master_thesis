Rotor_Solidity=0.042;                                                      % solidaty ratio for the upstream rotor 0.042 0.027
Up_root_Pitch=2;                                                           % pitch angle of upstream rotor at the rotor
H=0.2;                                                                     % Non-dimensional axial distance between the two rotors for a coaxial
del_pitch=2;                                                               % difference between the pitch angle at the root for upstream and downstream rotor
rpm = 1240;                                                                % revoltuions per min
C_root = 0.062                                                             % chord length at the root in m 
taper_ratio = 0.637                                                        % taper ratio 
R=  0.762;                                                                 % rotor radius in m
Nb =2                                                                      % Number of blades for single rotor
R_cut_factor = 0.167                                                       % Dimensionless distance between root and hub of the propeller 0.167 0.13
blade_twist=0;                                                             % twist angle
roh = 1.225                                                                % freestream density kg/m3
Mioh = 0.0000181206 

% Dynamic viscosity  kg/m/s
out = test_sum(Rotor_Solidity, Up_root_Pitch, H, del_pitch, rpm, C_root, taper_ratio, R, R_cut_factor, blade_twist, Nb, roh, Mioh)
