function A_dic = Array_dictionary(N, G)

    % make dictionary matrix from p.2954
    % AoA & AoD are assumed to be distributed independently 
    % and uniformly in [0,2pi]
    theta_g = (0:G-1)/G*2*pi; 
%     A_dic = sqrt(1/N)*exp(1i*2*pi*(0:N-1)'.*cos(theta_g));
    A_dic = sqrt(1/N)*exp(1i*pi*(0:N-1)'.*cos(theta_g));
end