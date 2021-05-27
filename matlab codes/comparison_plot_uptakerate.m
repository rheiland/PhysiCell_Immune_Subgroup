
	
Vvoxel = 20*20*20;

m_half = 10;
rho_max = 0.125;
rho_half = 0.125;
u_Evirus = 0.0001;
IC_50_IFN = 3;
    
IFN_internal = 0;
u_Evirus_IFN = u_Evirus*(1-IFN_internal/(IC_50_IFN+IFN_internal));


rho_virus_vec = linspace(0,1,100);
m_i_vec = linspace(0,10,100);

for ii = 1:100
    rho_virus = rho_virus_vec(ii);
    for jj = 1:100
        m_i = m_i_vec(jj);
        if rho_virus<1e-6
            uptake_rates(ii,jj) = 0;
        else
            uptake_rates(ii,jj) = u_Evirus_IFN*(rho_virus/(rho_half+rho_virus))*(m_half/(m_i+m_half));
        end

    end
end

[m_i_mat,rho_virus_mat] = meshgrid(m_i_vec,rho_virus_vec);

figure
hold on 
surf(m_i_mat,rho_virus_mat, uptake_rates,'EdgeColor','none')
colorbar
ylabel('\rho_{virus}')
xlabel('m_i')
set(gca,'FontSize',18)