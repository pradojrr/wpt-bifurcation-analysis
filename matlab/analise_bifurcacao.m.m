clear; clc; close all;

%% Parametros do projeto
f0 = 50e3;
w0 = 2*pi*f0;
Pout = 1000;
Vout = 60;
Rc = Vout^2 / Pout; % 3.6 ohms
k = 0.25;
Rp = 0.1; 
Rs = 0.1;

f = linspace(30e3, 70e3, 1000); 
w = 2*pi*f;

%% Fig 1 - SS Simetrico (Bifurcado)
Lp1 = 267e-6; 
Ls1 = 267e-6;
Vp1 = 350;

M1 = k * sqrt(Lp1 * Ls1); % Eq. 2
Cs1 = 1 / (w0^2 * Ls1);   % Eq. 27
Cp1 = 1 / (w0^2 * Lp1);   % Eq. 29

Zs1 = Rs + Rc + 1i*(w*Ls1 - 1./(w*Cs1));
Zr1 = (w.^2 * M1^2) ./ Zs1;
Zeq1 = Rp + 1i*(w*Lp1 - 1./(w*Cp1)) + Zr1; % Eq. 25

Ip1 = Vp1 ./ Zeq1;
Is1 = (1i*w*M1 .* Ip1) ./ Zs1;
Pout1 = (abs(Is1).^2) * Rc;
Pin1 = real(Vp1 * conj(Ip1));
eta1 = (Pout1 ./ Pin1) * 100;

plotar_2D(1, f, Zeq1, Ip1, Pout1, eta1, 'Fig 1 - SS Original');

%% Fig 2 - Superficies SS
[Lp_grid, Ls_grid] = meshgrid(linspace(50e-6, 1500e-6, 40), linspace(20e-6, 300e-6, 40));
M_grid = k * sqrt(Lp_grid .* Ls_grid);

Is_req = sqrt(Pout / Rc);
Ip_req = (Is_req * (Rs + Rc)) ./ (w0 * M_grid);
Zref_grid = (w0 * M_grid).^2 / (Rs + Rc);
Zin_grid = Rp + Zref_grid;
Vp_req_SS = Ip_req .* Zin_grid;

Ploss_SS = Rp * (Ip_req.^2) + Rs * (Is_req^2);
Eta_SS = (Pout ./ (Pout + Ploss_SS)) * 100;

plotar_3D(2, Ls_grid, Lp_grid, Vp_req_SS, Eta_SS, 'Fig 2 - Superficie SS', [0 800], [80 100]);

%% Fig 3 - SS Otimizado
Lp3 = 1450e-6; 
Ls3 = 60e-6;
Vp3 = 380;

M3 = k * sqrt(Lp3 * Ls3);
Cs3 = 1 / (w0^2 * Ls3);
Cp3 = 1 / (w0^2 * Lp3);

Zs3 = Rs + Rc + 1i*(w*Ls3 - 1./(w*Cs3));
Zr3 = (w.^2 * M3^2) ./ Zs3;
Zeq3 = Rp + 1i*(w*Lp3 - 1./(w*Cp3)) + Zr3;

Ip3 = Vp3 ./ Zeq3;
Is3 = (1i*w*M3 .* Ip3) ./ Zs3;
Pout3 = (abs(Is3).^2) * Rc;
Pin3 = real(Vp3 * conj(Ip3));
eta3 = (Pout3 ./ Pin3) * 100;

plotar_2D(3, f, Zeq3, Ip3, Pout3, eta3, 'Fig 3 - SS Otimizado');

%% Fig 4 - SP Simetrico
Lp4 = 400e-6; 
Ls4 = 400e-6;
Vp4 = 250;

M4 = k * sqrt(Lp4 * Ls4);
Cs4 = 1 / (w0^2 * Ls4);
Cp4 = Ls4 / (w0^2 * (Lp4*Ls4 - M4^2)); % Eq. 44

Z_tank4 = Rc ./ (1 + 1i*w*Rc*Cs4);
Zs4 = Rs + 1i*w*Ls4 + Z_tank4;
Zr4 = (w.^2 * M4^2) ./ Zs4;
Zeq4 = Rp + 1i*(w*Lp4 - 1./(w*Cp4)) + Zr4; % Eq. 40

Ip4 = Vp4 ./ Zeq4;
Is4 = (1i*w*M4 .* Ip4) ./ Zs4;
Vout4 = Is4 .* Z_tank4;
Pout4 = (abs(Vout4).^2) / Rc;
Pin4 = real(Vp4 * conj(Ip4));
eta4 = (Pout4 ./ Pin4) * 100;

plotar_2D(4, f, Zeq4, Ip4, Pout4, eta4, 'Fig 4 - SP Original');

%% Fig 5 - Superficies SP
Vp_req_SP = zeros(size(Lp_grid));
Eta_SP = zeros(size(Lp_grid));

for i = 1:size(Lp_grid,1)
    for j = 1:size(Lp_grid,2)
        Lp_val = Lp_grid(i,j); 
        Ls_val = Ls_grid(i,j);
        M_val = k * sqrt(Lp_val * Ls_val);
        Cs_val = 1 / (w0^2 * Ls_val);
        
        Vout_val = sqrt(Pout * Rc);
        I_Rc = Vout_val / Rc;
        I_Cs = Vout_val * 1i * w0 * Cs_val;
        Is_val = I_Rc + I_Cs;
        
        V_ind = Is_val * (Rs + 1i*w0*Ls_val) + Vout_val;
        Ip_val = V_ind / (1i * w0 * M_val);
        
        Ztank = Rc / (1 + 1i*w0*Rc*Cs_val);
        Zsec = Rs + 1i*w0*Ls_val + Ztank;
        Zref_val = (w0 * M_val)^2 / Zsec;
        
        Zin_res = Rp + real(Zref_val);
        Vp_req_SP(i,j) = abs(Ip_val) * Zin_res;
        
        Ploss = Rp * abs(Ip_val)^2 + Rs * abs(Is_val)^2;
        Eta_SP(i,j) = (Pout / (Pout + Ploss)) * 100;
    end
end

plotar_3D(5, Ls_grid, Lp_grid, Vp_req_SP, Eta_SP, 'Fig 5 - Superficie SP', [0 800], [30 100]);

%% Fig 6 - SP Otimizado
Lp6 = 800e-6; 
Ls6 = 25e-6;
Vp6 = 360;

M6 = k * sqrt(Lp6 * Ls6);
Cs6 = 1 / (w0^2 * Ls6);
Cp6 = Ls6 / (w0^2 * (Lp6*Ls6 - M6^2)); % Eq. 44

Z_tank6 = Rc ./ (1 + 1i*w*Rc*Cs6);
Zs6 = Rs + 1i*w*Ls6 + Z_tank6;
Zr6 = (w.^2 * M6^2) ./ Zs6;
Zeq6 = Rp + 1i*(w*Lp6 - 1./(w*Cp6)) + Zr6;

Ip6 = Vp6 ./ Zeq6;
Is6 = (1i*w*M6 .* Ip6) ./ Zs6;
Vout6 = Is6 .* Z_tank6;
Pout6 = (abs(Vout6).^2) / Rc;
Pin6 = real(Vp6 * conj(Ip6));
eta6 = (Pout6 ./ Pin6) * 100;

plotar_2D(6, f, Zeq6, Ip6, Pout6, eta6, 'Fig 6 - SP Otimizado');

%% Funcoes de plot
function plotar_2D(num, f, Zeq, Ip, Pout, eta, nome)
    figure(num);
    set(gcf, 'Name', nome, 'Color', 'k', 'Position', [100, 100, 850, 650]);
    fk = f / 1000;
    
    subplot(2,2,1);
    set(gca, 'Color', [0.08 0.08 0.08], 'XColor', 'w', 'YColor', 'w');
    yyaxis left;
    plot(fk, abs(Zeq), 'b-', 'LineWidth', 1.5);
    ylabel('|Z_{eq}| (\Omega)');
    yyaxis right;
    plot(fk, abs(Ip), 'r--', 'LineWidth', 1.5);
    ylabel('Corrente de entrada (A)');
    xlabel('Frequência (kHz)');
    title('(a)', 'Color', 'w');
    grid on;
    
    subplot(2,2,2);
    set(gca, 'Color', [0.08 0.08 0.08], 'XColor', 'w', 'YColor', 'w');
    plot(fk, rad2deg(angle(Zeq)), 'w-', 'LineWidth', 1.5);
    ylabel('Ângulo da Impedância (graus)');
    xlabel('Frequência (kHz)');
    ylim([-90 90]);
    title('(b)', 'Color', 'w');
    grid on;
    
    subplot(2,2,3);
    set(gca, 'Color', [0.08 0.08 0.08], 'XColor', 'w', 'YColor', 'w');
    plot(fk, Pout, 'g-', 'LineWidth', 1.5);
    ylabel('Potência de saída (W)');
    xlabel('Frequência (kHz)');
    title('(c)', 'Color', 'w');
    grid on;
    
    subplot(2,2,4);
    set(gca, 'Color', [0.08 0.08 0.08], 'XColor', 'w', 'YColor', 'w');
    plot(fk, eta, 'm-', 'LineWidth', 1.5);
    ylabel('Rendimento (%)');
    xlabel('Frequência (kHz)');
    ylim([0 100]);
    title('(d)', 'Color', 'w');
    grid on;
end

function plotar_3D(num, X, Y, Z1, Z2, nome, zlim1, zlim2)
    figure(num);
    set(gcf, 'Name', nome, 'Color', 'k', 'Position', [150, 150, 950, 420]);
    
    subplot(1,2,1);
    set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    surf(X*1e6, Y*1e6, Z1, 'EdgeColor', 'none');
    colormap jet; colorbar('Color', 'w');
    xlabel('L_s (\muH)'); ylabel('L_p (\muH)'); zlabel('V_p (V)');
    title('(a) Tensão de Entrada Necessária', 'Color', 'w');
    zlim(zlim1); view(-45, 30); grid on;
    
    subplot(1,2,2);
    set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    surf(X*1e6, Y*1e6, Z2, 'EdgeColor', 'none');
    colormap jet; colorbar('Color', 'w');
    xlabel('L_s (\muH)'); ylabel('L_p (\muH)'); zlabel('Rendimento (%)');
    title('(b) Rendimento Projetado', 'Color', 'w');
    zlim(zlim2); view(-45, 30); grid on;
end