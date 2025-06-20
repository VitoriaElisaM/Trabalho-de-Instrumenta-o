%---------------------------------------------------------%
% Calibração Estática do Sensor LDR com Múltiplas Medições
% Ajuste Linear: Lux = f(Tensão Média do LDR)
%---------------------------------------------------------%

clear; clc; close all;

%---------------------------------------------------------%
% Parâmetros do sistema
%---------------------------------------------------------%
V_ref = 3.3;         % Tensão de alimentação do divisor de tensão (V)
ADC_max = 1023;      % Resolução do ADC (10 bits)

%---------------------------------------------------------%
% Carregamento dos dados
%---------------------------------------------------------%
dados = readtable('ldr.csv', 'VariableNamingRule', 'preserve');

% Verificação de valores ausentes
if any(ismissing(dados), 'all')
    error('Há valores ausentes na tabela ldr.csv. Verifique os dados.');
end

% Leitura das colunas
adc1 = dados{:,1};
adc2 = dados{:,2};
adc3 = dados{:,3};
lux  = dados{:,4};  % valor de referência

% Conversão das 3 medidas para tensão
tensao1 = V_ref * adc1 / ADC_max;
tensao2 = V_ref * adc2 / ADC_max;
tensao3 = V_ref * adc3 / ADC_max;

% Média das tensões
tensao_media = mean([tensao1, tensao2, tensao3], 2);

%---------------------------------------------------------%
% Ajuste linear: lux = a·tensao + b
%---------------------------------------------------------%
p = polyfit(tensao_media, lux, 1);
lux_estimado = polyval(p, tensao_media);

% Cálculo do coeficiente de determinação (R²)
lux_media = mean(lux);
SST = sum((lux - lux_media).^2);     % Soma total dos quadrados
SSE = sum((lux - lux_estimado).^2);  % Soma dos quadrados dos erros
R2 = 1 - SSE/SST;

% Cálculo dos resíduos
residuos = lux - lux_estimado;

%---------------------------------------------------------%
% Cálculo da incerteza dos coeficientes
%---------------------------------------------------------%
n = length(tensao_media);
X = [tensao_media ones(n,1)];
res_var = var(residuos);
cov_p = res_var * inv(X' * X);
inc_p = sqrt(diag(cov_p));  % Incertezas dos coeficientes

%---------------------------------------------------------%
% Plotagem dos resultados
%---------------------------------------------------------%
figure('Name', 'Calibração do LDR');

% Gráfico do ajuste
subplot(2,1,1)
plot(tensao_media, lux, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
hold on
plot(tensao_media, lux_estimado, 'r-', 'LineWidth', 2)
grid on
xlabel('Tensão média do LDR (V)')
ylabel('Iluminação (lux)')
title(sprintf('Ajuste Linear'))
legend('Dados experimentais', 'Ajuste linear', 'Location', 'best')

% Gráfico dos resíduos
subplot(2,1,2)
plot(tensao_media, residuos, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
hold on
yline(0, 'r--')
grid on
xlabel('Tensão média do LDR (V)')
ylabel('Resíduo (lux)')
title(sprintf('Resíduos'))

%---------------------------------------------------------%
% Exibição no console
%---------------------------------------------------------%
fprintf('\nModelo Ajustado:\n');
fprintf('lux = %.4f·V + %.4f\n', p(1), p(2));
fprintf('Incerteza da inclinação: ±%.4f\n', inc_p(1));
fprintf('Incerteza do offset: ±%.4f\n', inc_p(2));
fprintf('R² = %.4f\n', R2);
fprintf('Desvio padrão dos resíduos: %.4f lux\n', std(residuos));

%---------------------------------------------------------%
% Exportação opcional
%---------------------------------------------------------%
T = table(tensao_media, lux, lux_estimado, residuos);
writetable(T, 'resultados_calibracao_ldr.csv');

%%

%---------------------------------------------------------%
% Calibração Estática do Sensor Aplicativo de Smartphone com Múltiplas Medições
%---------------------------------------------------------%

clc; close all;

%---------------------------------------------------------%
% Carregamento dos dados
%---------------------------------------------------------%
dados = readtable('app.csv', 'VariableNamingRule', 'preserve');

% Verificação de valores ausentes
if any(ismissing(dados), 'all')
    error('Há valores ausentes na tabela ldr.csv. Verifique os dados.');
end

% Leitura das colunas
lx1 = dados{:,1};
lx2 = dados{:,2};
lx3 = dados{:,3};
lux  = dados{:,4};  % valor de referência

% Média das tensões
lx_media = mean([lx1, lx2, lx3], 2);

%---------------------------------------------------------%
% Ajuste linear: lux = a·tensao + b
%---------------------------------------------------------%
p = polyfit(lx_media, lux, 1);
lux_estimado = polyval(p, lx_media);

% Cálculo do coeficiente de determinação (R²)
lux_media = mean(lux);
SST = sum((lux - lux_media).^2);     % Soma total dos quadrados
SSE = sum((lux - lux_estimado).^2);  % Soma dos quadrados dos erros
R2 = 1 - SSE/SST;

% Cálculo dos resíduos
residuos = lux - lux_estimado;

%---------------------------------------------------------%
% Cálculo da incerteza dos coeficientes
%---------------------------------------------------------%
n = length(lx_media);
X = [lx_media ones(n,1)];
res_var = var(residuos);
cov_p = res_var * inv(X' * X);
inc_p = sqrt(diag(cov_p));  % Incertezas dos coeficientes

%---------------------------------------------------------%
% Plotagem dos resultados
%---------------------------------------------------------%
figure('Name', 'Calibração do Aplicativo');

% Gráfico do ajuste
subplot(2,1,1)
plot(lx_media, lux, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
hold on
plot(lx_media, lux_estimado, 'r-', 'LineWidth', 2)
grid on
xlabel('Leitura do aplicativo (lux)')
ylabel('Iluminação (lux)')
title(sprintf('Ajuste Linear'))
legend('Dados experimentais', 'Ajuste linear', 'Location', 'best')

% Gráfico dos resíduos
subplot(2,1,2)
plot(lx_media, residuos, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
hold on
yline(0, 'r--')
grid on
xlabel('Leitura média do aplicativo (lux)')
ylabel('Resíduo (lux)')
title(sprintf('Resíduos'))

%---------------------------------------------------------%
% Exibição no console
%---------------------------------------------------------%
fprintf('\nModelo Ajustado:\n');
fprintf('lux = %.4f·V + %.4f\n', p(1), p(2));
fprintf('Incerteza da inclinação: ±%.4f\n', inc_p(1));
fprintf('Incerteza do offset: ±%.4f\n', inc_p(2));
fprintf('R² = %.4f\n', R2);
fprintf('Desvio padrão dos resíduos: %.4f lux\n', std(residuos));

%---------------------------------------------------------%
% Exportação opcional
%---------------------------------------------------------%
T = table(lx_media, lux, lux_estimado, residuos);
writetable(T, 'resultados_calibracao_ldr.csv');

%% Análise temporal

%---------------------------------------------------------%
% Análise de variação de iluminação com LDR calibrado e app
%---------------------------------------------------------%

close all; clc;

%---------------------------------------------------------%
% Parâmetros do sistema
%---------------------------------------------------------%
V_ref = 3.3;         % Tensão de alimentação do divisor de tensão (V)
ADC_max = 1023;      % Resolução do ADC (10 bits)

% Coeficientes da calibração do LDR (substitua pelos seus)
a = p(1);          % Exemplo: lux = a·V + b
b = p(2);

%---------------------------------------------------------%
% Carregamento dos dados do LDR
%---------------------------------------------------------%
dados_ldr = readtable('ldr_temporal.csv');  % Nome do arquivo com os dados contínuos
adc_ldr = dados_ldr.("adc_ldr");

% Número de amostras por grupo (15 amostras = 15 s)
N = 15;
num_blocos = floor(length(adc_ldr) / N);

adc_media = zeros(num_blocos, 1);
for k = 1:num_blocos
    idx_ini = (k-1)*N + 1;
    idx_fim = k*N;
    adc_media(k) = mean(adc_ldr(idx_ini:idx_fim));
end

% Conversão para tensão média
tensao_media = V_ref * adc_media / ADC_max;

% Estimativa de lux com base no modelo calibrado
lux_estimado = a * tensao_media + b;

tempo_15s = (0:num_blocos-1)' * 15;  % vetor de tempo (em segundos)

%---------------------------------------------------------%
% Carregamento dos dados do aplicativo
%---------------------------------------------------------%
lux_app = readmatrix('app_temporal.csv');  % Valores anotados manualmente

% Verificação de consistência
if length(lux_app) ~= length(lux_estimado)
    warning('O número de pontos do app e do LDR não coincide.');
end

%---------------------------------------------------------%
% Gráfico comparativo
%---------------------------------------------------------%
figure('Name', 'Comparação Temporal: LDR calibrado vs App');

plot(tempo_15s, lux_estimado, '-ob', 'LineWidth', 2, 'DisplayName', 'LDR (estimado)');
hold on;
plot(tempo_15s, lux_app, '-sr', 'LineWidth', 2, 'DisplayName', 'App (medido)');
xlim([0 450])
grid on

xlabel('Tempo (s)')
ylabel('Iluminação (lux)')
title('Iluminação estimada pelo LDR calibrado vs App')
legend('Location', 'best')

%---------------------------------------------------------%
% Exportação opcional
%---------------------------------------------------------%
T = table(tempo_15s, lux_estimado, lux_app);
writetable(T, 'comparacao_lux_estimado_app.csv');

%% Implementação do Filtro de Kalman

%---------------------------------------------------------%
% Fusão de sensores com Filtro de Kalman
%---------------------------------------------------------%

close all; clear; clc;

%---------------------------------------------------------%
% Carregamento dos dados
%---------------------------------------------------------%
T = readtable('comparacao_lux_estimado_app.csv');
z1 = T.lux_estimado;  % Sensor LDR
z2 = T.lux_app;       % Sensor App
N = length(z1);

%---------------------------------------------------------%
% Parâmetros do modelo
%---------------------------------------------------------%
Q = 50;        % Variância do processo (ajustável)
R1 = 912.99;      % Variância da medição do LDR
R2 = 371.37;       % Variância da medição do App

%---------------------------------------------------------%
% Inicialização
%---------------------------------------------------------%
x_est = zeros(N, 1);        % Estimativa do estado
P = zeros(N, 1);            % Covariância do erro
K1 = zeros(N, 1);           % Ganho do Kalman para o sensor 1
K2 = zeros(N, 1);           % Ganho do Kalman para o sensor 2

x_est(1) = (z1(1) + z2(1)) / 2;  % Estimativa inicial: média dos sensores
P(1) = 1;                        % Incerteza inicial

%---------------------------------------------------------%
% Loop do Filtro de Kalman
%---------------------------------------------------------%
for k = 2:N
    % Predição
    x_pred = x_est(k-1);          % x_k|k-1 = x_{k-1}
    P_pred = P(k-1) + Q;          % P_k|k-1 = P_{k-1} + Q

    % Inovação para cada sensor
    y1 = z1(k) - x_pred;
    y2 = z2(k) - x_pred;

    % Ganhos de Kalman (dois sensores)
    K1(k) = P_pred / (P_pred + R1);
    K2(k) = P_pred / (P_pred + R2);

    % Fusão com média ponderada das inovações
    x_update = x_pred + (K1(k)*y1 + K2(k)*y2)/2;
    P_update = (1 - (K1(k) + K2(k))/2) * P_pred;

    % Atualização
    x_est(k) = x_update;
    P(k) = P_update;
end

tempo = T.tempo_15s;

figure('Name','Fusão com Filtro de Kalman');
plot(tempo, z1, '--ob', 'DisplayName', 'LDR (estimado)');
hold on;
plot(tempo, z2, '--sr', 'DisplayName', 'App (medido)');
plot(tempo, x_est, '-k', 'LineWidth', 2, 'DisplayName', 'Kalman (fusão)');
grid on;
xlabel('Tempo (s)');
ylabel('Iluminação (lux)');
legend('Location', 'best');
title('Fusão dos Sensores com Filtro de Kalman');
xlim([0 450])

figure('Name', 'Resíduos');
plot(tempo, x_est - z1, '-ob', 'DisplayName', 'Resíduo Kalman - LDR');
hold on;
plot(tempo, x_est - z2, '-sr', 'DisplayName', 'Resíduo Kalman - App');
yline(0, '--k');
xlabel('Tempo (s)');
ylabel('Resíduo (lux)');
legend('Location', 'best');
title('Resíduos do Filtro de Kalman');
xlim([0 450])
grid on;

std_ldr = std(z1);
std_app = std(z2);
std_kalman = std(x_est);

fprintf('Desvio padrão do LDR estimado: %.2f lux\n', std_ldr);
fprintf('Desvio padrão do App: %.2f lux\n', std_app);
fprintf('Desvio padrão do Kalman (fusão): %.2f lux\n', std_kalman);
