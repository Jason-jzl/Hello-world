%% run_PAH_project.m
% Main script for PAH endothelial signaling project
% 课程 project 的主脚本：三种场景 + 参数扫描 + 简单"治疗"模拟

clear; clc; close all;

%% 1. 三种基线场景参数

% (1) normal BMPR2 + normal shear
param_normal.BMPR2_level  = 1.0;
param_normal.shear_level  = 1.0;
param_normal.tspan        = [0 100];

% (2) low BMPR2 + normal shear
param_lowBMPR2            = param_normal;
param_lowBMPR2.BMPR2_level = 0.2;

% (3) abnormal (low) shear + normal BMPR2
param_lowShear            = param_normal;
param_lowShear.shear_level = 0.2;

%% 2. 仿真三种场景

out_normal   = simulate_PAH_model(param_normal);
out_lowBMPR2 = simulate_PAH_model(param_lowBMPR2);
out_lowShear = simulate_PAH_model(param_lowShear);

names = out_normal.speciesNames;
idxNO     = find(strcmpi(names,'NO'),1);
idxProlif = find(strcmpi(names,'PASMC_Proliferation'),1);

%% 3. 画时间曲线：NO 和 PASMC_Proliferation

figure;
subplot(2,1,1);
plot(out_normal.t,   out_normal.y(:,idxNO),   'LineWidth',1.5); hold on;
plot(out_lowBMPR2.t, out_lowBMPR2.y(:,idxNO), 'LineWidth',1.5);
plot(out_lowShear.t, out_lowShear.y(:,idxNO), 'LineWidth',1.5);
ylabel('NO');
legend({'normal','low BMPR2','low shear'},'Location','best');
title('Time courses of NO');

subplot(2,1,2);
plot(out_normal.t,   out_normal.y(:,idxProlif),   'LineWidth',1.5); hold on;
plot(out_lowBMPR2.t, out_lowBMPR2.y(:,idxProlif), 'LineWidth',1.5);
plot(out_lowShear.t, out_lowShear.y(:,idxProlif), 'LineWidth',1.5);
ylabel('PASMC proliferation');
xlabel('Time');
legend({'normal','low BMPR2','low shear'},'Location','best');
title('Time courses of PASMC proliferation');

%% 4. 稳态条形图（NO & Proliferation）

NO_vals   = [out_normal.NO_steady,   out_lowBMPR2.NO_steady,   out_lowShear.NO_steady];
PRO_vals  = [out_normal.INFL_steady, out_lowBMPR2.INFL_steady, out_lowShear.INFL_steady]; % 用 INFL_steady 存的 PASMC_Proliferation

labels = {'normal','low BMPR2','low shear'};

figure;
subplot(1,2,1);
bar(NO_vals);
set(gca,'XTickLabel',labels);
ylabel('NO (steady)');
title('NO steady state');

subplot(1,2,2);
bar(PRO_vals);
set(gca,'XTickLabel',labels);
ylabel('PASMC proliferation (steady)');
title('PASMC proliferation steady state');

sgtitle('Steady-state comparison across scenarios');

%% 5. 参数扫描：BMPR2_level × shear_level

BMPR2_vec = linspace(0.1, 1.0, 10);
shear_vec = linspace(0.1, 1.0, 10);

NO_matrix   = zeros(length(BMPR2_vec), length(shear_vec));
PRO_matrix  = zeros(size(NO_matrix));

for i = 1:length(BMPR2_vec)
    for j = 1:length(shear_vec)
        param_scan             = param_normal;
        param_scan.BMPR2_level = BMPR2_vec(i);
        param_scan.shear_level = shear_vec(j);

        out_scan = simulate_PAH_model(param_scan);

        NO_matrix(i,j)  = out_scan.NO_steady;
        PRO_matrix(i,j) = out_scan.INFL_steady;
    end
end

% NO heatmap
figure;
imagesc(shear_vec, BMPR2_vec, NO_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('shear level');
ylabel('BMPR2 level');
title('Steady-state NO');

% PASMC proliferation heatmap
figure;
imagesc(shear_vec, BMPR2_vec, PRO_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('shear level');
ylabel('BMPR2 level');
title('Steady-state PASMC proliferation');

%% 6. 简单"治疗"模拟：disease vs therapy

param_disease             = param_normal;
param_disease.BMPR2_level = 0.2;
param_disease.shear_level = 0.2;

param_therapy             = param_disease;
param_therapy.BMPR2_level = 1.0;
param_therapy.shear_level = 1.0;

out_dis   = simulate_PAH_model(param_disease);
out_thera = simulate_PAH_model(param_therapy);

figure;
bar( [out_dis.NO_steady,   out_thera.NO_steady;
      out_dis.INFL_steady, out_thera.INFL_steady] );
set(gca,'XTickLabel',{'NO','PASMC prolifer.'});
legend({'disease','therapy'},'Location','best');
ylabel('steady-state value');
title('Effect of simple in silico therapy');
