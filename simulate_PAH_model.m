function outputs = simulate_PAH_model(param)
% SIMULATE_PAH_MODEL
% Wrapper for a logic-based PAH endothelial model specified in Netflux Excel.
% 用 Netflux Excel 里那套节点/反应，自己在 MATLAB 里实现的连续逻辑 ODE 模型。
%
% INPUT (param struct, 所有量都 0–1):
%   .BMPR2_level  : 0–1, 1 = 正常 BMPR2, 0 = 完全失活/抑制
%   .shear_level  : 0–1, 1 = 规则层流, 0 = 完全振荡/异常剪切
%   .tspan        : [t0 tf], 仿真时间，不给就默认 [0 100]
%
% OUTPUT:
%   outputs.t            - 时间向量
%   outputs.y            - 状态矩阵 (time x species)
%   outputs.speciesNames - 物种名字（和 Excel 一致）
%   outputs.NO_steady    - NO 稳态
%   outputs.ET1_steady   - 这里暂时没有 ET-1 节点, 先返回 NaN
%   outputs.INFL_steady  - 这里用 PASMC_Proliferation 当"坏结局"，返回其稳态
%
% 物种顺序 (12 个，完全按照你 Excel 的 model):
%   1  LaminarFlow
%   2  OscillatoryFlow
%   3  KLF2
%   4  BMP4
%   5  SMAD1_5
%   6  BMPR2
%   7  BMPR2_mutation
%   8  BMPR2_inhibitors
%   9  PKA
%   10 eNOS_Phos
%   11 NO
%   12 PASMC_Proliferation

    % ------------ 默认参数 ------------
    if ~isfield(param, 'BMPR2_level')
        param.BMPR2_level = 1.0;   % 正常
    end
    if ~isfield(param, 'shear_level')
        param.shear_level = 1.0;   % 规则层流
    end
    if ~isfield(param, 'tspan')
        param.tspan = [0 100];
    end

    % ------------ 把高层参数映射到具体"输入节点" ------------
    % shear_level 高 → LaminarFlow 高, OscillatoryFlow 低
    p.Laminar_input = param.shear_level;          % => LaminarFlow
    p.Osc_input     = 1 - param.shear_level;      % => OscillatoryFlow

    % BMPR2_level 低 → mutation / inhibitor 高
    mut_in = 1 - param.BMPR2_level;
    p.Mut_input   = mut_in;   % => BMPR2_mutation
    p.Inhib_input = mut_in;   % => BMPR2_inhibitors

    % 时间
    tspan = param.tspan;

    % 初值：全部从 0 开始（对应 Excel 里的 Yinit=0）
    y0 = zeros(12,1);

    % ------------ ODE 仿真 ------------
    odeFun = @(t,y) pah_logic_ode(t, y, p);
    opts   = odeset('RelTol',1e-6,'AbsTol',1e-8);
    [t, y] = ode15s(odeFun, tspan, y0, opts);

    % 物种名字（和 Excel 对齐）
    speciesNames = { ...
        'LaminarFlow', ...
        'OscillatoryFlow', ...
        'KLF2', ...
        'BMP4', ...
        'SMAD1_5', ...
        'BMPR2', ...
        'BMPR2_mutation', ...
        'BMPR2_inhibitors', ...
        'PKA', ...
        'eNOS_Phos', ...
        'NO', ...
        'PASMC_Proliferation'};

    % 索引
    idxNO     = 11;
    idxProlif = 12;

    NO_steady   = y(end, idxNO);
    Prol_steady = y(end, idxProlif);

    % 这里暂时没有 ET-1 这个节点，用 NaN 占位，方便后面代码兼容
    ET1_steady = NaN;

    % ------------ 打包输出 ------------
    outputs.t            = t;
    outputs.y            = y;
    outputs.speciesNames = speciesNames;
    outputs.NO_steady    = NO_steady;
    outputs.ET1_steady   = ET1_steady;
    outputs.INFL_steady  = Prol_steady;   % 用 PASMC_Proliferation 作为"坏结局/炎症"代理
end


function dydt = pah_logic_ode(~, y, p)
% PAH_LOGIC_ODE
% 这里是根据你 Netflux Excel 里反应表，写成的连续逻辑 ODE：
%
%  Reactions (和 Excel 一一对应)：
%   r1  => LaminarFlow
%   r2  => OscillatoryFlow
%   r3  'LaminarFlow           => KLF2
%   r4  '!OscillatoryFlow      => KLF2
%   r5  '!KLF2                 => BMP4
%   r6  'OscillatoryFlow       => BMP4
%   r7  'BMP4                  => SMAD1_5
%   r8  '!LaminarFlow          => SMAD1_5
%   r9  => BMPR2_mutation
%   r10 => BMPR2_inhibitors
%   r11 '!BMPR2_mut & !BMPR2_inhib => BMPR2
%   r12 'BMPR2                => PKA
%   r13 'PKA                  => eNOS_Phos
%   r14 'eNOS_Phos            => NO
%   r15 '!NO                  => PASMC_Proliferation

    % 物种拆出来（方便写方程）
    Lam   = y(1);
    Osc   = y(2);
    KLF2  = y(3);
    BMP4  = y(4);
    SMAD  = y(5);
    BMPR2 = y(6);
    Mut   = y(7);
    Inhib = y(8);
    PKA   = y(9);
    eNOSP = y(10);
    NO    = y(11);
    Prol  = y(12);

    tau = 1.0;   % 所有节点统一 time constant，简单起见

    % -------- 输入节点按一阶动力学靠近"外源输入" --------
    Lam_in   = p.Laminar_input;
    Osc_in   = p.Osc_input;
    Mut_in   = p.Mut_input;
    Inhib_in = p.Inhib_input;

    dLam   = (Lam_in   - Lam)   / tau;
    dOsc   = (Osc_in   - Osc)   / tau;
    dMut   = (Mut_in   - Mut)   / tau;
    dInhib = (Inhib_in - Inhib) / tau;

    % -------- 中间节点采用简单的连续逻辑规则 --------
    % 1) KLF2 高：LaminarFlow 高 且 OscillatoryFlow 低
    %    用  f_KLF2 = Lam * (1 - Osc)
    f_KLF2 = Lam * (1 - Osc);
    dKLF2  = (f_KLF2 - KLF2) / tau;

    % 2) BMP4 高：Osc 高 或 KLF2 低
    %    用 "软 OR"：1 - (1 - Osc) * (KLF2)
    f_BMP4 = 1 - (1 - Osc) * (KLF2);
    dBMP4  = (f_BMP4 - BMP4) / tau;

    % 3) SMAD1_5 高：BMP4 高 或 Laminar 低
    %    同样用 1 - (1 - BMP4) * (Lam)
    f_SMAD = 1 - (1 - BMP4) * (Lam);
    dSMAD  = (f_SMAD - SMAD) / tau;

    % 4) BMPR2 高：mutation 低 且 inhibitor 低
    %    f_BMPR2 = (1 - Mut) * (1 - Inhib)
    f_BMPR2 = (1 - Mut) * (1 - Inhib);
    dBMPR2  = (f_BMPR2 - BMPR2) / tau;

    % 5) 下游 BMPR2 → PKA → eNOS_Phos → NO
    f_PKA   = BMPR2;
    dPKA    = (f_PKA   - PKA)   / tau;

    f_eNOSP = PKA;
    deNOSP  = (f_eNOSP - eNOSP) / tau;

    f_NO    = eNOSP;
    dNO     = (f_NO    - NO)    / tau;

    % 6) PASMC_Proliferation 高：NO 低
    f_Prol  = 1 - NO;
    dProl   = (f_Prol  - Prol)  / tau;

    dydt = [ dLam; dOsc; dKLF2; dBMP4; dSMAD; ...
             dBMPR2; dMut; dInhib; dPKA; deNOSP; dNO; dProl ];
end
