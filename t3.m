% ---- 問題3：連立常微分方程式のオイラー法による数値解 ----

% 区間 [0, 1] を 200 分割 ⇒ ステップ幅 h = 1/200
N = 200;           % ステップ数（分割数）ここを変える
h = 1 / N;         % ステップ幅
t = 0:h:1;         % 時間配列（確認用）

% 初期条件（x(0) = 0, y(0) = 1）
x = zeros(1, N+1); % x の値を格納
y = zeros(1, N+1); % y の値を格納
x(1) = 0;          % 初期条件設定(もしかしたら問題によって変える)
y(1) = 1;          % ここも 

% オイラー法のメインループ
for n = 1:N
    dxdt = x(n)^2 + 2*y(n);         % dx/dt の定義（問題によって変える）
    dydt = y(n) - 2*cos(x(n));      % dy/dt の定義（問題によって変える）

    % 次のステップの値を計算
    x(n+1) = x(n) + h * dxdt;
    y(n+1) = y(n) + h * dydt;
end

% x(1) を表示（指数表記で、小数点以下4桁）
fprintf('x(1) = %.4e\n', x(end));

% ---- グラフ描画 ----

% 1. x(t), y(t) の時間変化
figure;
plot(t, x, '-b', 'LineWidth', 1.5); hold on;
plot(t, y, '-r', 'LineWidth', 1.5);
xlabel('t'); ylabel('値');
legend('x(t)', 'y(t)');
title('オイラー法による x(t), y(t) の近似解');
grid on;

% 2. x-y 平面上の軌跡（位相図的なもの）
figure;
plot(x, y, '-k', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('x-y 平面上の軌跡');
grid on;