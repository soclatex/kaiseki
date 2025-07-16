% ---- 問題4：熱方程式の陽的差分法（Forward Euler法） ----

% 空間方向の設定：x ∈ [0, 1] を 30 等分 → 31 区間、32 点
Nx = 30; % xについての分割数を変更              
dx = 1 / Nx;
x = linspace(0, 1, Nx+1);   % xの離散点（Nx+1個）

% 時間方向の設定：t ∈ [0, 1] を 3000 等分 → 3001 個の時刻点
Nt = 3000; %tについての分割数を変更
dt = 1 / Nt;
t = linspace(0, 1, Nt+1);   % tの離散点（Nt+1個）

% 安定条件（陽解法の場合）： dt/dx^2 ≦ 0.5 を満たす必要がある
alpha = dt / dx^2;
if alpha > 0.5
    warning('安定条件を満たしていません！ alpha = %.3f > 0.5', alpha);
end

% u(i,n)：x_i における時刻 t_n の値を格納する配列（初期化）
u = zeros(Nx+1, Nt+1);

% 初期条件 u(x,0) = f(x) = -7x(x-1)　<-このxの部分をx(i)として問題によって実装↓
for i = 1:Nx+1
    u(i,1) = -7 * x(i) * (x(i) - 1 ) ;  % i番目のxに対してf(x_i)を代入
end

% 境界条件：u(0,t) = u(1,t) = 0（明示的に毎時刻に設定）
u(1,:)   = 0;
u(end,:) = 0;

% 陽的差分法（forward Euler 法）
% u_i^{n+1} = u_i^n + alpha*(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
for n = 1:Nt
    for i = 2:Nx   % i=2〜Nx の範囲（端は境界条件）
        u(i,n+1) = u(i,n) + alpha * (u(i+1,n) - 2*u(i,n) + u(i-1,n));
    end
end

% x_11 = x(12番目) ⇒ MATLABではu(12,:)が対応
% t_3001 = 最後の時刻 ⇒ MATLABではNt+1番目 = u(:, end)
fprintf('u(x11, t3001) = %.4e\n', u(12, end));
