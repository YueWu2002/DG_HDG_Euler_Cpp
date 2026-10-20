# 1D / 2D Euler DG 与凝聚 HDG 的数学定义

1D 和 2D 实现使用 Eigen，不链接 MFEM。数值定义参照
[MFEM 项目](https://github.com/YueWu2002/mfem_euler2d_cpp) 的一维简单波限制。
后文 2D 小节定义新增二维离散；既有 1D 数值协议保持不变。

## 方程与周期简单波

守恒变量 $U=(\rho,\rho u,E)^T$，$\gamma=3$ 为算例常量，
$p=(\gamma-1)(E-\rho u^2/2)$，
$F=(\rho u,\rho u^2+p,u(E+p))^T$。
区间为 $[-1,1]$，周期边界。令
$a=(3-\gamma)/2$，$b=(\gamma+1)/4$，通过
$x=\xi+t(a+b\sin(\pi\xi))$ 求足点，再令
$J_-=\sin(\pi\xi)$，$u=(2+J_-)/2$，
$c=(\gamma-1)(2-J_-)/4$，$\rho=(c^2/\gamma)^{1/(\gamma-1)}$，
$p=\rho^\gamma$。光滑区间为 $t<4/((\gamma+1)\pi)$。
初值采用高阶局部 L2 投影，必要时围绕单元平均值作保守正性缩放。

## Gauss 配点 DG

$n=k+1$ 个 Gauss-Legendre 节点，基函数 $\ell_i$、权重 $w_i$，
$M_{ei,ei}=h w_i/2$。主路径体积分使用 n 点 Gauss，质量对角；
`--quad Q` 支持 Q>=n 的诊断超积分，不能混用于默认性能比较。

时间阶段基值为 $U_b$，阶段步长 $\tau=a_{ss}\Delta t$。
实现存储物理弱残差

$$r_i=M_i(U_i-U_{b,i})/\tau-\sum_q w_q\ell'_i(x_q)F(U_q)
 +\ell_i(1)\widehat F_R-\ell_i(-1)\widehat F_L.$$

每个面冻结非负数 $a_f,b_f$（且 $a_f+b_f>0$）（Einfeldt/Roe HLL，安全系数 1.1）：

$$\widehat F={b_fF(U_L)+a_fF(U_R)+a_fb_f(U_L-U_R)\over a_f+b_f}.$$

CFL 速率为 $\max_e\alpha_e/h$，$\alpha_e$ 包括体积谱半径与两侧面 HLL 界；
$\Delta t=\mathrm{CFL}/\mathrm{rate}$，不额外乘次数因子，最后一步截断到 T。
每步开始冻结波速以选取步长；DIRK 后续阶段在前一阶段解处重新冻结，
阶段内 Newton 和线搜索不改变波速，通量导数随当前 Newton 状态更新。

## SSP-DIRK(4,3)

令 $\alpha=0.24169426078821$，$\beta=\alpha/4$，
$\eta=0.12915286960590$。系数与参考项目相同：

$$A=\begin{pmatrix}
\alpha&0&0&0\\-\alpha&\alpha&0&0\\0&1-\alpha&\alpha&0\\
\beta&\eta&1/2-\beta-\eta-\alpha&\alpha
\end{pmatrix},\qquad b=(0,1/6,1/6,2/3).$$

阶段基值 $U_b=U^n+\Delta t\sum_{j<s}a_{sj}K_j$，
求解 $r(U_s;U_b,a_{ss}\Delta t)=0$，
$K_s=(U_s-U_b)/(a_{ss}\Delta t)$，最终
$U^{n+1}=U^n+\Delta t\sum_s b_sK_s$。
`--integrator be` 为单阶段诊断和受控速度测试。

## HDG 与 DG Newton 更新等价

正向面两侧通量为
$F_L^H=F(U_L)+a_f(U_L-\lambda)$、
$F_R^H=F(U_R)+b_f(\lambda-U_R)$，骨架方程 $g=F_L^H-F_R^H=0$。
因此每个 Newton 状态重构

$$\lambda=\Phi(U)={a_fU_L+b_fU_R+F(U_L)-F(U_R)\over a_f+b_f}.$$

在该流形上 $r_H=r_{DG}$、$g=0$；线搜索和终止判据统一使用 DG 残差。
增广线性系统为 $[A\ B;C\ D](\delta U,\delta\lambda)^T=-(r,g)^T$，
$D_f=-(a_f+b_f)I$。按单元选主元 LU：

$$X=A^{-1}B,\quad y=A^{-1}r,\quad S=D-CX,$$
$$S\delta\lambda=-g+Cy,\qquad\delta U=-y-X\delta\lambda.$$

只通过回代得到 X 和 y，不显式存储逆。
缓存单元面块 $-CX$，Schur 乘法仅 gather、局部小块乘法、scatter。
Newton 下一次更新前重新投影 trace，不累积独立非线性 trace 误差。

## 与参考项目一致的残差测量

弱体积向量使用 $\|r\|_{M^{-1}}=\sqrt{r^TM^{-1}r}$。
物理 RMS 为 $\|r\|_{M^{-1}}/\sqrt{|\Omega|}$，这里 $|\Omega|=2$。
参考实现存储质量逆作用后的、乘阶段步长的方程；除去该行缩放后此定义相同。
阶段初始 RMS 记为 $r_0$，终止条件为
$\mathrm{RMS}\le\mathrm{ntol}\max(1,r_0)$。

GMRES 为固定右预条件 $JP^{-1}z=-r$，$\delta U=P^{-1}z$。
DG 的 Arnoldi 内积为 $M^{-1}$ 加权，1D 面为零维、trace 内积为欧氏。
两遍重正交、Givens 旋转、重启；每次重启后显式检查原方程真残差，
$\|b-Jx\|_W/\|b\|_W\le\mathrm{ltol}$。
不能用左预条件残差代替此数值。

HDG 还必须通过共同 DG 方程验收。令 Schur 误差
$e=S\delta\lambda-\mathrm{rhs}$，恢复体积残差为 $-BD^{-1}e$。
要求其 $M^{-1}$ 范数不超过
$\max(10\,\mathrm{ltol}\|r\|_{M^{-1}},\epsilon_{round})$，
$\epsilon_{round}=128\epsilon(\|U\|_M+\|U_b\|_M+
\tau\|r\|_{M^{-1}})/\tau$。
必要时将 Schur 容差乘 0.01 重解，最多四次；不降低 DG 验收标准。
Armijo 系数 $10^{-4}$，最多 30 次折半，带同一浮点舍入余量；
候选状态必须在节点、求积点与面保持密度及压力正性。

## 本次实验协议

收敛：DG/HDG，k=2，32/64/128/256/512，SSP-DIRK3，CFL=2，T=0.2，
安全系数 1.1，Newton/线性相对容差均为 1e-12，阶段至多 30 次更新。
比较全部三个分量 L2 误差、网格二分阶、DG/HDG 差及守恒。

速度：512 单元，k=0..8，BE，CFL=5，单个时间步固定 4 次 Newton 更新，
第 1 次预热，后 3 次平均；matrix-free、右 block-Jacobi、容差均为 1e-12。
固定更新工作负载不等同于完全收敛的计算，不能拿它报告终止时间精度。
加速比为同机串行 DG/HDG 的后 3 次平均耗时比，k=8 的 8 倍是目标，须实测。

## 2D 方程、旋转简单波和张量结构

$U=(\rho,\rho u,\rho v,E)$，$p=(\gamma-1)(E-\rho(u^2+v^2)/2)$。
令 $s=x\cos\alpha+y\sin\alpha$，将一维简单波的密度、压力、速度模长作用于 s，
速度方向取 $(\cos\alpha,\sin\alpha)$。计算域两个边长为
$L_x=2/\cos\alpha$、$L_y=2/\sin\alpha$，中心位于原点，分别周期。
此时冲击时间仍为 $4/((\gamma+1)\pi)$，不额外乘 $\sqrt2$。

参考单元是 $[-1,1]^2$，节点和体积点均为 n×n Gauss，n=k+1。
$M_{ij}=h_xh_yw_iw_j/4$。令 $D^w=WD$，体积弱项分别为
$-(h_y/2)(D^w)^T F_x W$ 和 $-(h_x/2)W F_y D^w$。
按两个一维收缩实现 O(n³) 单元作用，不构造 n²×n² 插值或微分矩阵。
CFL 速率为 $\max_e(\alpha_x/h_x+\alpha_y/h_y)$。
Newton 物理 RMS 为 $\|r\|_{M^{-1}}/\sqrt{L_xL_y}$；
trace 的真残差采用物理面质量逆权重。DIRK3、冻结策略、相对阈值与 1D 相同。

### 压缩 Jacobian 缓存

体积和每侧面的物理 Jacobian 仅由速度 u、v、总焓 H 和 gamma 决定。
对方向变量 $(d\rho,dm_x,dm_y,dE)$，
$dp=(\gamma-1)[dE-u dm_x-v dm_y+(u²+v²)d\rho/2]$，
$dF_x=(dm_x,2u dm_x-u²d\rho+dp,u dm_y+v dm_x-uv d\rho,
H(dm_x-u d\rho)+u(dE+dp))$，y 方向对称。
缓存三场而不是两个完整 4×4 Jacobian；这是解析导数的精确表示。

### 2D 精确 HDG 凝聚

每单元 m=4n²，四面合计 t=16n。矩阵 A 是 m×m，X 是 m×t。
使用一套 A/LU/右端工作区，一次求解 t+1 个右端得到 X 和 y，
仅按单元保留 X、y、局部 -CX。骨架 D 对角为 $-w_f(a_f+b_f)$。
令面限制矩阵为 $E_s$、$T_s=w_f(\mathrm{sign}_s F'_s+\tau_sI)$，
则 $C_s=T_sE_s$、$B_s=-w_f\tau_s E_s^T$。
对 X 先沿法向限制，再做面上 4×4 小矩阵作用，避免通用稠密 C×X。
恢复后的 DG 残差仍为 $-BD^{-1}e_S$，复用 GMRES 已形成的 Schur 真残差。

### 2D 精确块 Jacobi

配点下精确的质量归一化对角块为
$\widetilde A=I+\tau M^{-1}K_x+\tau M^{-1}K_y$，其中 K 含相应方向的体积项与对角面项。
Jacobi 单元块保留两个方向的完整耦合，直接对 $A=M\widetilde A/\tau$ 作带行选主元的 LU。
HDG 的面块同样作精确 LU。所有因子在每次 Newton 刷新后固定，供右预条件 GMRES 使用。

体积和面的张量收缩与原始离散代数等价，不替代完整块求解。
每单元 m=4n²，精确 LU 需 m²=16n⁴ 个 double。
实现仅全局保存一份 LU 和置换，一套装填矩阵和分解工作区在所有单元间复用；
不按单元再保存原始对角块或显式逆矩阵。
