import numpy as np
from numpy.typing import ArrayLike
from dataclasses import dataclass
from typing import Tuple


@dataclass
class DiracParams:
    v: float = 1.0  # 费米速度
    Omega: float = 1.0  # 驱动频率 Ω = 2π/T
    A0: float = 0.2  # 矢势幅值（圆偏振）：A(t) = A0 (cosΩt, sinΩt)
    valley_tau_z: int = 1  # 谷指标 τ_z = ±1
    trunc_M: int = 3  # Floquet 截断 |m| ≤ M
    Nt: int = 2048  # 时间积分采样数（越大越精确）


def pauli_matrices() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return sigma_x, sigma_y, sigma_z


def vector_potential(t: ArrayLike, A0: float, Omega: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    圆偏振光的 Peierls 替换：k → k + A(t)，A(t) = A0 (cosΩt, sinΩt)。
    返回 A_x(t), A_y(t)。
    """
    t = np.asarray(t) # t是一维数组，记载了时间序列。长度由Nt决定。
    Ax = A0 * np.cos(Omega * t) # Ax同样为一维数组，记载了在不同时刻下Ax取值
    Ay = A0 * np.sin(Omega * t) # Ay同样为一维数组，记载了在不同时刻下Ay取值

    # t，Ax，Ay的维度都相同。它们都是一个含时序列。
    # 因此，在这个程序中不是给一个时间，再去生成当前时刻的A。
    # 而是把时间设成list，同时让各个变量成为时序大数组


    return Ax, Ay


def dirac_h_k_t(kx: float, ky: float, t: ArrayLike, params: DiracParams) -> np.ndarray:
    """
    构造时间依赖哈密顿量 H(k,t) = τ_z v [(kx+Ax) σ_x + (ky+Ay) σ_y]。
    返回形状 (Nt, 2, 2) 的数组，对应每个时间点的 2x2 矩阵。
    """
    sx, sy, _ = pauli_matrices() # sigma_z用哑变量接收，但丢掉不用。
    tau = params.valley_tau_z
    v = params.v
    Ax, Ay = vector_potential(t, params.A0, params.Omega) # Ax.shape=(Nt,)  Ay.shape=(Nt,)
    kx_t = kx + Ax # kx_t.shape=(Nt,)
    ky_t = ky + Ay # ky_t.shape=(Nt,)
    # 按广播生成
    Ht = tau * v * (kx_t[:, None, None] * sx[None, :, :] + ky_t[:, None, None] * sy[None, :, :]) 

    # tau: int
    # v: float
    # kx_t[:, None, None].shape)=（Nt，1，1），等价于kx_t.reshape(N, 1, 1)。
    # 这种写法常常用于广播 (broadcasting)，比如把一维向量扩展成三维张量，以便和 (N, M, K) 这样的数组在运算时自动匹配。
    # sx[None, :, :].shape = (1,2,2)
    # Ht.shape=(Nt,2,2)

    # Nt把矢势撑起来了，矢势把H（k，t）撑起来了

    return Ht


def fourier_block_Hmn(kx: float, ky: float, m: int, n: int, params: DiracParams) -> np.ndarray:
    """
    通过时间积分计算 H_{mn}(k) = (1/T) ∫_0^T H(k,t) e^{i(m-n)Ω t} dt。
    数值积分采用均匀采样：t_j = j Δt, Δt = T/Nt。
    返回 2x2 复矩阵。
    """
    Omega = params.Omega #驱动频率，type(Omega)=<class 'float'>
    T = 2.0 * np.pi / Omega #时间周期, type(T)=<class 'float'>
    Nt = params.Nt #时间序列个数（时间积分采样数）,type(Nt)=<class 'int'>
    t = np.linspace(0.0, T, Nt, endpoint=False) #从0开始，到T结束，生成一个长度为Nt的时间序列.  t.shape=(Nt,)
    Ht = dirac_h_k_t(kx, ky, t, params)  # Ht.shape=(Nt,2,2)
    phase = np.exp(1.0j * (m - n) * Omega * t)  # phase.shape=(Nt,)

    # 复合梯形/矩形求和： (1/T) * Σ H(t_j) e^{i(m-n)Ω t_j} Δt
    dt = T / Nt # 小区间长度
    integrand = Ht * phase[:, None, None] # phase[:, None, None].shape=(Nt,1,1). Ht.shape=(Nt,2,2)

    # integrand 是一个 长度为 Nt 的列表（沿时间轴），每个元素是一个 2×2 矩阵
    # numpy.sum 会在指定的轴上做逐元素求和. 结果就是 (2,2) 的矩阵，每个元素等于该元素在时间维度上所有切片的和
    # 先求和,再乘以归一化系数。

    # 思考：如果积分当中引入权重函数，使得积分过程不平权，该怎么做？

    Hmn = integrand.sum(axis=0) * (dt / T)
    return Hmn


def build_floquet_sambe(kx: float, ky: float, params: DiracParams) -> np.ndarray:
    """
    使用 H_{mn} 的时间积分结果，构造 Sambe 空间 Floquet 哈密顿量。
    采用等式 ∑_n H_{mn}|u^n⟩ = (ε + mΩ)|u^m⟩ 的实现方式，
    即在对角块处加上 mΩ I_2（与现有脚本保持一致）。
    仅线性耦合的圆偏振会产生 0 与 ±1 次谐波，
    但这里仍按通式通过积分生成，以确保形式正确且可拓展。
    """
    M = params.trunc_M # 截断阶数
    dim_block = 2 # 每个 Floquet block 对应原始系统的 Hilbert 空间维度（Dirac 是 2x2）。
    dim_total = dim_block * (2 * M + 1) # H_F大矩阵的维度
    HF = np.zeros((dim_total, dim_total), dtype=complex)

    # 这个函数的作用：给定Floquet指标m，返回它在大矩阵中的切片范围（idx，idx+dim_block）
    # 在两个方向上都要用到这个函数的
    def sl(m: int) -> slice: 
        idx = (m + M) * dim_block
        return slice(idx, idx + dim_block)

    # 计算所有 H_{mn}，避免重复计算。否则每调用一次就得重新算一下。
    ms = range(-M, M + 1)
    Hmn_cache = {}
    for m in ms:
        for n in ms:
            Hmn_cache[(m, n)] = fourier_block_Hmn(kx, ky, m, n, params)

    # 装配 Sambe 矩阵，并在对角加上 mΩ I_2。注意符号的convention
    Omega = params.Omega
    for m in ms:
        for n in ms:
            HF[sl(m), sl(n)] = Hmn_cache[(m, n)]
        HF[sl(m), sl(m)] += (m * Omega) * np.eye(dim_block)

    return HF


# 对角化H_F，得到本征矢，及其对应的准能。
def solve_floquet(kx: float, ky: float, params: DiracParams) -> Tuple[np.ndarray, np.ndarray]:
    HF = build_floquet_sambe(kx, ky, params)
    evals, evecs = np.linalg.eigh(HF)
    # evals.shape = (2*trunc_M +1 ,)
    # evecs.shape = (2*trunc_M +1 , 2*trunc_M +1)
    return evals, evecs


# 计算静态权重。实际上这个函数能计算任意阶的权重。
# 在给定k点，有2*trunc_M +1个能量本征态。这些能量本征态用一个dim_block的列向量描述。

# Floquet阶： -m -m+1 -m+2 -m+3 -m+4  ……   0    1    ……   m
# Python索引： 0   1    2    3    4   ……   m   m+1   ……   2m

def static_component_weights(evecs: np.ndarray, M: int) -> np.ndarray:
    dim_block = 2
    start = (0 + M) * dim_block
    sl0 = slice(start, start + dim_block)
    comp = evecs[sl0, :] # comp.shape=(dim_block, 2*trunc_M +1)。取出0阶（或m阶）子空间的本征矢部分
    w = np.sum(np.abs(comp) ** 2, axis=0).real
    return np.clip(w, 0.0, 1.0) # 保证结果在 [0, 1] 区间内


def bands_along_path(ks: ArrayLike, params: DiracParams) -> Tuple[np.ndarray, np.ndarray]:
    ks = np.asarray(ks, dtype=float)
    dim = 2 * (2 * params.trunc_M + 1)
    energies = np.zeros((ks.shape[0], dim))
    weights = np.zeros((ks.shape[0], dim))
    for i, (kx, ky) in enumerate(ks):

        e, vecs = solve_floquet(kx, ky, params)

        # Floquet阶： -m -m+1 -m+2 -m+3 -m+4  ……   0    1    ……   m
        # Python索引： 0   1    2    3    4   ……   m   m+1   ……   2m
        w = static_component_weights(vecs, params.trunc_M)

        # 返回能量从小到大排序的索引。
        # 然后同时对能量 e 和权重 w 按照相同顺序排序，这样保证能量和对应的 Floquet 阶投影权重对齐
        idx = np.argsort(e)
        energies[i] = e[idx]
        weights[i] = w[idx]

        #print(f"Calculating:{i+1}/{len(ks)}，{(i+1)/len(ks)*100:.2f}%")

    return energies, weights


def demo():
    import importlib

    params = DiracParams(v=1.0, Omega=1.0, A0=0.3, valley_tau_z=1, trunc_M=4, Nt=4096)

    kmax = 2.0
    num_k = 181
    kxs = np.linspace(-kmax, kmax, num_k)
    ks = np.stack([kxs, np.zeros_like(kxs)], axis=1)
    bands, weights = bands_along_path(ks, params)

    np.save("floquet_bands.npy", bands)
    np.save("floquet_weights_m0.npy", weights)

    mpl = importlib.util.find_spec("matplotlib")
    if mpl is None:
        print("[info] 已保存数据到 floquet_bands.npy；若需绘图请安装 matplotlib。")
        return

    import matplotlib.pyplot as plt

    plt.figure(figsize=(7, 4))
    kk = np.repeat(kxs, bands.shape[1])
    ee = bands.reshape(-1)
    ww = weights.reshape(-1)
    sc = plt.scatter(kk, ee, c=ww, s=10, cmap="viridis", vmin=0.0, vmax=1.0)
    cbar = plt.colorbar(sc)
    cbar.set_label("weight of m=0 component")
    plt.xlabel("$k_x$")
    plt.ylabel("quasi-energy  ε  (unfolded)")
    plt.title(
        rf"2D massless Dirac Floquet bands (Hmn by time integration), M={params.trunc_M}, A0={params.A0}, Ω={params.Omega}, τ_z={params.valley_tau_z}"
    )
    plt.tight_layout()
    plt.savefig("floquet_bands.png", dpi=160)
    plt.show()


if __name__ == "__main__":
    demo()


