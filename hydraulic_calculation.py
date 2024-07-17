import numpy as np


# 定义colebrook函数，简易计算摩擦系数
def colebrook(Re, K):
    
    if Re < 2300:
        return 64 / Re
    else:
        # Colebrook方程的迭代求解
        def f(friction_factor):
            return -2.0 * np.log10((K / 3.7) + (2.51 / (Re * np.sqrt(friction_factor)))) - 1.0 / np.sqrt(friction_factor)
        
        # 初始猜测值
        friction_factor = 0.02
        for _ in range(10):
            friction_factor = friction_factor - f(friction_factor) / (2.0 / (friction_factor * np.log(10) * (K / 3.7 + 2.51 / (Re * np.sqrt(friction_factor)))))
        
        return friction_factor


# 定义管网参数
N = 3  # 热网节点总数
Np = 3  # 管道的数量
Nd = 2  # 负荷节点的数量
D = np.ones(3) * 150e-3  # 管道直径
ep = np.ones(3) * 1.25e-3  # 管道粗糙度
s = np.pi * D * D / 4  # 管道横截面积
len = np.array([400, 400, 600])  # 管道长度
mq = np.array([2, 3])  # 节点质量流量 kg/s
rho = 958.4  # 水在100度时的密度 kg/m^3
g = 9.81  # 重力加速度
viscosity = 0.294e-6  # 温度为100度时的动力粘度 m2/s

# 热网关联矩阵
A = np.array([[1, -1, 0], [0, 1, 1], [-1, 0, -1]])  # 节点-管道关联矩阵
B = np.array([[1, 1, -1]])  # 环路关联矩阵
dm = np.ones(Np)  # 初始化dm
err = 1
pre = 0

while err > 1e-3:
    # 计算管道流量dm
    m_node = np.dot(A, dm)  # 节点的流量注入
    dPhi = m_node[:N - 1] - mq  # 节点流量偏差值
    HJ0 = A[:N - 1, :]

    vel = dm / s / rho  # 单位 m kg/s, V m/s
    Re = abs(vel) * D / viscosity
    factor = np.zeros(Np)
    for ire in range(Np):
        if Re[ire] < 2300:
            factor[ire] = 64 / Re[ire]
        else:
            factor[ire] = colebrook(Re[ire], ep[ire] / D[ire])
    Kf = factor * len / D / s ** 2 / 2 / g / rho ** 2

    dpre = np.dot(B, Kf * abs(dm) * dm)  # 压力环方程
    HJpre = 2 * B * (Kf * abs(dm))

    dH = np.concatenate((dPhi, dpre))
    HJ = np.concatenate((HJ0, HJpre))

    dx = -np.linalg.solve(HJ, dH)
    err = max(abs(dH))
    dm += dx
    pre += 1

print("支路1流量：{:.3f}\n支路2流量：{:.3f}\n支路3流量：{:.3f}".format(dm[0], dm[1], dm[2]))
