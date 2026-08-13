import matplotlib.pyplot as plt
import numpy as np
import math
import time
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

def func1(x, k1, poolsize):
    return poolsize * (1 - np.exp(-k1 * x))

def func_r_b(x, kb, kc, k3, flux_in):
    a = kb + kc + k3
    b = kb * k3
    lamda1 = (-a + np.sqrt(a**2 - 4*b)) / 2
    lamda2 = (-a - np.sqrt(a**2 - 4*b)) / 2
    return flux_in * (k3 + kc) / (kb * k3) + flux_in / (lamda1 - lamda2) * (
        (lamda1 + k3 + kc) / lamda1 * np.exp(lamda1 * x) -
        (lamda2 + k3 + kc) / lamda2 * np.exp(lamda2 * x)
    )

def func_r_c(x, kb, kc, k3, flux_in):
    a = kb + kc + k3
    b = kb * k3
    lamda1 = (-a + np.sqrt(a**2 - 4*b)) / 2
    lamda2 = (-a - np.sqrt(a**2 - 4*b)) / 2
    return flux_in / k3 + flux_in * kb / (lamda1 - lamda2) * (
        1 / lamda1 * np.exp(lamda1 * x) - 1 / lamda2 * np.exp(lamda2 * x)
    )

# 设置随机种子
rng = np.random.Generator(np.random.PCG64(int(time.time())))

# 固定参数
k3 = 1
flux_in = 100
RSD = 0.1
delta_error = 0.3
Standard = 100

timepoint = np.array([0.0167, 0.05, 0.0833, 0.1667, 0.5, 1, 3, 6])
timepoint_extended = np.tile(timepoint, 9)

kb_values = np.logspace(-1, 2, 50)
kc_values = np.logspace(-1, 2, 50)

kb_grid, kc_grid = np.meshgrid(kb_values, kc_values)
n_kb = len(kb_values)
n_kc = len(kc_values)

f_single_b_grid = np.zeros((n_kc, n_kb))
f_group_grid = np.zeros((n_kc, n_kb))

total_points = n_kb * n_kc
processed = 0
start_time = time.time()

for i, kb in enumerate(kb_values):
    for j, kc in enumerate(kc_values):
        processed += 1
        if processed % 100 == 0:
            print(f"Processed {processed}/{total_points}, time: {time.time()-start_time:.1f}s")

        try:
            yb_list = []
            for _ in range(9):
                yb = np.array([
                    func_r_b(t, kb, kc, k3, flux_in) +
                    rng.normal(0, RSD * func_r_b(t, kb, kc, k3, flux_in))
                    for t in timepoint
                ])
                yb_list.append(yb)
            yb_combined = np.concatenate(yb_list)
        except Exception as e:
            print(f"Error generating B data (kb={kb:.3f}, kc={kc:.3f}): {e}")
            continue

        try:
            popt_b, _ = curve_fit(func1, timepoint_extended, yb_combined, p0=[0.5, 100], maxfev=5000)
            kb_app = popt_b[0]
            poolsize_b_app = popt_b[1]
        except:
            kb_app = 1
            poolsize_b_app = 100
        f_single_b = kb_app * poolsize_b_app

        try:
            yc_list = []
            for _ in range(9):
                yc = np.array([
                    func_r_c(t, kb, kc, k3, flux_in) +
                    rng.normal(0, RSD * func_r_c(t, kb, kc, k3, flux_in))
                    for t in timepoint
                ])
                yc_list.append(yc)
            yc_combined = np.concatenate(yc_list)
        except Exception as e:
            print(f"Error generating C data (kb={kb:.3f}, kc={kc:.3f}): {e}")
            continue

        try:
            popt_c, _ = curve_fit(func1, timepoint_extended, yc_combined, p0=[0.5, 100], maxfev=5000)
            kc_app = popt_c[0]
            poolsize_c_app = popt_c[1]
        except:
            kc_app = 1
            poolsize_c_app = 100
        f_single_c = kc_app * poolsize_c_app

        ybc_list = []
        for rep in range(9):
            yb = yb_list[rep]
            yc = yc_list[rep]
            ybc = yb + yc
            ybc_list.append(ybc)
        ybc_combined = np.concatenate(ybc_list)

        try:
            popt_bc, _ = curve_fit(func1, timepoint_extended, ybc_combined, p0=[0.5, 200], maxfev=5000)
            kbc_app = popt_bc[0]
            poolsize_bc_app = popt_bc[1]
            f_group = kbc_app * poolsize_bc_app
        except:
            f_group = 100

        f_single_b_grid[j, i] = f_single_b
        f_group_grid[j, i] = f_group

print(f"All points processed, total time: {time.time() - start_time:.1f}s")

# 准确度计算
def accuracy(value, Standard=100):
    return 100 - 100*abs(value - Standard)/Standard

acc_single_b = np.vectorize(accuracy)(f_single_b_grid)
acc_group = np.vectorize(accuracy)(f_group_grid)

acc_single_b = np.clip(acc_single_b, 0, 100)
acc_group = np.clip(acc_group, 0, 100)

# 自定义 colormap
def create_custom_cmap(base_cmap='viridis', gray_low=100*(1-delta_error), gray_high=100):
    base = plt.cm.get_cmap(base_cmap)
    colors = []
    for i in range(256):
        val = i / 255 * 100
        if gray_low <= val <= gray_high:
            colors.append((0.5, 0.5, 0.5, 1.0))
        else:
            colors.append(base(i / 255))
    return plt.matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', colors)

cmap_custom = create_custom_cmap('viridis', 100*(1-delta_error), 100)

# ==========绘图==========
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 40
fig = plt.figure(figsize=(18, 7))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.08])

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
cax = fig.add_subplot(gs[2])   # colorbar 的轴

pc1 = ax1.pcolormesh(np.log10(kb_grid), np.log10(kc_grid), acc_single_b,
                     cmap=cmap_custom, shading='auto', vmin=0, vmax=100)


pc2 = ax2.pcolormesh(np.log10(kb_grid), np.log10(kc_grid), acc_group,
                     cmap=cmap_custom, shading='auto', vmin=0, vmax=100)


cbar = fig.colorbar(pc1, cax=cax)
cbar.ax.tick_params(labelsize=36)
cbar.set_label('Accuracy (%)', fontsize=40)
plt.tight_layout()
plt.show()

# 使用 tight_layout 并增加边距和子图间距
plt.tight_layout()

plt.show()
