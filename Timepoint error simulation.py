import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd 
import time
from scipy.optimize import curve_fit
from matplotlib.font_manager import FontProperties

# 定义函数
def func(x, k1):
    return (1 - np.exp(-k1 * x))

def func1(x, k1, poolsize):
    return poolsize * (1 - np.exp(-k1 * x))

# 线性函数定义
def linear_func(x, k):
    return k * x

# 新的拟合函数：y = 100*ln(1+E*e^(kx)) / (kx)
def custom_fit_func(x, E, k):
    """
    自定义拟合函数：y = ln(1 + E * exp(k*x)) / (k*x)
    注意：当x接近0时，使用极限值避免除以0
    """
    # 当x非常接近0时，使用极限值
    if np.isscalar(x):
        if abs(x) < 1e-10:
            return 100*np.log(1 + E) / (k * 1e-10) if E > 0 else 0
        else:
            return 100*np.log(1 + E * np.exp(k * x)) / (k * x)
    else:
        # 处理数组
        result = np.zeros_like(x, dtype=float)
        for i, xi in enumerate(x):
            if abs(xi) < 1e-10:
                result[i] = 100*np.log(1 + E) / (k * 1e-10) if E > 0 else 0
            else:
                result[i] = 100*np.log(1 + E * np.exp(k * xi)) / (k * xi)
        return result

# 设置随机种子
rng = np.random.Generator(np.random.PCG64(int(time.time())))

# 参数设置
RSD = 0.05  # 相对标准误
poolsize = 100
threshold = 0.5  # 阈值

# 时间点、kapp值设置
time_OnePoint = np.append(np.arange(0.01, 2, 0.01), np.arange(2, 12, 0.02))
kapp_point = np.append(np.arange(0.01, 2, 0.02), np.arange(2, 10, 0.05))

# 存储所有结果的数据结构
all_results = []  # 每个元素是一个字典，包含一个kapp的所有时间点的拟合结果

# 定义最小二乘法计算斜率的函数
def least_squares_slope(x_points, y_points):
    """
    使用最小二乘法计算斜率 (y = kx)
    假设模型为过原点的直线
    """
    # 计算斜率 k = Σ(x*y) / Σ(x²)
    numerator = sum(x * y for x, y in zip(x_points, y_points))
    denominator = sum(x * x for x in x_points)
    
    if denominator == 0:
        return 0, 0
    
    k = numerator / denominator
    
    # 计算残差和标准差
    residuals = [y - k * x for x, y in zip(x_points, y_points)]
    n = len(residuals)
    
    if n <= 1:
        return k, 0
    
    # 计算残差平方和
    ss_res = sum(r * r for r in residuals)
    
    # 计算斜率的标准误差
    # 公式: SE(k) = sqrt(SS_res / ((n-1) * Σ(x²)))
    se_k = np.sqrt(ss_res / ((n - 1) * denominator))
    
    return k, se_k

# 存储最佳时间点
best_linear_times = []  # 线性拟合的最佳时间点
best_exp_times = []     # 指数拟合的最佳时间点
best_linear_kapps = []  # 对应的kapp值
best_exp_kapps = []     # 对应的kapp值
best_linear_diffs = []  # 最小差距
best_exp_diffs = []     # 最小差距

# 对每个kapp值进行循环
for kapp in kapp_point:
    print(f"\n处理 kapp = {kapp:.2f}")
    
    # 生成多组数据
    y_lists = []
    for _ in range(9):  # 9组数据
        y_lists.append([poolsize*(func(i, kapp) + rng.normal(loc=0.0, scale=RSD)) for i in time_OnePoint])
    
    # 存储这个kapp下每个时间点的结果
    kapp_results = {
        'kapp': kapp,
        'time_points': time_OnePoint.copy(),
        'linear_flux': [],
        'linear_flux_err': [],
        'linear_diff_percent': [],
        'exp_flux': [],
        'exp_flux_err': [],
        'exp_diff_percent': []
    }
    
    # 真实flux值
    true_flux = poolsize * kapp
    
    # 对每个时间点进行循环
    for i in range(len(time_OnePoint)):
        t = time_OnePoint[i]
        
        # 获取当前时间点的数据点
        x_points = [t] * 9
        y_points = [y_list[i] for y_list in y_lists]
        
        # 方法1：线性拟合（最小二乘法直接计算斜率）
        try:
            k_linear, se_k_linear = least_squares_slope(x_points, y_points)
            flux_linear = k_linear
            flux_linear_err = se_k_linear
        except:
            flux_linear = 0
            flux_linear_err = 0
        
        # 方法2：改进的func1拟合方法
        # 使用加权最小二乘法，考虑变换后的方差
        try:
            # 计算变换后的值：z = -ln(1-y/poolsize)
            z_points = []
            weights = []  # 权重，基于变换后的方差
            
            for y in y_points:
                # 避免对数中的负数或零
                ratio = y / poolsize
                
                # 确保ratio在合理的范围内
                if ratio >= 0.999:
                    ratio = 0.999
                elif ratio <= 0.001:
                    ratio = 0.001
                
                # 计算变换后的值
                z = -np.log(1 - ratio)
                z_points.append(z)
                
                # 计算变换后的近似方差
                # 使用误差传播公式：Var(z) ≈ (1/(poolsize*(1-ratio))²) * Var(y)
                # 假设Var(y) = (RSD * poolsize)²，但这里我们使用个体y值的方差估计
                # 简化：权重与 (1-ratio)² 成反比
                weight = (1 - ratio)**2
                weights.append(weight)
            
            # 加权最小二乘法计算斜率
            # 公式：k = Σ(w_i * x_i * z_i) / Σ(w_i * x_i²)
            w_sum_xy = sum(weights[j] * x_points[j] * z_points[j] for j in range(len(x_points)))
            w_sum_xx = sum(weights[j] * x_points[j] * x_points[j] for j in range(len(x_points)))
            
            if w_sum_xx == 0:
                k_exp = kapp
                se_k_exp = 0
            else:
                k_exp = w_sum_xy / w_sum_xx
                
                # 计算残差
                residuals = [z_points[j] - k_exp * x_points[j] for j in range(len(x_points))]
                
                # 计算加权残差平方和
                w_rss = sum(weights[j] * residuals[j]**2 for j in range(len(residuals)))
                
                # 计算标准误差
                n = len(x_points)
                if n > 1:
                    se_k_exp = np.sqrt(w_rss / ((n - 1) * w_sum_xx))
                else:
                    se_k_exp = 0
            
            flux_exp = poolsize * k_exp
            flux_exp_err = poolsize * se_k_exp
            
        except Exception as e:
            #print(f"指数拟合出错: {e}")
            flux_exp = 0
            flux_exp_err = 0
            k_exp = 0

        # 计算百分比差距
        linear_diff_percent = abs((flux_linear - true_flux) / true_flux) * 100
        exp_diff_percent = abs((k_exp - kapp) / kapp) * 100
        
        # 存储结果
        kapp_results['linear_flux'].append(flux_linear)
        kapp_results['linear_flux_err'].append(flux_linear_err)
        kapp_results['linear_diff_percent'].append(linear_diff_percent)
        
        kapp_results['exp_flux'].append(flux_exp)
        kapp_results['exp_flux_err'].append(flux_exp_err)
        kapp_results['exp_diff_percent'].append(exp_diff_percent)
    
    # 转换为numpy数组
    for key in ['linear_flux', 'linear_flux_err', 'linear_diff_percent', 
                'exp_flux', 'exp_flux_err', 'exp_diff_percent']:
        kapp_results[key] = np.array(kapp_results[key])
    
    # 存储这个kapp的结果
    all_results.append(kapp_results)
    
    # 找到线性拟合的最小差距时间点
    min_linear_idx = np.argmin(kapp_results['linear_diff_percent'])
    min_linear_diff = kapp_results['linear_diff_percent'][min_linear_idx]
    min_linear_time = kapp_results['time_points'][min_linear_idx]
    
    # 找到指数拟合的最小差距时间点
    min_exp_idx = np.argmin(kapp_results['exp_diff_percent'])
    min_exp_diff = kapp_results['exp_diff_percent'][min_exp_idx]
    min_exp_time = kapp_results['time_points'][min_exp_idx]
    
    # 只记录差距小于阈值的最佳时间点
    if min_linear_diff < threshold:
        best_linear_times.append(min_linear_time)
        best_linear_kapps.append(kapp)
        best_linear_diffs.append(min_linear_diff)
        print(f"  线性拟合最佳时间点: {min_linear_time:.3f}, 最小差距: {min_linear_diff:.2f}%")
    
    if min_exp_diff < threshold:
        best_exp_times.append(min_exp_time)
        best_exp_kapps.append(kapp)
        best_exp_diffs.append(min_exp_diff)
        print(f"  指数拟合最佳时间点: {min_exp_time:.3f}, 最小差距: {min_exp_diff:.2f}%")

# 第一个图：特定kapp下拟合差距随时间的变化
kapps_to_plot = [0.1, 0.2, 0.5, 1, 5]
fig1, axes1 = plt.subplots(len(kapps_to_plot), 1, figsize=(12, 2*len(kapps_to_plot)))

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

if len(kapps_to_plot) == 1:
    axes1 = [axes1]

# 存储拟合参数
fit_params_linear = {}  # 线性拟合差距的拟合参数
fit_params_exp = {}     # 指数拟合差距的拟合参数

for idx, target_kapp in enumerate(kapps_to_plot):
    # 找到对应的kapp结果
    kapp_result = None
    for result in all_results:
        if abs(result['kapp'] - target_kapp) < 0.1:  # 允许小的浮点误差
            kapp_result = result
            break
    
    if kapp_result is None:
        print(f"未找到kapp={target_kapp}的结果")
        continue
    
    ax = axes1[idx]
    
    # 获取数据
    x_data = kapp_result['time_points']
    y_linear_data = kapp_result['linear_diff_percent']
    y_exp_data = kapp_result['exp_diff_percent']
    
    # 绘制散点图
    ax.scatter(x_data, y_linear_data, c='blue', s=6, label='linear fit', alpha=0.2)
    ax.scatter(x_data, y_exp_data, c='red', s=6, label='exponential fit', alpha=0.2)
    
    # 对线性拟合差距数据进行拟合
    try:
        # 初始参数猜测
        p0_linear = [0.5, 10]
        
        # 进行曲线拟合
        # 注意：我们需要过滤掉无效数据（NaN或无限大）
        valid_indices_linear = np.isfinite(y_linear_data) & (y_linear_data >= 0)
        if np.sum(valid_indices_linear) > 10:  # 至少需要10个有效点
            x_fit_linear = x_data[valid_indices_linear]
            y_fit_linear = y_linear_data[valid_indices_linear]
            
            # 设置参数范围：E>0, k可以为负数
            bounds_linear = ([0.00001, 0], [0.99999, 10000])
            
            popt_linear, pcov_linear = curve_fit(custom_fit_func, x_fit_linear, y_fit_linear, 
                                                 p0=p0_linear, bounds=bounds_linear, maxfev=5000)
            
            E_linear, k_linear = popt_linear
            
            # 生成拟合曲线
            x_fit_smooth = np.linspace(x_data.min(), x_data.max(), 200)
            y_fit_smooth = custom_fit_func(x_fit_smooth, E_linear, k_linear)
            
            # 绘制拟合曲线
            ax.plot(x_fit_smooth, y_fit_smooth, 'b--', linewidth=2)
            
            # 存储拟合参数
            fit_params_linear[target_kapp] = (E_linear, k_linear)
            
            # 计算R²
            y_pred = custom_fit_func(x_fit_linear, E_linear, k_linear)
            ss_res = np.sum((y_fit_linear - y_pred) ** 2)
            ss_tot = np.sum((y_fit_linear - np.mean(y_fit_linear)) ** 2)
            r2_linear = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            print(f"kapp={target_kapp:.1f} 线性拟合差距拟合: E={E_linear:.4f}, k={k_linear:.4f}, R²={r2_linear:.4f}")
        else:
            print(f"kapp={target_kapp:.1f} 线性拟合差距有效数据不足")
    except Exception as e:
        print(f"kapp={target_kapp:.1f} 线性拟合差距拟合失败: {e}")
    
    # 对指数拟合差距数据进行拟合
    try:
        # 初始参数猜测
        p0_exp = [0.5, 10]
        
        # 进行曲线拟合
        valid_indices_exp = np.isfinite(y_exp_data) & (y_exp_data >= 0)
        if np.sum(valid_indices_exp) > 10:  # 至少需要10个有效点
            x_fit_exp = x_data[valid_indices_exp]
            y_fit_exp = y_exp_data[valid_indices_exp]
            
            # 设置参数范围：E>0, k可以为负数
            bounds_exp = ([0.00001, 0], [0.99999, 10000])
            
            popt_exp, pcov_exp = curve_fit(custom_fit_func, x_fit_exp, y_fit_exp, 
                                           p0=p0_exp, bounds=bounds_exp, maxfev=5000)
            
            E_exp, k_exp = popt_exp
            
            # 生成拟合曲线
            x_fit_smooth = np.linspace(x_data.min(), x_data.max(), 200)
            y_fit_smooth = custom_fit_func(x_fit_smooth, E_exp, k_exp)
            
            # 绘制拟合曲线
            ax.plot(x_fit_smooth, y_fit_smooth, 'r--', linewidth=2)
            
            # 存储拟合参数
            fit_params_exp[target_kapp] = (E_exp, k_exp)
            
            # 计算R²
            y_pred = custom_fit_func(x_fit_exp, E_exp, k_exp)
            ss_res = np.sum((y_fit_exp - y_pred) ** 2)
            ss_tot = np.sum((y_fit_exp - np.mean(y_fit_exp)) ** 2)
            r2_exp = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            print(f"kapp={target_kapp:.1f} 指数拟合差距拟合: E={E_exp:.4f}, k={k_exp:.4f}, R²={r2_exp:.4f}")
        else:
            print(f"kapp={target_kapp:.1f} 指数拟合差距有效数据不足")
    except Exception as e:
        print(f"kapp={target_kapp:.1f} 指数拟合差距拟合失败: {e}")
    
    ax.set_xlabel('timepoint(h)')
    ax.set_ylabel('Fitting deviation (%)')
    ax.set_title(f'kapp = {target_kapp:.1f}, Fitting deviation vs timepoint(h)', pad=20)
    
    # 将图例放在图形外部，避免遮挡数据
    ax.legend(loc='upper right', fontsize='small', bbox_to_anchor=(1, 1))
    
    ax.grid(True, alpha=0.3)
    
    # 设置y轴范围
    y_max = max(np.nanmax(kapp_result['linear_diff_percent']), 
                np.nanmax(kapp_result['exp_diff_percent'])) * 1.1
    ax.set_ylim([0, min(y_max, 200)])

# 调整子图之间的间距，避免标题重叠
plt.subplots_adjust(hspace=0.4)
plt.tight_layout(rect=[0, 0, 0.9, 1])  # 调整布局，为图例留出空间
plt.show()

# 第二个图：最佳时间点与kapp的关系
fig2, ax2 = plt.subplots(figsize=(12, 8))
font = FontProperties()
font.set_family('sans-serif') 
font.set_name('Arial')
font.set_size(40)

# 绘制散点图
if len(best_linear_kapps) > 0:
    ax2.scatter(best_linear_kapps, best_linear_times, alpha=0.6, s=30, 
               color='blue')
    
if len(best_exp_kapps) > 0:
    ax2.scatter(best_exp_kapps, best_exp_times, alpha=0.6, s=30, 
               color='red', marker='s')

# 进行反比例函数拟合（t = a/k，无截距）
# 线性拟合的最佳时间点拟合
if len(best_linear_kapps) > 2:
    try:
        # 直接计算 a，使用列表推导式而不是生成器
        a_linear = np.mean([best_linear_times[i] * best_linear_kapps[i] for i in range(len(best_linear_times))])
        
        # 生成拟合曲线
        k_fit_linear = np.linspace(0.07, max(best_linear_kapps), 200)
        t_fit_linear = a_linear / k_fit_linear
        
        ax2.plot(k_fit_linear, t_fit_linear, 'b-', alpha=0.7, linewidth=2,
                label=f'linear fit: t = {a_linear:.2f}/k')
        
        # 计算R²
        t_pred_linear = a_linear / np.array(best_linear_kapps)
        ss_res = sum([(best_linear_times[i] - t_pred_linear[i]) ** 2 for i in range(len(best_linear_times))])
        ss_tot = sum([(best_linear_times[i] - np.mean(best_linear_times)) ** 2 for i in range(len(best_linear_times))])
        r2_linear = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        print(f"\n线性拟合反比例函数拟合结果:")
        print(f"  拟合参数: a = {a_linear:.4f}")
        print(f"  R² = {r2_linear:.4f}")
    except Exception as e:
        print(f"线性拟合反比例函数拟合失败: {e}")

# 指数拟合的最佳时间点拟合
if len(best_exp_kapps) > 2:
    try:
        # 直接计算 a，使用列表推导式而不是生成器
        a_exp = np.mean([best_exp_times[i] * best_exp_kapps[i] for i in range(len(best_exp_times))])
        
        # 生成拟合曲线
        k_fit_exp = np.linspace(0.07, max(best_exp_kapps), 200)
        t_fit_exp = a_exp / k_fit_exp
        
        ax2.plot(k_fit_exp, t_fit_exp, 'r-', alpha=0.7, linewidth=2,
                label=f'exponential fit: t = {a_exp:.2f}/k')
        
        # 计算R²
        t_pred_exp = a_exp / np.array(best_exp_kapps)
        ss_res = sum([(best_exp_times[i] - t_pred_exp[i]) ** 2 for i in range(len(best_exp_times))])
        ss_tot = sum([(best_exp_times[i] - np.mean(best_exp_times)) ** 2 for i in range(len(best_exp_times))])
        r2_exp = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        print(f"\n指数拟合反比例函数拟合结果:")
        print(f"  拟合参数: a = {a_exp:.4f}")
        print(f"  R² = {r2_exp:.4f}")
    except Exception as e:
        print(f"指数拟合反比例函数拟合失败: {e}")

ax2.xaxis.set_tick_params(labelsize=40)
ax2.yaxis.set_tick_params(labelsize=40)
ax2.set_xlabel('kapp value', fontproperties=font)
ax2.set_ylabel('Best timepoint(h)', fontproperties=font)
ax2.set_title('Best timepoint vs kapp value', pad=15, fontproperties=font)
ax2.grid(True, alpha=0.3)

# 将图例放在图形外部，避免遮挡数据
ax2.legend(loc='upper right', bbox_to_anchor=(1, 1), prop=font)

plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.show()

# 打印拟合参数汇总
print("\n=== 第一个图的拟合参数汇总 ===")
print("线性拟合差距的拟合参数:")
for kapp, params in fit_params_linear.items():
    E, k = params
    print(f"  kapp={kapp}: E={E:.4f}, k={k:.4f}")

print("\n指数拟合差距的拟合参数:")
for kapp, params in fit_params_exp.items():
    E, k = params
    print(f"  kapp={kapp}: E={E:.4f}, k={k:.4f}")

# 计算并显示整体统计信息
print("\n=== 整体统计信息 ===")
print(f"kapp范围: {kapp_point.min():.2f} - {kapp_point.max():.2f}")
print(f"时间点范围: {time_OnePoint.min():.2f} - {time_OnePoint.max():.2f}")
print(f"阈值: {threshold}%")

print("\n模拟完成!")