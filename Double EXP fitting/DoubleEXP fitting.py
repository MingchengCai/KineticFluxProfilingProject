import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd 
from scipy.optimize import curve_fit

# 定义函数（使用 np.exp 支持数组运算）
def func1(x, k1):
    return 1 - np.exp(-k1 * x)

def func2(x, k1, k2):
    # 双指数模型：1 + k2/(k1-k2)*exp(-k1*x) - k1/(k1-k2)*exp(-k2*x)
    # 注意处理 k1 == k2 的情况（此处假设不等）
    return 1 + k2/(k1-k2) * np.exp(-k1*x) - k1/(k1-k2) * np.exp(-k2*x)

# 读取数据
data = pd.read_csv(r'k_value_fitting.csv')

# 读取表头
head = data.columns.values.tolist()

# 忽略除以0的警告
np.seterr(divide='ignore', invalid='ignore')

n = 6           # 每组平行的数量
m = 1
col = 0         # 代谢物计数
len_col = len(data.iloc[0])         # 总列数
ms_mtx = pd.DataFrame()
column = []

# 计算每 n 列的均值和标准差
while m < len_col:
    if m % n == 1:
        mean = pd.DataFrame(np.mean(data.iloc[:, m:m+n], axis=1))
        std = pd.DataFrame(np.std(data.iloc[:, m:m+n], axis=1))
        column.append((head[m], "mean"))
        column.append((head[m], "std"))
        mtx = pd.concat([mean, std], axis=1)
        ms_mtx = pd.concat([ms_mtx, mtx], axis=1)
        col += 1
    m += 1

ms_mtx.columns = pd.MultiIndex.from_tuples(column)
ms_mtx.insert(loc=0, column='time(h)', value=data.iloc[:, 0])   # 插入时间列
outputpath = 'result_mean_std.csv'
ms_mtx.to_csv(outputpath, sep=',', index=True, header=True)


time = ms_mtx['time(h)'].values   # 时间点数组
num_metabolites = col             # 代谢物数量

# 用于存储拟合结果的列表
results = []  # 每个元素为 (代谢物名, k1, k2)

# 遍历每个代谢物（ms_mtx 中每两列：mean 和 std，我们只使用 mean 列）
for i in range(num_metabolites):
    # 均值列位于第 2*i+1 列（因第0列为 time）
    mean_col = ms_mtx.iloc[:, 2*i+1]
    met_name = ms_mtx.columns[2*i+1][0]   # 代谢物原始名称
    
    y_data = mean_col.values
    
    # 拟合前检查数据：若全为0或全为1等极端情况，跳过拟合
    if np.all(y_data == 0) or np.all(y_data == 1):
        print(f"警告：代谢物 {met_name} 数据无变化，跳过拟合")
        results.append([met_name, np.nan, np.nan])
        continue
    
    # 初始猜测：k1=1, k2=10（可根据实际情况调整）
    p0 = [1.0, 0.5]
    # 边界：k1 > 0, k2 > 0
    bounds = (0, np.inf)
    
    try:
        popt, pcov = curve_fit(func2, time, y_data, p0=p0, bounds=bounds, maxfev=5000)
        k1, k2 = popt
        results.append([met_name, k1, k2])
    except Exception as e:
        print(f"拟合失败：代谢物 {met_name}，错误：{e}")
        results.append([met_name, np.nan, np.nan])

# 将结果保存为 DataFrame 并输出 CSV
fit_df = pd.DataFrame(results, columns=['Metabolite', 'k1', 'k2'])
fit_df.to_csv('fit_k_values.csv', sep=',', index=False)

print("拟合完成，结果已保存至 fit_k_values.csv")
