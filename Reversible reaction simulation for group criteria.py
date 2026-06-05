import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd 
from scipy.optimize import curve_fit
import time
from scipy import stats
from matplotlib.font_manager import FontProperties

# Define functions
def func1(x, k1, poolsize):
    return poolsize * (1 - np.exp(-k1 * x))

def func2(x, k2, k3, poolsize):
    if abs(k2 - k3) < 1e-6:
        return poolsize * (1 - (1 + k2 * x) * np.exp(-k2 * x))
    else:
        return poolsize * (1 - k3/(k3-k2)*np.exp(-k2*x) - k2/(k2-k3)*np.exp(-k3*x))

def func3(x,k2,k3,k4,poolsize):
    return poolsize*(1-k3*k4/((k3-k2)*(k4-k2))*(math.e)**(-k2*x)-k2*k4/((k2-k3)*(k4-k3))*(math.e)**(-k3*x)-k2*k3/((k2-k4)*(k3-k4))*(math.e)**(-k4*x))

def func_r_b(x,kb,kc,k3,flux_in):
    a = kb+kc+k3; b = kb*k3
    lamda1 = (-a+np.sqrt(a**2-4*b))/2; lamda2= (-a-np.sqrt(a**2-4*b))/2
    return flux_in*(k3+kc)/(kb*k3)+flux_in/(lamda1-lamda2)*((lamda1+k3+kc)/lamda1*math.e**(lamda1*x)-(lamda2+k3+kc)/lamda2*math.e**(lamda2*x))

def func_r_c(x,kb,kc,k3,flux_in):
    a = kb+kc+k3; b = kb*k3
    lamda1 = (-a+(a**2-4*b)**0.5)/2; lamda2= (-a-(a**2-4*b)**0.5)/2
    return flux_in/k3+flux_in*kb/(lamda1-lamda2)*(1/lamda1*math.e**(lamda1*x)-1/lamda2*math.e**(lamda2*x))

# Set random seed
rng = np.random.Generator(np.random.PCG64(int(time.time())))

# Parameter settings
k3 = 1
flux_in = 100
RSD = 0.1  # Relative standard error
delta_error = 0.2 # Acceptable error range
delta_sg = 0.8  # Ratio threshold: f_single_B vs f_group
kc_threshold_log = -1  # kc threshold (log10 value)

# Time points setting
timepoint = np.array([0.0167, 0.05, 0.0833, 0.1667, 0.5, 1, 3, 6])
timepoint_extended = np.tile(timepoint, 9)  # Extend time point array to match data point count

# Set kb and kc value ranges
kb_values = np.logspace(-1, 2, 50)  # From 0.1 to 100, take 40 logarithmically spaced values
kc_values = np.logspace(-1, 2, 50)  # From 1 to 100, take 10 logarithmically spaced values


# Create grid
kb_grid, kc_grid = np.meshgrid(kb_values, kc_values)

# Store results in 2D arrays
f_group_grid = np.zeros_like(kb_grid)
f_single_b_grid = np.zeros_like(kb_grid)
f_single_c_grid = np.zeros_like(kb_grid)
kb_app_grid = np.zeros_like(kb_grid)
kc_app_grid = np.zeros_like(kc_grid)

# Store boundary combinations
boundary_combinations = []  # Store (kb, kc) boundary combinations
boundary_app_combinations = []  # Store (kb_app, kc_app) boundary combinations

# For storing all points for new scatter plots
all_points_kb_kc = []  # Store all (kb, kc) points
all_points_kb_app_kc_app = []  # Store all (kb_app, kc_app) points
all_points_conditions = []  # Store condition for each point: -1 for f_single_B < f_group*delta_sg, 1 for f_single_B > f_group*delta_sg

# For storing all points for the new scatter plots (无视kc_threshold_log限制)
all_points_all = []  # Store all (kb_app, kc_app, f_single_B, f_single_C, f_group) for all points

# Process each kb and kc value
total_points = kb_values.size * kc_values.size
processed_points = 0
start_time = time.time()

for i, kb in enumerate(kb_values):
    for j, kc in enumerate(kc_values):
        processed_points += 1
        if processed_points % 100 == 0:
            elapsed_time = time.time() - start_time
            print(f"Processed {processed_points}/{total_points} points, time elapsed: {elapsed_time:.1f} seconds")
        
        # Generate B data - nine replicates
        try:
            yb_1 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_2 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_3 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_4 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_5 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_6 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_7 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_8 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
            yb_9 = np.array([func_r_b(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_b(t, kb, kc, k3, flux_in)) for t in timepoint])
        except Exception as e:
            print(f"Error generating B data (kb={kb:.3f}, kc={kc:.3f}): {e}")
            continue
        
        # Fit B data - combine nine sets of data for fitting
        yb_combined = np.concatenate([yb_1, yb_2, yb_3, yb_4, yb_5, yb_6, yb_7, yb_8, yb_9])
        
        try:
            popt_b, pcov_b = curve_fit(func1, timepoint_extended, yb_combined, p0=[0.5, 100], maxfev=5000)
            kb_app = popt_b[0]
            poolsize_b_app = popt_b[1]
        except Exception as e:
            # If fitting fails, use default values
            kb_app = 1
            poolsize_b_app = 100
        
        # Calculate f_single for B
        f_single_b = kb_app * poolsize_b_app
        
        # Generate C data - nine replicates
        try:
            yc_1 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_2 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_3 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_4 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_5 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_6 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_7 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_8 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
            yc_9 = np.array([func_r_c(t, kb, kc, k3, flux_in) + rng.normal(loc=0.0, scale=RSD*func_r_c(t, kb, kc, k3, flux_in)) for t in timepoint])
        except Exception as e:
            print(f"Error generating C data (kb={kb:.3f}, kc={kc:.3f}): {e}")
            continue

        # Fit C data - combine nine sets of data for fitting
        yc_combined = np.concatenate([yc_1, yc_2, yc_3, yc_4, yc_5, yc_6, yc_7, yc_8, yc_9])
        
        try:
            popt_c, pcov_c = curve_fit(func1, timepoint_extended, yc_combined, p0=[0.5, 100], maxfev=5000)
            kc_app = popt_c[0]
            poolsize_c_app = popt_c[1]
        except Exception as e:
            # If fitting fails, use default values
            kc_app = 1
            poolsize_c_app = 100
        
        # Calculate f_single for C
        f_single_c = kc_app * poolsize_c_app
        
        # Add B and C data
        ybc_1 = yb_1 + yc_1
        ybc_2 = yb_2 + yc_2
        ybc_3 = yb_3 + yc_3
        ybc_4 = yb_4 + yc_4
        ybc_5 = yb_5 + yc_5
        ybc_6 = yb_6 + yc_6
        ybc_7 = yb_7 + yc_7
        ybc_8 = yb_8 + yc_8
        ybc_9 = yb_9 + yc_9
        
        # Fit B+C data - combine nine sets of data for fitting
        ybc_combined = np.concatenate([ybc_1, ybc_2, ybc_3, ybc_4, ybc_5, ybc_6, ybc_7, ybc_8, ybc_9])
        
        try:
            popt_bc, pcov_bc = curve_fit(func1, timepoint_extended, ybc_combined, p0=[0.5, 200], maxfev=5000)
            kbc_app = popt_bc[0]
            poolsize_bc_app = popt_bc[1]
            f_group = kbc_app * poolsize_bc_app
            
        except Exception as e:
            # If fitting fails, use default value
            f_group = 100
        
        # Store to grid
        f_group_grid[j, i] = f_group
        f_single_b_grid[j, i] = f_single_b
        f_single_c_grid[j, i] = f_single_c
        kb_app_grid[j, i] = kb_app
        kc_app_grid[j, i] = kc_app
        
        # Store all points for new scatter plots (only if log10(kc) > kc_threshold_log)
        if np.log10(kc) > kc_threshold_log:
            all_points_kb_kc.append((kb, kc))
            all_points_kb_app_kc_app.append((kb_app, kc_app))
            
            # Determine condition: -1 for f_single_B < f_group*delta_sg, 1 for f_single_B > f_group*delta_sg
            if f_single_b < f_group * delta_sg:
                all_points_conditions.append(-1)  # Gray
            else:
                all_points_conditions.append(1)   # Light red
        
        # Store all points for the new scatter plots (无视kc_threshold_log限制)
        all_points_all.append((kb_app, kc_app, f_single_b, f_single_c, f_group))

print(f"All points processed, total time: {time.time() - start_time:.1f} seconds")

# Calculate boundary value combinations
print(f"\nCalculating boundary value combinations (only for log10(kc) > {kc_threshold_log})...")
boundary_start_time = time.time()

# Calculate indices of kc that meet the condition
valid_kc_indices = [j for j, kc in enumerate(kc_values) if np.log10(kc) > kc_threshold_log]
print(f"Number of kc values meeting condition: {len(valid_kc_indices)}/{len(kc_values)}")

for j in valid_kc_indices:
    kc = kc_values[j]
    # For each kc, find kb that minimizes |f_single_B - f_group * delta_sg|
    min_diff = float('inf')
    best_kb_idx = -1
    
    for i, kb in enumerate(kb_values):
        f_single_b = f_single_b_grid[j, i]
        f_group_val = f_group_grid[j, i]
        diff = abs(f_single_b - f_group_val * delta_sg)
        
        if diff < min_diff:
            min_diff = diff
            best_kb_idx = i
    
    if best_kb_idx >= 0:
        # Record boundary combination
        kb_boundary = kb_values[best_kb_idx]
        kc_boundary = kc
        kb_app_boundary = kb_app_grid[j, best_kb_idx]
        kc_app_boundary = kc_app_grid[j, best_kb_idx]
        
        boundary_combinations.append((kb_boundary, kc_boundary))
        boundary_app_combinations.append((kb_app_boundary, kc_app_boundary))

print(f"Boundary value calculation complete, found {len(boundary_combinations)} boundary points, time: {time.time() - boundary_start_time:.1f} seconds")

# Set font
plt.rcParams['font.family'] = 'sans-serif'           
plt.rcParams['font.sans-serif'] = ['Arial']    
plt.rcParams['font.size'] = 40

# Create custom colormap function
def create_custom_cmap(base_cmap, vmin, vmax, lower_thresh, upper_thresh):
    """Create custom colormap, set areas within threshold range to gray"""
    # Get base colormap
    base = plt.cm.get_cmap(base_cmap)
    
    # Create new colormap
    colors = []
    # Use linear space to define colormap
    for i in range(256):
        # Map 0-255 to data range
        value = vmin + (vmax - vmin) * i / 255
        
        # If value is within threshold range, use gray
        if lower_thresh <= value <= upper_thresh:
            # Gray, RGBA format
            colors.append((0.5, 0.5, 0.5, 1.0))
        else:
            # Otherwise use base colormap
            colors.append(base(i))
    
    return plt.matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', colors)

# Define thresholds
lower_thresh = (1 - delta_error) * flux_in
upper_thresh = (1 + delta_error) *flux_in

# Calculate global min and max of all data
f_single_bc_grid = f_single_b_grid + f_single_c_grid
all_data = [f_single_b_grid, f_single_c_grid, f_single_bc_grid, f_group_grid]
global_min = min(np.nanmin(data) for data in all_data)
global_max = max(np.nanmax(data) for data in all_data)

# Create unified colormap
base_cmap = 'viridis'
unified_cmap = create_custom_cmap(base_cmap, global_min, global_max, lower_thresh, upper_thresh)

# 创建一个更大的画布
fig = plt.figure(figsize=(14, 7))

# 第一个子图：f_single_B 热图
ax1 = fig.add_subplot(121)
pc_b = ax1.pcolormesh(np.log10(kb_grid), np.log10(kc_grid), f_single_b_grid,
                       cmap=unified_cmap, shading='auto', vmin=global_min, vmax=global_max)
ax1.set_xlabel('log10(kb)', fontsize=40)
ax1.set_ylabel('log10(kc)', fontsize=40)
ax1.set_title(f'f_single_B', fontsize=40)
ax1.tick_params(axis='x', labelsize=36)
ax1.tick_params(axis='y', labelsize=36)


# 第二个子图：f_group 热图
ax4 = fig.add_subplot(122)
pc_group = ax4.pcolormesh(np.log10(kb_grid), np.log10(kc_grid), f_group_grid,
                           cmap=unified_cmap, shading='auto', vmin=global_min, vmax=global_max)
ax4.set_xlabel('log10(kb)', fontsize=40)
ax4.set_ylabel('log10(kc)', fontsize=40)
ax4.set_title(f'f_group', fontsize=40)
ax4.tick_params(axis='x', labelsize=40)
ax4.tick_params(axis='y', labelsize=40)

cbar = fig.colorbar(pc_b, ax=[ax1, ax4], shrink=0.6)
cbar.ax.tick_params(labelsize=36)
cbar.set_label('Flux (e.g., f value)', fontsize=40)

# 使用 tight_layout 并增加边距和子图间距
plt.tight_layout()

plt.show()