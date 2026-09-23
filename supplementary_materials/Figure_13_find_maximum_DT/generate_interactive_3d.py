#!/usr/bin/env python3
"""
完全独立的交互式 3D 图表生成脚本
使用 plotly 生成可在浏览器中旋转的 3D 网格图
"""

import re
from math import log2
from collections import defaultdict
import plotly.graph_objects as go

# 手动解析数据（避免依赖复杂的导入）
def parse_data():
    # 解析 trivial strategy
    trivial_data = []
    with open('trivial_strategy_vary_q.log', 'r') as f:
        content = f.read()
        for block in content.split('Build LWE  instance'):
            if not block.strip():
                continue
            
            n_match = re.search(r'n=(\d+)', block)
            q_match = re.search(r'q=(\d+)', block)
            g_match = re.search(r'G\s*=\s*([\d.]+)', block)
            
            if q_match and g_match:
                n = int(n_match.group(1)) if n_match else 240
                trivial_data.append({
                    'n': n,
                    'q': int(q_match.group(1)),
                    'G': float(g_match.group(1))
                })
    
    # 解析 optimized strategy
    optimized_data = []
    current = {}
    with open('optimized_strategy_vary_q.log', 'r') as f:
        for line in f:
            if line.startswith('n='):
                m = re.search(r'n=(\d+)\s+m=\d+\s+q=(\d+)', line)
                if m:
                    current = {'n': int(m.group(1)), 'q': int(m.group(2))}
            elif 'succ_prob' in line or 'success_rate' in line:
                m = re.search(r'(?:succ_prob|success_rate)\s*=\s*([\d\.eE\-\+]+)', line)
                if m and current:
                    current['sr'] = float(m.group(1))
            elif 'Min Cost =' in line:
                if current and 'sr' in current:
                    m = re.search(r'Min Cost = ([\d.]+)', line)
                    if m:
                        current['G'] = float(m.group(1))
                        optimized_data.append(current.copy())
                current = {}
    
    return trivial_data, optimized_data

print("正在加载数据...")
trivial_list, optimized_list = parse_data()

# 过滤 success_rate >= 0.99
opt_filtered = [d for d in optimized_list if d.get('sr', 0) >= 0.99]
print(f"找到 {len(opt_filtered)} 个优化策略数据点 (SR >= 0.99)")

# 合并数据并创建完整网格（缺失数据填充为0）
# 首先，找出所有唯一的 n 和 q 值
all_n = sorted(set([d['n'] for d in trivial_list]))
all_q = sorted(set([d['q'] for d in trivial_list]))

print(f"网格维度: {len(all_n)} 个 n 值 × {len(all_q)} 个 q 值")

# 创建数据字典以便快速查找
trivial_dict = {(d['n'], d['q']): d['G'] for d in trivial_list}
optimized_dict = {(d['n'], d['q']): d['G'] for d in opt_filtered}

# 构建完整网格
complete_grid = []
for n_val in all_n:
    for q_val in all_q:
        key = (n_val, q_val)
        
        # 如果两个策略都有数据，计算差值
        if key in trivial_dict and key in optimized_dict:
            diff = trivial_dict[key] - optimized_dict[key]
            complete_grid.append({
                'n': n_val,
                'q': q_val,
                'diff': diff,
                'has_data': True
            })
        # 如果数据缺失，填充为0
        else:
            complete_grid.append({
                'n': n_val,
                'q': q_val,
                'diff': 0.0,
                'has_data': False
            })

print(f"完整网格点数: {len(complete_grid)} (包括 {sum(1 for d in complete_grid if not d['has_data'])} 个填充点)")
print(f"实际数据点数: {sum(1 for d in complete_grid if d['has_data'])}")

# 提取坐标（使用完整网格）
X = [d['n'] for d in complete_grid]
Y = [log2(d['q']) for d in complete_grid]
Z = [d['diff'] for d in complete_grid]
has_data = [d['has_data'] for d in complete_grid]

# 找最大值（仅在真实数据中）
real_data_points = [(X[i], Y[i], Z[i]) for i in range(len(X)) if has_data[i]]
if real_data_points:
    max_point = max(real_data_points, key=lambda p: p[2])
    print(f"最大差值: {max_point[2]:.2f} 在 n={max_point[0]}, log2(q)={max_point[1]:.2f}")

# 创建 plotly 图表
fig = go.Figure()

# Create proper mesh grid lines
# Group data by n and q values
by_n = defaultdict(list)
by_q = defaultdict(list)

for i, (x, y, z) in enumerate(zip(X, Y, Z)):
    by_n[x].append((y, z, i))  # For fixed n, varying q
    by_q[y].append((x, z, i))  # For fixed q, varying n

# Draw lines along q direction (fixed n)
for n_val, points in by_n.items():
    if len(points) > 1:
        # Sort by q (which is y in our plot)
        points_sorted = sorted(points, key=lambda p: p[0])
        y_line = [p[0] for p in points_sorted]
        z_line = [p[1] for p in points_sorted]
        x_line = [n_val] * len(points_sorted)
        
        fig.add_trace(go.Scatter3d(
            x=x_line, y=y_line, z=z_line,
            mode='lines',
            line=dict(color='blue', width=3),
            showlegend=False,
            hoverinfo='skip'
        ))

# Draw lines along n direction (fixed q)
for q_val, points in by_q.items():
    if len(points) > 1:
        # Sort by n (which is x in our plot)
        points_sorted = sorted(points, key=lambda p: p[0])
        x_line = [p[0] for p in points_sorted]
        z_line = [p[1] for p in points_sorted]
        y_line = [q_val] * len(points_sorted)
        
        fig.add_trace(go.Scatter3d(
            x=x_line, y=y_line, z=z_line,
            mode='lines',
            line=dict(color='blue', width=3),
            showlegend=False,
            hoverinfo='skip'
        ))

# Data points (scatter only, no connecting line)
# 分别绘制真实数据点和填充点
X_real = [X[i] for i in range(len(X)) if has_data[i]]
Y_real = [Y[i] for i in range(len(Y)) if has_data[i]]
Z_real = [Z[i] for i in range(len(Z)) if has_data[i]]

X_fill = [X[i] for i in range(len(X)) if not has_data[i]]
Y_fill = [Y[i] for i in range(len(Y)) if not has_data[i]]
Z_fill = [Z[i] for i in range(len(Z)) if not has_data[i]]

# 真实数据点（彩色）
if X_real:
    fig.add_trace(go.Scatter3d(
        x=X_real, y=Y_real, z=Z_real,
        mode='markers',
        marker=dict(
            size=12,
            color=Z_real,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="成本差异"),
            line=dict(color='black', width=1)
        ),
        name='成本差异数据'
    ))

# 填充数据点（灰色，半透明）
if X_fill:
    fig.add_trace(go.Scatter3d(
        x=X_fill, y=Y_fill, z=Z_fill,
        mode='markers',
        marker=dict(
            size=8,
            color='lightgray',
            opacity=0.3,
            line=dict(color='gray', width=0.5)
        ),
        name='填充点 (缺失数据)'
    ))

# 添加垂直支撑线（连接到底面）
for xi, yi, zi in zip(X, Y, Z):
    fig.add_trace(go.Scatter3d(
        x=[xi, xi],
        y=[yi, yi],
        z=[0, zi],
        mode='lines',
        line=dict(color='lightgray', width=2, dash='dot'),
        showlegend=False,
        hoverinfo='skip'
    ))

# 标记最大值点（仅在真实数据中查找）
if Z_real:
    max_idx_real = Z_real.index(max(Z_real))
    max_n = X_real[max_idx_real]
    max_logq = Y_real[max_idx_real]
    max_z = Z_real[max_idx_real]
    
    fig.add_trace(go.Scatter3d(
        x=[max_n],
        y=[max_logq],
        z=[max_z],
        mode='markers+text',
        marker=dict(
            size=20,
            color='red',
            symbol='diamond',
            line=dict(color='darkred', width=2)
        ),
        text=[f"最大值: {max_z:.2f}"],
        textposition="top center",
        name='最大值点',
        textfont=dict(size=14, color='red', family='Arial Black')
    ))

# 布局设置
fig.update_layout(
    title=dict(
        text='<b>成本差异 3D 可视化</b><br><sub>用鼠标拖动旋转 | 滚轮缩放</sub>',
        x=0.5,
        xanchor='center',
        font=dict(size=20)
    ),
    scene=dict(
        xaxis_title='<b>n (维度)</b>',
        yaxis_title='<b>log₂(q)</b>',
        zaxis_title='<b>差值 (G_trivial - G_optimized)</b>',
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.3)
        ),
        xaxis=dict(
            backgroundcolor="rgb(245,245,245)",
            gridcolor="white",
            showbackground=True,
            gridwidth=2
        ),
        yaxis=dict(
            backgroundcolor="rgb(245,245,245)",
            gridcolor="white",
            showbackground=True,
            gridwidth=2
        ),
        zaxis=dict(
            backgroundcolor="rgb(245,245,245)",
            gridcolor="white",
            showbackground=True,
            gridwidth=2
        )
    ),
    width=1400,
    height=1000,
    showlegend=True,
    legend=dict(
        x=0.02,
        y=0.98,
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="Black",
        borderwidth=1
    )
)

# 保存为 HTML
output_file = "cost_difference_3d_interactive.html"
fig.write_html(output_file)

print(f"\n{'='*60}")
print(f"✓ 交互式 3D 图表已生成: {output_file}")
print(f"{'='*60}")
print(f"  打开方式: 在浏览器中打开 {output_file}")
print(f"  操作说明:")
print(f"    • 拖动鼠标 - 旋转视角")
print(f"    • 滚轮 - 缩放")
print(f"    • 点击图例 - 切换显示/隐藏")
print(f"{'='*60}\n")
