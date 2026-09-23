"""
Standalone script to generate interactive 3D plot
Usage: python plot_3d_interactive.py
Output: cost_difference_3d_interactive.html (open in browser)
"""

import re
import pandas as pd
import numpy as np
from math import log2

try:
    import plotly.graph_objects as go
except ImportError:
    print("ERROR: plotly not installed. Install with: pip install plotly")
    exit(1)


def parse_trivial(filepath):
    data = []
    with open(filepath, 'r') as f:
        blocks = f.read().split('Build LWE  instance')
    
    for block in blocks:
        if not block.strip(): continue
        
        q_match = re.search(r'q=(\d+)', block)
        n_match = re.search(r'n=(\d+)', block)
        match = re.search(r'β_strategy=\[(.*?)\].*dsvp\s*=\s*(\d+).*G\s*=\s*([\d.]+)', block, re.DOTALL)
        
        if q_match and match:
            n_val = int(n_match.group(1)) if n_match else 240
            data.append({'n': n_val, 'q': int(q_match.group(1)), 'G': float(match.group(3))})
    
    return pd.DataFrame(data)


def parse_optimized(filepath):
    data = []
    current = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('n='):
                match = re.search(r'n=(\d+)\s+m=\d+\s+q=(\d+)', line)
                if match:
                    current = {'n': int(match.group(1)), 'q': int(match.group(2))}
            elif 'succ_prob' in line or 'success_rate' in line:
                match = re.search(r'(?:succ_prob|success_rate)\s*=\s*([\d\.eE\-\+]+)', line)
                if match and current:
                    current['success_rate'] = float(match.group(1))
            elif 'Min Cost =' in line:
                if current and 'success_rate' in current:
                    match = re.search(r'Min Cost = ([\d.]+)', line)
                    if match:
                        current['G'] = float(match.group(1))
                        data.append(current.copy())
                current = {}
    
    return pd.DataFrame(data)


# Load data
print("Loading data...")
trivial_df = parse_trivial('trivial_strategy_vary_q.log')
optimized_df = parse_optimized('optimized_strategy_vary_q.log')

# Filter SR >= 0.99
opt_df_filtered = optimized_df[optimized_df['success_rate'] >= 0.99].copy()
print(f"Found {len(opt_df_filtered)} optimized data points with SR >= 0.99")

# Merge
t_cols = trivial_df[['n', 'q', 'G']].rename(columns={'G': 'G_trivial'})
o_cols = opt_df_filtered[['n', 'q', 'G']].rename(columns={'G': 'G_optimized'})
merged_df = pd.merge(t_cols, o_cols, on=['n', 'q'])

print(f"Matched {len(merged_df)} data points")

# Calculate difference
merged_df['diff'] = merged_df['G_trivial'] - merged_df['G_optimized']
merged_df['log2_q'] = merged_df['q'].apply(log2)

X = merged_df['n'].values
Y = merged_df['log2_q'].values
Z = merged_df['diff'].values

max_idx = np.argmax(Z)
print(f"Maximum: {Z[max_idx]:.2f} at n={X[max_idx]}, log2(q)={Y[max_idx]:.2f}")

# Create 3D plot
fig = go.Figure()

# Data points with mesh lines
fig.add_trace(go.Scatter3d(
    x=X, y=Y, z=Z,
    mode='markers+lines',
    marker=dict(
        size=12,
        color=Z,
        colorscale='Viridis',
        showscale=True,
        colorbar=dict(title="Cost<br>Difference"),
        line=dict(color='black', width=1)
    ),
    line=dict(color='blue', width=4),
    name='Cost Difference'
))

# Stems to floor
for xi, yi, zi in zip(X, Y, Z):
    fig.add_trace(go.Scatter3d(
        x=[xi, xi], y=[yi, yi], z=[0, zi],
        mode='lines',
        line=dict(color='lightgray', width=2, dash='dot'),
        showlegend=False,
        hoverinfo='skip'
    ))

# Maximum point
fig.add_trace(go.Scatter3d(
    x=[X[max_idx]],
    y=[Y[max_idx]],
    z=[Z[max_idx]],
    mode='markers+text',
    marker=dict(size=20, color='red', symbol='diamond', line=dict(color='darkred', width=2)),
    text=[f"最大值: {Z[max_idx]:.2f}"],
    textposition="top center",
    name='Maximum',
    textfont=dict(size=14, color='red', family='Arial Black')
))

# Layout
fig.update_layout(
    title=dict(
        text='成本差异：Trivial vs Optimized 策略 (成功率=1.0)<br><sub>用鼠标拖动旋转 | 滚轮缩放</sub>',
        x=0.5,
        font=dict(size=18)
    ),
    scene=dict(
        xaxis_title='n (维度)',
        yaxis_title='log2(q)',
        zaxis_title='差值 (G_trivial - G_optimized)',
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.3)),
        xaxis=dict(backgroundcolor="rgb(240,240,240)", gridcolor="#white"),
        yaxis=dict(backgroundcolor="rgb(240,240,240)", gridcolor="white"),
        zaxis=dict(backgroundcolor="rgb(240,240,240)", gridcolor="white")
    ),
    width=1300,
    height=900,
    showlegend=True,
    legend=dict(x=0.02, y=0.98)
)

# Save
output_file = "cost_difference_3d_interactive.html"
fig.write_html(output_file)

print(f"\n✓ Interactive 3D plot saved: {output_file}")
print(f"✓ Open {output_file} in a web browser to view and interact")
print(f"  - Drag to rotate")
print(f"  - Scroll to zoom")
print(f"  - Click legend items to toggle visibility")
