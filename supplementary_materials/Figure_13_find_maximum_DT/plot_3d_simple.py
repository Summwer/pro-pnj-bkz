#!/usr/bin/env python3
"""
Self-contained interactive 3D plot generator.
Uses plotly to create a rotatable 3D surface plot in the browser.
"""

import re
from math import log2, log
from collections import defaultdict
import plotly.graph_objects as go


def parse_logq_from_dim_line(block):
    """Extract log2(q) from a 'dim = ..., dvol = ...' block.
    Since dvol = n * log2(q), log2(q) = dvol / n.
    Returns log2q or None if parsing fails.
    """
    n_match = re.search(r'n=(\d+)', block)
    dvol_match = re.search(r'dvol\s*=\s*([\d.]+)', block)
    if n_match and dvol_match:
        n = int(n_match.group(1))
        dvol = float(dvol_match.group(1))
        if n > 0:
            log2q = dvol / n
            return log2q
    return None


def parse_data():

    # Parse trivial strategy log file
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
                q = int(q_match.group(1))
                # Compute log2(q) exactly from integer q
                log2q_from_q = round(log2(q))
                log2q = log2q_from_q
                trivial_data.append({
                    'n': n,
                    'q': q,
                    'log2q': log2q,
                    'G': float(g_match.group(1))
                })
    
    # Parse optimized strategy log file
    optimized_data = []
    current = {}
    with open('optimized_strategy_vary_q.log', 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        
        if line.startswith('n = '):
            m = re.search(r'n\s*=\s*(\d+)\s+m\s*=\s*(\d+)\s+log2\(q\)\s*=\s*(\d+)', line)
            if m:
                n_val = int(m.group(1))
                log2q_val = int(m.group(3))
                current = {
                    'n': n_val,
                    'q': pow(2, log2q_val),
                    'log2q': log2q_val
                }
        elif ('succ_prob' in line or 'success_rate' in line) and current:
            m = re.search(r'(?:succ_prob|success_rate)\s*=\s*([\d\.eE\-\+]+)', line)
            if m:
                current['sr'] = float(m.group(1))
        elif 'Min Cost =' in line:
            if current and 'sr' in current:
                m = re.search(r'Min Cost = ([\d.]+)', line)
                if m:
                    current['G'] = float(m.group(1))
                    optimized_data.append(current.copy())
            current = {}
        i += 1
    
    return trivial_data, optimized_data

print("Loading data...")
trivial_list, optimized_list = parse_data()
opt_filtered = [d for d in optimized_list if d.get('sr', 0) >= 0.99]
print(f"Found {len(opt_filtered)} optimized strategy data points (SR >= 0.99)")

# Merge data and build the full grid, filling missing entries with zeros
# Collect unique n and log2(q) values from the trivial data set
all_n = sorted(set([d['n'] for d in trivial_list if d['n'] >= 450]))
all_q = sorted(set([d['log2q'] for d in trivial_list if d['log2q'] >= 30]))

print(f"Grid dimensions: {len(all_n)} n values × {len(all_q)} q values")
print(f"Display range: n >= 250, log2(q) >= 20")

# Create lookup tables for quick access
trivial_dict = {(d['n'], d['log2q']): d['G'] for d in trivial_list}
optimized_dict = {(d['n'], d['log2q']): d['G'] for d in opt_filtered}

# Build the full grid and compute differences where data is available
complete_grid = []
for n_val in all_n:
    for q_val in all_q:
        key = (n_val, q_val)
        
        if key in trivial_dict and key in optimized_dict:
            diff = trivial_dict[key] - optimized_dict[key]
            complete_grid.append({
                'n': n_val,
                'q': q_val,
                'diff': diff,
                'has_data': True
            })
        else:
            complete_grid.append({
                'n': n_val,
                'q': q_val,
                'diff': 0.0,
                'has_data': False
            })

print(f"Total grid points: {len(complete_grid)} (including {sum(1 for d in complete_grid if not d['has_data'])} fill points)")
print(f"Actual data points: {sum(1 for d in complete_grid if d['has_data'])}")

# Extract coordinates from the complete grid
X = [d['n'] for d in complete_grid]
Y = [d['q'] for d in complete_grid]
Z = [d['diff'] for d in complete_grid]
has_data = [d['has_data'] for d in complete_grid]

# Find maximum (only in real data)
real_data_points = [(X[i], Y[i], Z[i]) for i in range(len(X)) if has_data[i]]
if real_data_points:
    max_point = max(real_data_points, key=lambda p: p[2])
    print(f"Maximum difference: {max_point[2]:.2f} at n={max_point[0]}, log2(q)={max_point[1]:.2f}")

# Export data to Excel
import pandas as pd

# Print debug information about parsed data
print(f"\nDebug Info:")
print(f"  trivial_dict has {len(trivial_dict)} entries")
print(f"  optimized_dict has {len(optimized_dict)} entries")
if trivial_dict:
    sample_key = list(trivial_dict.keys())[0]
    print(f"  Sample trivial entry: {sample_key} -> {trivial_dict[sample_key]}")
if optimized_dict:
    sample_key = list(optimized_dict.keys())[0]
    print(f"  Sample optimized entry: {sample_key} -> {optimized_dict[sample_key]}")

# Collect rows for export with original q restored from log2(q)
export_data = []
has_data_count = 0
for d in complete_grid:
    n_val = d['n']
    q_val = d['q']
    key = (n_val, q_val)
    
    if d['has_data']:
        s0 = trivial_dict.get(key)
        sop = optimized_dict.get(key)
        diff = d['diff']
        has_data_count += 1
    else:
        s0 = None
        sop = None
        diff = None
    
    actual_q = 2 ** q_val
    export_data.append({
        'n': n_val,
        'q': actual_q,
        'log2(q)': q_val,
        'S0 (Trivial)': s0,
        'Sop (Optimized)': sop,
        'Difference (S0 - Sop)': diff
    })

print(f"  Rows with actual data: {has_data_count}")

# Create DataFrame and export
df = pd.DataFrame(export_data)
# Sort by n and then by q
df = df.sort_values(['n', 'q'])
excel_file = 'cost_comparison_data.xlsx'
df.to_excel(excel_file, index=False, sheet_name='Cost Comparison')
print(f"\n✓ Excel file generated: {excel_file}")
print(f"  Columns: n, q, log2(q), S0 (Trivial), Sop (Optimized), Difference (S0 - Sop)")
print(f"  Total rows: {len(df)}")


# Create plotly chart
fig = go.Figure()

# Reshape data into a 2D grid for the surface plot
import numpy as np

n_vals = sorted(set(X))
q_vals = sorted(set(Y))

# Create meshgrid
Z_grid = np.zeros((len(n_vals), len(q_vals)))
for i in range(len(X)):
    n_idx = n_vals.index(X[i])
    q_idx = q_vals.index(Y[i])
    Z_grid[n_idx, q_idx] = Z[i]

# Add colored surface
fig.add_trace(go.Surface(
    x=n_vals,
    y=q_vals,
    z=Z_grid.T,  # Transpose to match orientation
    colorscale=[
        [0.0, 'rgb(0, 0, 255)'],      # Blue (low values)
        [0.25, 'rgb(0, 191, 255)'],   # Cyan
        [0.5, 'rgb(0, 255, 127)'],    # Green-cyan
        [0.75, 'rgb(173, 255, 47)'],  # Yellow-green
        [1.0, 'rgb(255, 255, 0)']     # Yellow (high values)
    ],
    showscale=True,
    colorbar=dict(title="Cost Difference", thickness=20, len=0.7),
    opacity=0.5,  # 50% transparency
    name='Surface',
    hovertemplate='n: %{x}<br>log2(q): %{y}<br>Diff: %{z:.2f}<extra></extra>'
))

# Create proper mesh grid lines
# Group data by n and q values
by_n = defaultdict(list)
by_q = defaultdict(list)

for i, (x, y, z) in enumerate(zip(X, Y, Z)):
    by_n[x].append((y, z, i))  # For fixed n, varying q
    by_q[y].append((x, z, i))  # For fixed q, varying n

# Draw lines along q direction (fixed n) - uniform black color
for n_val, points in by_n.items():
    if len(points) > 1:
        # Sort by q (which is y in our plot)
        points_sorted = sorted(points, key=lambda p: p[0])
        
        # Draw each segment with uniform black color
        for i in range(len(points_sorted) - 1):
            y_seg = [points_sorted[i][0], points_sorted[i+1][0]]
            z_seg = [points_sorted[i][1], points_sorted[i+1][1]]
            x_seg = [n_val, n_val]
            
            fig.add_trace(go.Scatter3d(
                x=x_seg, y=y_seg, z=z_seg,
                mode='lines',
                line=dict(color='rgb(0, 0, 0)', width=2),  # Black
                showlegend=False,
                hoverinfo='skip'
            ))

# Draw lines along n direction (fixed q) - uniform black color
for q_val, points in by_q.items():
    if len(points) > 1:
        # Sort by n (which is x in our plot)
        points_sorted = sorted(points, key=lambda p: p[0])
        
        # Draw each segment with uniform black color
        for i in range(len(points_sorted) - 1):
            x_seg = [points_sorted[i][0], points_sorted[i+1][0]]
            z_seg = [points_sorted[i][1], points_sorted[i+1][1]]
            y_seg = [q_val, q_val]
            
            fig.add_trace(go.Scatter3d(
                x=x_seg, y=y_seg, z=z_seg,
                mode='lines',
                line=dict(color='rgb(0, 0, 0)', width=2),  # Black
                showlegend=False,
                hoverinfo='skip'
            ))

# Use a clean surface style without explicit data markers.
# Grid lines already visualize the data structure.

# Mark the maximum point for each fixed n using actual data values
for n_val, points in by_n.items():
    valid_points = [(y, z, i) for (y, z, i) in points if has_data[i]]
    if not valid_points:
        continue
    max_q, max_z, _ = max(valid_points, key=lambda p: p[1])
    fig.add_trace(go.Scatter3d(
        x=[n_val],
        y=[max_q],
        z=[max_z],
        mode='markers',
        marker=dict(
            size=4,
            color='rgb(255, 0, 0)',
            symbol='circle',
            opacity=0.9
        ),
        showlegend=False,
        hovertemplate='n: %{x}<br>log2(q): %{y}<br>Max for n: %{z:.2f}<extra></extra>'
    ))

# Add parametric relationship curve: diff = 0.01887 * n^(-0.98)
curve_n_vals = np.linspace(min(n_vals), max(n_vals), 100)
curve_diff_vals = 0.01887 * curve_n_vals  -0.98
curve_q_vals =   0.198 * curve_n_vals**0.832

fig.add_trace(go.Scatter3d(
    x=curve_n_vals,
    y=curve_q_vals,
    z=curve_diff_vals,
    mode='lines',
    line=dict(color='rgb(255, 128, 0)', width=8),
    name=r'$\Delta(\log_2 T) = 0.01887\; n - 0.98,\; \log_2 q = 0.198\; n^{0.832}$',
    hovertemplate='n: %{x:.2f}<br>diff: %{z:.4f}<extra></extra>'
))

# Mark the overall maximum point with a distinct label
X_real = [X[i] for i in range(len(X)) if has_data[i]]
Y_real = [Y[i] for i in range(len(Y)) if has_data[i]]
Z_real = [Z[i] for i in range(len(Z)) if has_data[i]]

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
            size=7,  # Slightly larger for visibility
            color='rgb(255, 0, 0)',  
            line=dict(color='rgb(0, 0, 0)', width=2),
            opacity=1.0
        ),
        text=[f"Max: {max_z:.2f}"],
        textposition="top center",
        name='Maximum',
        textfont=dict(size=10, color='rgb(0, 0, 0)', family='Arial'),
        hovertemplate='n: %{x}<br>log2(q): %{y}<br>Max for n: %{z:.2f}<extra></extra>'
    ))

# Layout settings
fig.update_layout(
    title=dict(
        text='<b>Cost Difference: Trivial vs Optimized Strategy</b><br><sub>Drag to rotate | Scroll to zoom</sub>',
        x=0.5,
        xanchor='center',
        font=dict(size=20)
    ),
    scene=dict(
        xaxis_title='<b> n </b>',
        yaxis_title='<b>log₂q</b>',
        zaxis_title='<b> log₂(G(S0)) - log₂(G(Sop))</b>',
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.3)
        ),
        aspectmode='manual',  # Use manual aspect ratio
        aspectratio=dict(x=1, y=1, z=0.5),  # Equal spacing for x and y, smaller for z
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
    height=1400,
    showlegend=True,
    legend=dict(
        x=0.2,
        y=0.7,
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="Black",
        borderwidth=1
    )
)

# Save as HTML
output_file = "cost_difference_3d_interactive.html"
fig.write_html(output_file, include_mathjax='cdn')

print(f"\n{'='*60}")
print(f"✓ Interactive 3D plot generated: {output_file}")
print(f"{'='*60}")
print(f"  How to view: Open {output_file} in a web browser")
print(f"  Controls:")
print(f"    • Drag mouse - Rotate view")
print(f"    • Scroll wheel - Zoom")
print(f"    • Click legend - Toggle visibility")
print(f"{'='*60}\n")
