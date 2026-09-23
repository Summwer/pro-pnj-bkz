"""Plot comparison charts for 4 strategies - generates both q-based and success_rate-based plots"""
import re
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import log2
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not available. Interactive 3D plot will be skipped.")


def parse_log_file_to_df(filepath):
    """Parse log file and extract data, return complete DataFrame
    
    Args:
        filepath: Path to log file
    Returns:
        DataFrame with columns: q, success_rate, G, dsvp, max_beta
    """
    data = []
    
    with open(filepath, 'r') as f:
        blocks = f.read().split('Build LWE  instance')
    
    for block in blocks:
        if not block.strip():
            continue
        
        # Extract q value
        q_match = re.search(r'q=(\d+)', block)
        if not q_match:
            continue
            
        # Extract n value
        n_match = re.search(r'n=(\d+)', block)
        n_val = int(n_match.group(1)) if n_match else 240 # Default to 240 if not found, or skip? Better to default or let it be None
        
        # Extract success_rate (if present)
        sr_match = re.search(r'success_rate\s*=\s*([\d\.eE\-\+]+)', block)
        
        match = re.search(
            r'β_strategy=\[(.*?)\].*dsvp\s*=\s*(\d+).*G\s*=\s*([\d.]+)',
            block, re.DOTALL
        )
        
        if match:
            beta_str, dsvp, G = match.groups()
            max_beta = int(beta_str.split(',')[-1]) if '...' in beta_str else \
                       max([int(x) for x in beta_str.split(',') if x.strip().isdigit()], default=0)
            
            entry = {
                'n': n_val,
                'q': int(q_match.group(1)),
                'G': float(G),
                'dsvp': int(dsvp),
                'max_beta': max_beta,
                'success_rate': float(sr_match.group(1)) if sr_match else None
            }
            
            data.append(entry)
    
    if not data:
        return pd.DataFrame()
    
    return pd.DataFrame(data)


def parse_log_file(filepath, x_axis='q'):
    """Parse log file and extract data (backward compatible)
    
    Args:
        filepath: Path to log file
        x_axis: X-axis type, 'q' or 'success_rate'
    """
    df = parse_log_file_to_df(filepath)
    
    if df.empty:
        return [], [], [], []
    
    df = df.sort_values(x_axis)
    x_values = df[x_axis].tolist()
    
    # Adjust G values if plotting by success_rate
    if x_axis == 'success_rate':
        Gs = [G - log2(sr) for G, sr in zip(df['G'].tolist(), df['success_rate'].tolist())]
    else:
        Gs = df['G'].tolist()
    
    return x_values, Gs, df['dsvp'].tolist(), df['max_beta'].tolist()


def parse_optimized_log_to_df(filepath):
    """Parse optimized_strategy_vary_*.log specifically, return complete DataFrame
    
    Args:
        filepath: Path to log file
    Returns:
        DataFrame with columns: q, success_rate, G, dsvp, max_beta
    """
    data = []
    current = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('n='):
                match = re.search(r'n=(\d+)\s+m=\d+\s+q=(\d+)', line)
                if match:
                    current = {
                        'n': int(match.group(1)),
                        'q': int(match.group(2))
                    }
            
            elif 'succ_prob' in line or 'success_rate' in line:
                match = re.search(r'(?:succ_prob|success_rate)\s*=\s*([\d\.eE\-\+]+)', line)
                if match:
                    val = float(match.group(1))
                    current['success_rate'] = val
                    if val > 0.95:
                        print(f"DEBUG: Found high success rate: {val} for n={current.get('n')}, q={current.get('q')}")
            
            elif 'S(beta,jump,tours):' in line:
                beta_matches = re.findall(r'\((\d+),\d+,\d+\)', line)
                if beta_matches:
                    current['max_beta'] = int(beta_matches[-1])
            
            elif 'pump-{' in line:
                match = re.search(r'pump-\{\d+,(\d+),\d+\}', line)
                if match:
                    current['dsvp'] = int(match.group(1))
            
            elif 'Min Cost =' in line:
                match = re.search(r'Min Cost = ([\d.]+)', line)
                if match and current:
                    current['G'] = float(match.group(1))
                    if 'success_rate' not in current:
                        current['success_rate'] = None
                    data.append(current.copy())
    
    if not data:
        return pd.DataFrame()
    
    return pd.DataFrame(data)


def parse_optimized_log(filepath, x_axis='q'):
    """Parse optimized_strategy_vary_*.log specifically (backward compatible)
    
    Args:
        filepath: Path to log file
        x_axis: X-axis type, 'q' or 'success_rate'
    """
    df = parse_optimized_log_to_df(filepath)
    
    if df.empty:
        return [], [], [], []
    
    df = df.sort_values(x_axis)
    x_values = df[x_axis].tolist()
    
    # Adjust G values if plotting by success_rate
    if x_axis == 'success_rate':
        Gs = [G - log2(sr) for G, sr in zip(df['G'].tolist(), df['success_rate'].tolist())]
    else:
        Gs = df['G'].tolist()
    
    return x_values, Gs, df['dsvp'].tolist(), df['max_beta'].tolist()


def plot_comparison(datasets, output_path, x_axis='q'):
    """Plot comparison charts
    
    Args:
        datasets: List of datasets
        output_path: Output file path
        x_axis: X-axis type, 'q' or 'success_rate'
    """
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'PingFang SC']
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # Color scheme
    g_colors = ['tab:red', 'tab:purple', 'tab:orange', 'tab:pink']
    g_markers = ['o', 'D', 's', 'v']
    g_linestyles = ['-', '-', '-', '-']
    dim_colors = ['tab:blue', 'tab:green', 'tab:brown', 'tab:cyan']
    
    for idx, (name, xs, Gs, dsvps, max_betas) in enumerate(datasets):
        if not xs:
            continue
        
        # Subplot 1: G cost comparison
        ax1.plot(xs, Gs, color=g_colors[idx], marker=g_markers[idx], 
                linestyle=g_linestyles[idx], label=f'{name} - G')
        
        # Subplot 2: Dimension parameter comparison (same color for same dataset, different shapes)
        color = dim_colors[idx]
        ax2.plot(xs, max_betas, color=color, marker='o', 
                linestyle=g_linestyles[idx], label=f'{name} - Max Beta')
        ax2.plot(xs, dsvps, color=color, marker='^', 
                linestyle=g_linestyles[idx], label=f'{name} - dsvp')
        
        # Fill semi-transparent area between max_beta and dsvp
        ax2.fill_between(xs, max_betas, dsvps, color=color, alpha=0.15)
    
    # Configure subplot 1
    ax1.set_ylabel('log2(gate)')
    ax1.set_title('Cost Comparison (G)')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='best')
    
    # Configure subplot 2
    xlabel = 'Success Rate' if x_axis == 'success_rate' else 'q'
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Dimension')
    ax2.set_title('Dimension Comparison (Max Beta, dsvp)')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f"Plot saved to {output_path}")


def plot_q_comparison_by_success_rate(dataframes_dict, output_path):
    """Plot q comparison chart grouped by success rate
    
    Args:
        dataframes_dict: Dictionary {strategy_name: DataFrame}
        output_path: Output file path
    """
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'PingFang SC']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Collect all unique success rate values
    all_success_rates = set()
    for df in dataframes_dict.values():
        if not df.empty and 'success_rate' in df.columns:
            all_success_rates.update(df['success_rate'].dropna().unique())
    
    if not all_success_rates:
        print("Warning: No success rate data found, using default plotting")
        return False
    
    success_rates = sorted(all_success_rates)
    print(f"  Found {len(success_rates)} different success rates: {success_rates}")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Assign colors for each strategy
    strategy_colors = {
        "Two-step(S0)": 'tab:red',
        "Two-step(S_G6K)": 'tab:purple', 
        "Default G6K": 'tab:orange',
        "Two-step(S_OP)": 'tab:pink'
    }
    
    # Assign line styles and markers for each success rate
    sr_linestyles = ['-', '--', '-.', ':']
    sr_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    
    # Plot curves for each strategy and success rate combination
    for strategy_idx, (strategy_name, df) in enumerate(dataframes_dict.items()):
        if df.empty or 'success_rate' not in df.columns:
            continue
        
        color = strategy_colors.get(strategy_name, f'C{strategy_idx}')
        
        for sr_idx, sr in enumerate(success_rates):
            # Filter data for specific success rate
            df_sr = df[df['success_rate'] == sr].sort_values('q')
            
            if df_sr.empty:
                continue
            
            qs = df_sr['q'].tolist()
            qs = [log2(q) for q in qs]
            dsvps = df_sr['dsvp'].tolist()
            max_betas = df_sr['max_beta'].tolist()
            
            # Critical fix: Normalize G values by subtracting log2(sr) when success_rate < 1, regardless of x-axis
            if sr < 1:
                Gs = [G - log2(sr) for G in df_sr['G'].tolist()]
            else:
                Gs = df_sr['G'].tolist()
            
            linestyle = sr_linestyles[sr_idx % len(sr_linestyles)]
            marker = sr_markers[sr_idx % len(sr_markers)]
            
            # Subplot 1: G cost comparison
            label = f'{strategy_name} (sr={sr})'
            ax1.plot(qs, Gs, color=color, linestyle=linestyle, marker=marker,
                    markersize=5, label=label, alpha=0.8)
            
            # Subplot 2: Dimension parameter comparison
            ax2.plot(qs, max_betas, color=color, linestyle=linestyle, marker='o',
                    markersize=4, label=f'{strategy_name} (sr={sr}) - Max Beta', alpha=0.7)
            ax2.plot(qs, dsvps, color=color, linestyle=linestyle, marker='^',
                    markersize=4, label=f'{strategy_name} (sr={sr}) - dsvp', alpha=0.7)
    
    # Configure subplot 1
    ax1.set_ylabel('log2(gate)', fontsize=12)
    ax1.set_title('Cost Comparison by Success Rate (G vs q)', fontsize=14, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    
    # Configure subplot 2
    ax2.set_xlabel('log2(q)', fontsize=12)
    ax2.set_ylabel('Dimension', fontsize=12)
    ax2.set_title('Dimension Comparison by Success Rate', fontsize=14, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # Add vertical lines at q = n^2 and q = n^3 (where n=240)
    n = 240
    log2_n_squared = 2 * log2(n)
    log2_n_cubed = 3 * log2(n)
    ax1.axvline(x=log2_n_squared, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label=f'q = n² (log2(q) = {log2_n_squared:.2f})')
    ax1.axvline(x=log2_n_cubed, color='darkgray', linestyle='-.', linewidth=1.5, alpha=0.7, label=f'q = n³ (log2(q) = {log2_n_cubed:.2f})')
    ax2.axvline(x=log2_n_squared, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.axvline(x=log2_n_cubed, color='darkgray', linestyle='-.', linewidth=1.5, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f"  Plot saved: {output_path}")
    return True


def plot_cost_difference_3d(trivial_df, optimized_df, output_path):
    """Plot 3D surface of cost difference between two strategies
    
    Args:
        trivial_df: DataFrame for trivial strategy
        optimized_df: DataFrame for optimized strategy
        output_path: Output file path
    """
    print(f"\n{'='*60}")
    print("Generating 3D Cost Difference Plot...")
    print(f"{'='*60}")
    
    # Filter for success_rate = 1.0 only (user doesn't want SR < 1)
    # Use >= 0.99 to handle floating point precision
    if 'success_rate' in optimized_df.columns:
        opt_df_filtered = optimized_df[optimized_df['success_rate'] >= 0.99].copy()
        print(f"Filtered for success_rate >= 0.99: {len(opt_df_filtered)} data points")
    else:
        print("Warning: No success_rate column in optimized_df, using all data")
        opt_df_filtered = optimized_df.copy()
    
    # # Previous code that plotted ALL points (kept for reference):
    # print(f"Total optimized data points available: {len(optimized_df)}")
    # opt_df_filtered = optimized_df.copy()
    
    if opt_df_filtered.empty or trivial_df.empty:
        print("Error: Insufficient data for comparison")
        return False
        
    # 2. Merge DataFrames on n and q
    # Ensure 'n' exists in trivial_df (older logs might strictly be 240 if not parsed)
    if 'n' not in trivial_df.columns:
        print("Warning: 'n' column missing in trivial_df, assuming n=240")
        trivial_df['n'] = 240
    
    # Select relevant columns and merge
    t_cols = trivial_df[['n', 'q', 'G']].rename(columns={'G': 'G_trivial'})
    o_cols = opt_df_filtered[['n', 'q', 'G']].rename(columns={'G': 'G_optimized'})
    
    merged_df = pd.merge(t_cols, o_cols, on=['n', 'q'], how='inner')
    print(f"Matched {len(merged_df)} data points for comparison")
    
    if merged_df.empty:
        print("Error: No matching (n, q) pairs found")
        return False
        
    # 3. Calculate difference
    merged_df['diff'] = merged_df['G_trivial'] - merged_df['G_optimized']
    merged_df['log2_q'] = merged_df['q'].apply(log2)
    
    # 4. Find maximum difference point
    max_diff_row = merged_df.loc[merged_df['diff'].idxmax()]
    max_diff_val = max_diff_row['diff']
    max_diff_n = max_diff_row['n']
    max_diff_logq = max_diff_row['log2_q']
    
    print(f"Maximum difference found: {max_diff_val:.4f} at n={max_diff_n}, log2(q)={max_diff_logq:.2f}")
    
    # 5. Create 3D Plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Convert data for plotting
    X = merged_df['n']
    Y = merged_df['log2_q']
    Z = merged_df['diff']
    
    unique_n = merged_df['n'].unique()
    unique_q = merged_df['q'].unique()
    unique_points = len(merged_df)
    
    print(f"Data distribution: {len(unique_n)} unique n values, {len(unique_q)} unique q values")
    
    # Check if we have success_rate variation to create a mesh surface
    if 'success_rate' in merged_df.columns:
        unique_sr = merged_df['success_rate'].unique()
        print(f"Success rate variation: {len(unique_sr)} unique values")
    else:
        unique_sr = []
    
    # Try to create mesh surface using (n, success_rate) as grid if q is constant
    if len(unique_q) < 2 and len(unique_sr) >= 2 and len(unique_n) >= 2:
        print(f"Creating mesh surface using (n, success_rate) grid...")
        
        # Create mesh grid from data
        try:
            # Pivot data to create grid
            pivot_df = merged_df.pivot_table(
                index='success_rate', 
                columns='n', 
                values='diff',
                aggfunc='mean'
            )
            
            # Create mesh grid
            n_vals = np.sort(merged_df['n'].unique())
            sr_vals = np.sort(merged_df['success_rate'].unique())
            
            N_grid, SR_grid = np.meshgrid(n_vals, sr_vals)
            Z_grid = pivot_df.values
            
            # Plot mesh surface (Matlab-style wireframe)
            surf = ax.plot_wireframe(N_grid, SR_grid, Z_grid, 
                                    color='blue', linewidth=1, alpha=0.7,
                                    rstride=1, cstride=1)
            
            # Also add surface with color for better visualization
            surf_colored = ax.plot_surface(N_grid, SR_grid, Z_grid,
                                          cmap='viridis', alpha=0.3, 
                                          edgecolor='blue', linewidth=0.5)
            
            fig.colorbar(surf_colored, ax=ax, shrink=0.5, aspect=5, label='Cost Difference (log2 gate)')
            
            # Update axis labels for this view
            ax.set_xlabel('n (Dimension)', fontsize=10)
            ax.set_ylabel('Success Rate', fontsize=10)
            ax.set_zlabel('Difference (G_trivial - G_optimized)', fontsize=10)
            
            # Mark maximum point
            max_idx = np.unravel_index(np.argmax(Z_grid), Z_grid.shape)
            max_n = N_grid[max_idx]
            max_sr = SR_grid[max_idx]
            max_z = Z_grid[max_idx]
            
            ax.scatter([max_n], [max_sr], [max_z], color='red', s=150, 
                      marker='*', zorder=10, edgecolors='black', linewidths=1)
            
            label_text = f"Max: {max_z:.2f}\n(n={int(max_n)}, SR={max_sr:.1f})"
            ax.text(max_n, max_sr, max_z + (max_z * 0.05), label_text, 
                   color='black', fontweight='bold', fontsize=9)
            
            ax.set_title('Cost Difference: Trivial vs Optimized Strategy\n(Mesh Surface: n × Success Rate)', 
                        fontsize=12, fontweight='bold')
            
        except Exception as e:
            print(f"Warning: Could not create mesh surface ({e}). Using fallback.")
            # Fallback to scatter
            colors = merged_df['success_rate'].values
            scatter = ax.scatter(X, Y, Z, c=colors, cmap='RdYlGn', s=100, alpha=0.8, 
                               depthshade=True, vmin=0, vmax=1, edgecolors='black', linewidths=0.5)
            fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='Success Rate')
            
    elif len(unique_n) < 2 or len(unique_q) < 2:
        print(f"Insufficient grid dimensions. Plotting as 3D scatter/line.")
        
        # Plot as scatter with color-coded success rate
        if 'success_rate' in merged_df.columns:
            colors = merged_df['success_rate'].values
            scatter = ax.scatter(X, Y, Z, c=colors, cmap='RdYlGn', s=100, alpha=1.0, 
                               depthshade=False, vmin=0, vmax=1, edgecolors='black', linewidths=0.5)
            fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='Success Rate')
        else:
            ax.scatter(X, Y, Z, c=Z, cmap='viridis', s=100, alpha=1.0, depthshade=False)
        
        # Connect points with lines
        if unique_points > 1:
            sorted_df = merged_df.sort_values(['n', 'log2_q'])
            X_sorted = sorted_df['n']
            Y_sorted = sorted_df['log2_q']
            Z_sorted = sorted_df['diff']
            ax.plot(X_sorted, Y_sorted, Z_sorted, color='blue', linewidth=2, alpha=0.6)
            
        # Add stems
        for xi, yi, zi in zip(X, Y, Z):
            ax.plot([xi, xi], [yi, yi], [0, zi], 'k--', alpha=0.2, linewidth=0.5)
            
    elif unique_points < 4:
        print(f"Only {unique_points} data point(s) found. Plotting as 3D scatter.")
        
        # Plot as scatter with color-coded success rate
        if 'success_rate' in merged_df.columns:
            colors = merged_df['success_rate'].values
            scatter = ax.scatter(X, Y, Z, c=colors, cmap='RdYlGn', s=100, alpha=1.0, 
                               depthshade=False, vmin=0, vmax=1, edgecolors='black', linewidths=0.5)
            fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='Success Rate')
        else:
            ax.scatter(X, Y, Z, c=Z, cmap='viridis', s=100, alpha=1.0, depthshade=False)
        
        if unique_points > 1:
            sorted_indices = Y.argsort()
            ax.plot(X.iloc[sorted_indices], Y.iloc[sorted_indices], Z.iloc[sorted_indices], 
                    color='blue', label='Difference Curve', linewidth=2)
            
        # Add stems
        for xi, yi, zi in zip(X, Y, Z):
            ax.plot([xi, xi], [yi, yi], [0, zi], 'k--', alpha=0.3)
            
    else:
        # Try to create surface plot
        try:
            # Mesh plot style with edgecolor to look like "wang tu" (net/mesh)
            surf = ax.plot_trisurf(X, Y, Z, cmap='viridis', edgecolor='k', linewidth=0.3, alpha=0.8, antialiased=True)
            
            # Add colorbar only for surface
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Cost Difference (log2 gate)')
        except Exception as e:
            print(f"Warning: Could not create surface plot ({e}). Falling back to scatter.")
            # Fallback to colored scatter
            if 'success_rate' in merged_df.columns:
                colors = merged_df['success_rate'].values
                scatter = ax.scatter(X, Y, Z, c=colors, cmap='RdYlGn', s=80, alpha=0.8, 
                                   depthshade=True, vmin=0, vmax=1, edgecolors='black', linewidths=0.5)
                fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='Success Rate')
            else:
                ax.scatter(X, Y, Z, c=Z, cmap='viridis', s=50, alpha=0.8)
    
    # Annotate maximum point
    # Find max point separately to ensure we mark the global max
    ax.scatter([max_diff_n], [max_diff_logq], [max_diff_val], color='red', s=150, label='Max Difference', marker='*', zorder=10)
    
    # Add text label for max point
    label_text = f"Max: {max_diff_val:.2f}\n(n={int(max_diff_n)}, logq={max_diff_logq:.1f})"
    ax.text(max_diff_n, max_diff_logq, max_diff_val + (max_diff_val * 0.05), label_text, color='black', fontweight='bold')
    

    
    # Labels
    ax.set_xlabel('n (Dimension)')
    ax.set_ylabel('log2(q)')
    ax.set_zlabel('Difference (G_trivial - G_optimized)')
    ax.set_title('Cost Difference: Trivial vs Optimized Strategy (SR=1)')
    
    # View angle (optional, can adjust)
    ax.view_init(elev=30, azim=135)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f"3D Plot saved: {output_path}")
    return True


def load_and_process_data(x_axis):
    """Load and process data for specified x-axis type"""
    # Define data sources
    files = [
        ("Two-step(S0)", f"trivial_strategy_vary_{x_axis}.log", parse_log_file),
        ("Two-step(S_G6K)", f"trivial_strategy_g6k_vary_{x_axis}.log", parse_log_file),
        ("Default G6K", f"default_g6k_vary_{x_axis}.log", parse_log_file),
        ("Two-step(S_OP)", f"optimized_strategy_vary_{x_axis}.log", parse_optimized_log)
    ]
    
    # Load data
    datasets = []
    print(f"\n{'='*60}")
    print(f"Loading {x_axis.upper()} mode data...")
    print(f"{'='*60}")
    
    for name, filename, parser_func in files:
        try:
            xs, Gs, dsvps, max_betas = parser_func(filename, x_axis=x_axis)
            print(f"✓ {name}: Parsed {len(xs)} data points ({filename})")
            datasets.append((name, xs, Gs, dsvps, max_betas))
        except FileNotFoundError:
            print(f"✗ {name}: File not found ({filename})")
        except Exception as e:
            print(f"✗ {name}: Parse error ({filename}): {e}")
    
    # Align data (if plotting by success_rate)
    if x_axis == 'success_rate' and datasets:
        non_empty_datasets = [(name, xs) for name, xs, _, _, _ in datasets if len(xs) > 0]
        
        if non_empty_datasets:
            # Find dataset with minimum data points
            min_dataset_name, min_xs = min(non_empty_datasets, key=lambda x: len(x[1]))
            reference_xs = set(min_xs)
            print(f"\nAligning to {min_dataset_name}'s {len(reference_xs)} success_rate values")
            
            # Align all datasets
            aligned_datasets = []
            for name, xs, Gs, dsvps, max_betas in datasets:
                if len(xs) == 0:
                    aligned_datasets.append((name, [], [], [], []))
                    continue
                
                # Create data dictionary for easy lookup
                data_dict = {x: (G, dsvp, beta) for x, G, dsvp, beta in zip(xs, Gs, dsvps, max_betas)}
                
                # Keep only values in reference_xs
                aligned_xs = []
                aligned_Gs = []
                aligned_dsvps = []
                aligned_max_betas = []
                
                for x in sorted(reference_xs):
                    if x in data_dict:
                        aligned_xs.append(x)
                        G, dsvp, beta = data_dict[x]
                        aligned_Gs.append(G)
                        aligned_dsvps.append(dsvp)
                        aligned_max_betas.append(beta)
                
                print(f"  After alignment {name}: {len(aligned_xs)} data points")
                aligned_datasets.append((name, aligned_xs, aligned_Gs, aligned_dsvps, aligned_max_betas))
            
            datasets = aligned_datasets
    
    return datasets


if __name__ == "__main__":
    print("="*60)
    print("Strategy Comparison Plot Generator")
    print("="*60)
    
    # Define data sources
    files = [
        ("Two-step(S0)", "trivial_strategy_vary_q.log", parse_log_file_to_df),
        ("Two-step(S_G6K)", "trivial_strategy_g6k_vary_q.log", parse_log_file_to_df),
        ("Default G6K", "default_g6k_vary_q.log", parse_log_file_to_df),
        ("Two-step(S_OP)", "optimized_strategy_vary_q.log", parse_optimized_log_to_df)
    ]
    
    # === Generate Q mode plots ===
    print("\n" + "="*60)
    print("Loading Q mode data...")
    print("="*60)
    
    # Load DataFrames
    dataframes_dict = {}
    has_success_rate_data = False
    
    for name, filename, parser_func in files:
        try:
            df = parser_func(filename)
            if not df.empty:
                print(f"✓ {name}: Parsed {len(df)} data points ({filename})")
                dataframes_dict[name] = df
                
                # Check if success_rate data exists
                if 'success_rate' in df.columns and df['success_rate'].notna().any():
                    has_success_rate_data = True
            else:
                print(f"✓ {name}: Parsed 0 data points ({filename})")
        except FileNotFoundError:
            print(f"✗ {name}: File not found ({filename})")
        except Exception as e:
            print(f"✗ {name}: Parse error ({filename}): {e}")
    
    # Choose plotting method based on whether success_rate data exists
    if dataframes_dict:
        if has_success_rate_data:
            print("\nDetected success rate data, using success-rate-grouped plotting")
            success = plot_q_comparison_by_success_rate(
                dataframes_dict, 
                "vary_q_by_success_rate_comparison_plot.png"
            )
            if success:
                print("✓ Q mode plot generated successfully\n")
        else:
            print("\nNo success rate data detected, using traditional plotting")
            datasets = load_and_process_data('q')
            if datasets and any(len(xs) > 0 for _, xs, _, _, _ in datasets):
                plot_comparison(datasets, "vary_q_comparison_plot.png", x_axis='q')
                print("✓ Plot saved: vary_q_comparison_plot.png\n")
    
    # === Generate Success Rate mode plots ===
    datasets = load_and_process_data('success_rate')
    
    if datasets and any(len(xs) > 0 for _, xs, _, _, _ in datasets):
        plot_comparison(datasets, "vary_success_rate_comparison_plot.png", x_axis='success_rate')
        print("✓ Plot saved: vary_success_rate_comparison_plot.png\n")
    
    print("="*60)
    print("Complete! Generated plots:")
    if has_success_rate_data:
        print("  1. vary_q_by_success_rate_comparison_plot.png (Q comparison, grouped by success rate)")
    else:
        print("  1. vary_q_comparison_plot.png (Q comparison)")
    print("  2. vary_success_rate_comparison_plot.png (Success rate comparison)")
    
    # === Generate 3D Cost Difference Plot ===
    # We need "Two-step(S0)" as trivial and "Two-step(S_OP)" as optimized
    if "Two-step(S0)" in dataframes_dict and "Two-step(S_OP)" in dataframes_dict:
        trivial_df = dataframes_dict["Two-step(S0)"]
        optimized_df = dataframes_dict["Two-step(S_OP)"]
        
        plot_cost_difference_3d(
            trivial_df, 
            optimized_df, 
            "cost_difference_3d_plot.png"
        )
        print("  3. cost_difference_3d_plot.png (3D Cost Difference Surface)")
        
        # Also generate interactive HTML version
        if PLOTLY_AVAILABLE:
            plot_cost_difference_3d_interactive(
                trivial_df,
                optimized_df,
                "cost_difference_3d_interactive.html"
            )
            print("  4. cost_difference_3d_interactive.html (Interactive 3D - Open in Browser)")
        
    print("="*60)
