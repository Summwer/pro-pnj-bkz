\"\"\"Add this function to plot_vary_q.py after the plot_cost_difference_3d function\"\"\"

def plot_cost_difference_3d_interactive(trivial_df, optimized_df, output_path):
    \"\"\"Create interactive 3D plot with plotly (rotatable in browser)
    
    Args:
        trivial_df: DataFrame for trivial strategy
        optimized_df: DataFrame for optimized strategy  
        output_path: Output HTML file path
    \"\"\"
    if not PLOTLY_AVAILABLE:
        print("Plotly not available. Skipping interactive plot.")
        return False
        
    print(f"\\n{'='*60}")
    print("Generating Interactive 3D Plot (Plotly)...")
    print(f"{'='*60}")
    
    # Filter for SR >= 0.99
    if 'success_rate' in optimized_df.columns:
        opt_df_filtered = optimized_df[optimized_df['success_rate'] >= 0.99].copy()
    else:
        opt_df_filtered = optimized_df.copy()
    
    if opt_df_filtered.empty or trivial_df.empty:
        print("Error: Insufficient data")
        return False
    
    # Ensure 'n' exists in trivial_df
    if 'n' not in trivial_df.columns:
        trivial_df['n'] = 240
    
    # Merge
    t_cols = trivial_df[['n', 'q', 'G']].rename(columns={'G': 'G_trivial'})
    o_cols = opt_df_filtered[['n', 'q', 'G']].rename(columns={'G': 'G_optimized'})
    
    merged_df = pd.merge(t_cols, o_cols, on=['n', 'q'], how='inner')
    
    if merged_df.empty:
        print("Error: No matching data")
        return False
    
    # Calculate difference
    merged_df['diff'] = merged_df['G_trivial'] - merged_df['G_optimized']
    merged_df['log2_q'] = merged_df['q'].apply(log2)
    
    # Extract data
    X = merged_df['n'].values
    Y = merged_df['log2_q'].values
    Z = merged_df['diff'].values
    
    # Find max
    max_idx = np.argmax(Z)
    max_n = X[max_idx]
    max_logq = Y[max_idx]
    max_z = Z[max_idx]
    
    print(f"Maximum difference: {max_z:.4f} at n={max_n}, log2(q)={max_logq:.2f}")
    
    # Create figure
    fig = go.Figure()
    
    # Add scatter points with mesh grid lines
    fig.add_trace(go.Scatter3d(
        x=X, y=Y, z=Z,
        mode='markers+lines',
        marker=dict(
            size=10,
            color=Z,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Cost Difference"),
            line=dict(color='black', width=1)
        ),
        line=dict(color='blue', width=3),
        name='Data Points'
    ))
    
    # Add stems to floor (connecting to z=0)
    for xi, yi, zi in zip(X, Y, Z):
        fig.add_trace(go.Scatter3d(
            x=[xi, xi],
            y=[yi, yi],
            z=[0, zi],
            mode='lines',
            line=dict(color='gray', width=2, dash='dash'),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Mark maximum point
    fig.add_trace(go.Scatter3d(
        x=[max_n],
        y=[max_logq],
        z=[max_z],
        mode='markers+text',
        marker=dict(size=18, color='red', symbol='diamond', line=dict(color='darkred', width=2)),
        text=[f"Max: {max_z:.2f}"],
        textposition="top center",
        name='Maximum',
        textfont=dict(size=14, color='red', family='Arial Black')
    ))
    
    # Layout
    fig.update_layout(
        title=dict(
            text='Cost Difference: Trivial vs Optimized Strategy (SR=1.0)<br><sub>Interactive 3D Mesh Plot - Drag to Rotate</sub>',
            x=0.5,
            xanchor='center',
            font=dict(size=16)
        ),
        scene=dict(
            xaxis_title='n (Dimension)',
            yaxis_title='log2(q)',
            zaxis_title='Difference (G_trivial - G_optimized)',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.3)
            ),
            xaxis=dict(backgroundcolor="rgb(230, 230,230)", gridcolor="white", showbackground=True),
            yaxis=dict(backgroundcolor="rgb(230, 230,230)", gridcolor="white", showbackground=True),
            zaxis=dict(backgroundcolor="rgb(230, 230,230)", gridcolor="white", showbackground=True),
        ),
        width=1200,
        height=900,
        showlegend=True
    ))
    
    # Save as HTML
    fig.write_html(output_path)
    print(f"Interactive 3D plot saved: {output_path}")
    print(f"✓ Open '{output_path}' in a web browser to rotate and interact")
    return True
