import numpy as np
import matplotlib.pyplot as plt
# compares the twist angle results from all 5 sequences in BSC1
# data from process_twist.py and block_average.py scripts
labels = [
    "Sequence 1", 
    "Sequence 2", 
    "Sequence 3", 
    "Sequence 4", 
    "Sequence 5"
]

# hard coded final averages and std errors from block_average.py outputs
final_averages = [
    35.191,  # Sequence 1
    34.991,  # Sequence 2
    35.189,  # Sequence 3
    35.196,  # Sequence 4
    35.202   # Sequence 5
]

final_std_errors = [
    0.109,   # Sequence 1
    0.044,   # Sequence 2
    0.178,   # Sequence 3
    0.143,   # Sequence 4
    0.104    # Sequence 5
]


print("Plotting comparison of all sequences...")

# creates the figure and axis with higher resolution
fig, ax = plt.subplots(figsize=(10, 6), dpi=100)

# define x positions for bars
x_positions = np.arange(len(labels))

# create a color gradient for visual appeal (not important, just aesthetics)
colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(labels)))

# create the bar plot with enhanced styling, this all just defines how the bars look
bars = ax.bar(x_positions, final_averages, 
               yerr=final_std_errors,
               align='center', 
               alpha=0.85, 
               capsize=8,
               color=colors,
               edgecolor='black',
               linewidth=1.2,
               error_kw={'linewidth': 2, 'ecolor': 'dimgray', 'alpha': 0.8})

# Add value labels on top of bars with standard error
for i, (bar, avg, err) in enumerate(zip(bars, final_averages, final_std_errors)):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + err + 0.05,
            f'{avg:.3f}° ± {err:.3f}°',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_xlabel('DNA Sequence', fontsize=13, fontweight='bold', labelpad=10)
ax.set_ylabel('Overall Average Internal Twist (degrees)', 
              fontsize=13, fontweight='bold', labelpad=10)
ax.set_title('BSC1: Comparison of Average Twist Across 5 DNA Sequences', 
             fontsize=13, fontweight='bold', pad=20)

# Customize x-axis
ax.set_xticks(x_positions)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)

# Customize y-axis
ax.tick_params(axis='y', labelsize=11)
ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=0.8)
ax.set_axisbelow(True)  # Grid behind bars

# Add subtle background color
ax.set_facecolor('#f8f9fa')

# Adjust spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)

# Set y-axis limits with appropriate padding for tight error bars
y_max = max([avg + err for avg, err in zip(final_averages, final_std_errors)])
ax.set_ylim(33.9, y_max + 0.3)  # Adjusted to show meaningful differences

# Calculate overall statistics
overall_mean = np.mean(final_averages)
overall_std = np.std(final_averages, ddof=1)  # Sample standard deviation

# Calculate helical repeat and its error
helical_repeat = 360 / overall_mean

# Error propagation for helical repeat: if h = 360/θ, then δh = (360/θ²) * δθ
# Use standard error of the mean across sequences
overall_sem = overall_std / np.sqrt(len(final_averages))
helical_repeat_error = (360 / (overall_mean**2)) * overall_sem

# Alternative: calculate helical repeat for each sequence and take std dev
helical_repeats = [360 / avg for avg in final_averages]
helical_repeat_std = np.std(helical_repeats, ddof=1)

# Add summary statistics text box
textstr = (f'Overall Statistics:\n'
           f'Mean Twist = {overall_mean:.3f}° ± {overall_std:.3f}°\n'
           f'Helical Repeat = {helical_repeat:.2f} ± {helical_repeat_std:.2f} bp/turn')
props = dict(boxstyle='round,pad=0.8', facecolor='white', edgecolor='gray', 
             linewidth=1.5, alpha=0.95)
ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='right',
        bbox=props, family='monospace')

plt.tight_layout()

# Save with high quality
plt.savefig('twist_all_sequences_comparison.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')

print("Success! Enhanced plot saved to 'twist_all_sequences_comparison.png'")
print(f"\nSummary:")
print(f"Overall Mean Twist: {overall_mean:.3f}° ± {overall_std:.3f}°")
print(f"Helical Repeat: {helical_repeat:.2f} ± {helical_repeat_std:.2f} bp/turn")
print(f"\nIndividual helical repeats:")
for i, (label, hr) in enumerate(zip(labels, helical_repeats)):
    print(f"  {label}: {hr:.2f} bp/turn")

plt.show()
