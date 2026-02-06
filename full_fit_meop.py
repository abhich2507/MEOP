# %%
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np
import ROOT
from array import array

# %%
file_path="/Users/snip/Documents/MEOP/Diameter runs for fitting/1_999_trial.csv"

# df= pd.read_excel(file_path)
df= pd.read_csv(file_path)
name=file_path.split("/")
# name[-1]=name[-1].replace("_trial.xlsx","")
name[-1]=name[-1].replace("_trial.csv","")
diameter=float(name[-1].replace("_","."))
print(diameter)
# df.drop(columns=["Time - Voltage (Formula Result)","Polarization (%) - Voltage (Formula Result)"], inplace=True)
# df.rename(columns={"Time - Voltage (Formula Result).1":"time","Polarization (%) - Voltage (Formula Result).1":"polarization"}, inplace=True)


# %%
# Convert time string to seconds
# Handle both "MM:SS" and "HH:MM:SS" formats
def time_to_seconds(time_str):
    parts = str(time_str).split(':')
    if len(parts) == 2:  # MM:SS format
        return int(parts[0]) * 60 + float(parts[1])
    elif len(parts) == 3:  # HH:MM:SS format
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    else:
        return float(time_str)  # Already in seconds

# Convert all times to seconds
time_seconds = df["time"].apply(time_to_seconds).values

# Detect and fix rollovers (when time resets from 59:59 to 00:00)
# Add 3600 seconds (1 hour) each time we detect a backwards jump
corrected_time = [time_seconds[0]]
cumulative_offset = 0

for i in range(1, len(time_seconds)):
    # If current time is significantly less than previous time, we've had a rollover
    if time_seconds[i] < time_seconds[i-1] - 10:  # -10 threshold to avoid false positives from noise
        cumulative_offset += 3600  # Add 1 hour (60 minutes * 60 seconds)
        print(f"Detected time rollover at index {i}: {time_seconds[i-1]:.1f}s -> {time_seconds[i]:.1f}s (adding {cumulative_offset}s offset)")
    
    corrected_time.append(time_seconds[i] + cumulative_offset)

df["time"] = corrected_time

# Rescale time to start from zero
df["time"] = df["time"] - df["time"].min()
# df=df[df["time"]>200].reset_index(drop=True)
plt.figure(figsize=(10,6))
plt.plot(df["time"], df["polarization"], label="Raw Data", alpha=0.5)
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Polarization")
plt.show()


# %%



# max_index = df.loc[:decay_start_index, "smoothed"].idxmax()
max_pol = df["smoothed"][max_index]
min_index = df["smoothed"].idxmin()
min_pol = df["smoothed"][min_index]
print(f"Max index: {max_index}, Min index: {min_index}")
target_value = max_pol / np.e 
print(f"Target value for relaxation: {target_value:.2f}%")
relaxation_data = df.iloc[max_index:min_index+1]
relaxation_data_full = relaxation_data.copy()

# Filter for points above target value
relaxation_data = relaxation_data[relaxation_data["smoothed"] > target_value]
relaxation_data_full = relaxation_data_full[relaxation_data_full["smoothed"] > min_pol]
max_time = df["time"][max_index]

# Check if we have valid relaxation data after filtering
if len(relaxation_data) == 0:
    print(f"⚠ Warning: No data points above target value (1/e point)")
    print(f"Using full decay range for analysis")
    relaxation_data = relaxation_data_full
    if len(relaxation_data) < 2:
        print(f"Error: Insufficient data points for relaxation analysis")
        relaxation_time = 0
    else:
        relaxation_time = relaxation_data["time"].iloc[-1] - max_time
else:
    relaxation_time = relaxation_data["time"].iloc[-1] - max_time

print(f"Global max index: {global_max_index}, Decay start index: {decay_start_index}, Selected max index: {max_index}")
print(f"Max polarization (before decay): {max_pol:.2f}%")
print(f"Relaxation time: {relaxation_time:.3f} seconds")



# %%

from array import array

time_data = relaxation_data_full["time"].values
#pol_data = relaxation_data_full["smoothed"].values
pol_data = relaxation_data_full["polarization"].values


downsample_factor = 1 # use every nth point to downsample
time_data = time_data[::downsample_factor]
pol_data = pol_data[::downsample_factor]
print(f"Downsampled data by factor of {downsample_factor}")
print(f"Number of points after downsampling: {len(time_data)}")

fit_range_min = 75  # minimum time for fit (seconds)
fit_range_max = 350  # maximum time for fit (seconds)

# Initial parameter guesses (will be set after data is loaded)
# p0_init will be calculated as (first_point - baseline)
# tau_init will use the calculated relaxation time
# c_init will be the last point value

# Apply fit range if specified
if fit_range_min is not None or fit_range_max is not None:
    range_min = fit_range_min if fit_range_min is not None else time_data.min()
    range_max = fit_range_max if fit_range_max is not None else time_data.max()
    mask = (time_data >= range_min) & (time_data <= range_max)
    time_data = time_data[mask]
    pol_data = pol_data[mask]
    print(f"Fitting in range: [{range_min}, {range_max}] seconds")
    print(f"Number of points in fit range: {len(time_data)}")
else:
    range_min = time_data.min()
    range_max = time_data.max()
    print(f"Fitting full data range: [{range_min:.1f}, {range_max:.1f}] seconds")

n_points = len(time_data)

#shift time to start from zero for the fit 
time_shifted = time_data - time_data.min()
x_arr = array('d', time_shifted.tolist())
y_arr = array('d', pol_data.tolist())

print(f"Time range for fit: {time_data.min():.1f} to {time_data.max():.1f} seconds")
print(f"Shifted time range: 0 to {time_shifted.max():.1f} seconds (for fit)")

#create TGraph with shifted time
graph = ROOT.TGraph(n_points, x_arr, y_arr)
graph.SetTitle("Relaxation Fit;Time (s);Polarization (%)")
graph.SetMarkerStyle(20)
graph.SetMarkerSize(0.8)
graph.SetMarkerColor(ROOT.kBlue)

#fit function: P0*exp(-t/tau) + c
# parameters: [0]=P0, [1]=tau, [2]=c
# Using shifted time starting from 0
fit_func = ROOT.TF1("fit_func", "[0]*exp(-x/[1]) + [2]", 0, time_shifted.max())

# Set initial parameter estimates
P0_init = 68  # Initial amplitude (difference between start and end)
tau_init = relaxation_time  # Use calculated relaxation time as initial guess
c_init = pol_data[-1]  # Baseline offset

print(f"\nInitial parameter guesses:")
print(f"  P0_init  = {P0_init:.2f} % (amplitude)")
print(f"  tau_init = {tau_init:.2f} s (relaxation time)")
print(f"  c_init   = {c_init:.2f} % (baseline)\n")

fit_func.SetParameter(0, P0_init)
fit_func.SetParameter(1, tau_init)
fit_func.SetParameter(2, c_init)

# Set parameter limits to prevent overfitting
fit_func.SetParLimits(0, 55, 76)  # P0 must be positive and reasonable
fit_func.SetParLimits(1, 100, 150)  # tau between 20s and 200s
fit_func.SetParLimits(2, -50, 100)  # baseline constraint

# Set parameter names
fit_func.SetParName(0, "P0")
fit_func.SetParName(1, "tau")
fit_func.SetParName(2, "c")

#buildup function initial parameters
P_inf= 70.0  
P_int=3.0
tau_b= 35.

buildup = ROOT.TF1("buildup", "[0] - ([0] - [1]) * exp(-x/[2])", 125, 205)
buildup.SetParameter(0, P_inf)
buildup.SetParameter(1, P_int)
buildup.SetParameter(2, tau_b)

buildup.SetParName(0, "P_inf")
buildup.SetParName(1, "P_int")
buildup.SetParName(2, "tau_b")





# Perform the fit with improved options
# "R" = use function range, "S" = return fit result, "Q" = quiet mode
# "W" = ignore errors in points (for better chi2), "E" = use Minos errors
fit_result = graph.Fit(fit_func, "RS")  # "RS" for basic fit with result

# Extract fit parameters
P0_fit = fit_func.GetParameter(0)
tau_fit = fit_func.GetParameter(1)
c_fit = fit_func.GetParameter(2)

P0_err = fit_func.GetParError(0)
tau_err = fit_func.GetParError(1)
c_err = fit_func.GetParError(2)

# Get correlation matrix
corr_matrix = fit_result.GetCorrelationMatrix()
n_params = 3

# Print results
print("=" * 60)
print("ROOT Fit Results: P0*exp(-t/tau) + c")
print("=" * 60)
print(f"P0  = {P0_fit:.3f} ± {P0_err:.3f} %")
print(f"tau = {tau_fit:.3f} ± {tau_err:.3f} s  (relaxation time)")
print(f"c   = {c_fit:.3f} ± {c_err:.3f} %  (baseline)")
print(f"\nChi2/NDF = {fit_result.Chi2():.2f} / {fit_result.Ndf()}")
chi2_ndf = fit_result.Chi2() / fit_result.Ndf() if fit_result.Ndf() > 0 else 0
print(f"Chi2/NDF ratio = {chi2_ndf:.3f}")

# Quality check
if chi2_ndf < 0.1:
    print(" WARNING: Chi2/NDF < 0.1 suggests possible overfitting!")
    print("   Consider: reducing data points, using simpler function, or checking for systematic errors")
elif chi2_ndf > 10:
    print("WARNING: Chi2/NDF > 10 suggests poor fit quality!")
    print("   Consider: different fit function or checking data quality")

# Print correlation matrix
print("\nParameter Correlation Matrix:")
print("=" * 60)
print(f"{'':>8} {'P0':>12} {'tau':>12} {'c':>12}")
print("-" * 60)
for i in range(n_params):
    param_names = ['P0', 'tau', 'c']
    row = f"{param_names[i]:>8}"
    for j in range(n_params):
        corr_val = corr_matrix[i][j]
        row += f"{corr_val:>12.4f}"
    print(row)
print("=" * 60)
print("Correlation interpretation:")
print("  |ρ| < 0.3  : weak correlation")
print("  0.3 ≤ |ρ| < 0.7 : moderate correlation")
print("  |ρ| ≥ 0.7  : strong correlation (may indicate overfitting)")
print("=" * 60)

# Check for high correlations
high_corr = []
for i in range(n_params):
    for j in range(i+1, n_params):
        corr_val = abs(corr_matrix[i][j])
        if corr_val > 0.7:
            param_names = ['P0', 'tau', 'c']
            high_corr.append(f"{param_names[i]}-{param_names[j]}: {corr_matrix[i][j]:.3f}")

if high_corr:
    print("\n⚠ HIGH CORRELATIONS DETECTED:")
    for corr_str in high_corr:
        print(f"  • {corr_str}")
    print("\nSuggestions to reduce correlation:")
    print("  1. Fix one parameter (e.g., fix baseline 'c' if known)")
    print("  2. Use a narrower fit range")
    print("  3. Increase downsample_factor further")
    print("  4. Consider if baseline is truly needed (try P0*exp(-t/tau) only)")
print("=" * 60)

# Create canvas and plot
c = ROOT.TCanvas("c_fit", "Exponential Decay Fit", 1000, 700)
c.SetGrid()

graph.Draw("AP")
graph.GetXaxis().SetTitle("Time since decay start (s)")
graph.GetYaxis().SetTitle("Polarization (%)")
fit_func.SetLineColor(ROOT.kRed)
fit_func.SetLineWidth(2)
fit_func.SetRange(0, time_shifted.max())  # Limit fit line to fit range only
fit_func.Draw("SAME")

# Add legend with fit parameters
legend = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
legend.SetTextSize(0.03)
legend.AddEntry(graph, "Data", "p")
legend.AddEntry(fit_func, "Fit: P_{0}e^{-t/#tau} + c", "l")
legend.AddEntry(0, f"P_{0} = {P0_fit:.2f} #pm {P0_err:.2f}", "")
legend.AddEntry(0, f"#tau = {tau_fit:.2f} #pm {tau_err:.2f} s", "")
legend.AddEntry(0, f"c = {c_fit:.2f} #pm {c_err:.2f}", "")
legend.AddEntry(0, f"#chi^{2}/NDF = {fit_result.Chi2():.1f}/{fit_result.Ndf()}", "")
legend.Draw()

c.Draw()
c.Update()

# Save plot
diameter_str = str(diameter).replace(".", "_")
output_path = f"/Users/snip/Documents/MEOP/relaxation_plots/root_fit_{diameter_str}.png"
c.SaveAs(output_path)
print(f"\n✓ Plot saved to: {output_path}")

# %%


# Buildup phase fit range
buildup_range_min = 120.0
buildup_range_max = 210.0

# Constant/plateau phase fit range
constant_range_min = 215.0
constant_range_max = 470.0

# Decay phase fit range
decay_range_min = 488.0
decay_range_max = 1050.0

print("="*70)
print("FIT CONFIGURATION")
print("="*70)
print(f"Buildup fit range:  [{buildup_range_min}, {buildup_range_max}] seconds")
print(f"Constant fit range: [{constant_range_min}, {constant_range_max}] seconds")
print(f"Decay fit range:    [{decay_range_min}, {decay_range_max}] seconds")
print("="*70)

# Prepare full data from df
time_full = df["time"].values
pol_full = df["polarization"].values

# 1. BUILDUP FIT 
print("\n" + "="*70)
print("1. BUILDUP PHASE FIT: P_inf - (P_inf - P_int) * exp(-t/tau_b)")
print("="*70)

# Initial parameter guesses for buildup fit
P_inf_init = 58.0   # Plateau polarization (%)
P_int_init = 10.0    # Initial polarization (%)
tau_b_init = 50.0   # Buildup time constant (s)
time_off = 100.0

# dilter data for buildup range
mask_buildup = (time_full >= buildup_range_min) & (time_full <= buildup_range_max)
time_buildup = time_full[mask_buildup]
pol_buildup = pol_full[mask_buildup]

if len(time_buildup) > 0:
    n_buildup = len(time_buildup)
    x_buildup = array('d', time_buildup.tolist())
    y_buildup = array('d', pol_buildup.tolist())
    graph_buildup = ROOT.TGraph(n_buildup, x_buildup, y_buildup)
    graph_buildup.SetMarkerStyle(21)
    graph_buildup.SetMarkerSize(0.8)
    graph_buildup.SetMarkerColor(ROOT.kGreen+2)
    
    # Define buildup function with manual range
    fit_buildup = ROOT.TF1("fit_buildup", "[0] - ([0] - [1]) * exp(-(x-[3]/[2])", buildup_range_min, buildup_range_max)
    fit_buildup.SetParName(0, "P_inf")
    fit_buildup.SetParName(1, "P_int")
    fit_buildup.SetParName(2, "tau_b")
    fit_buildup.SetParName(3,"time_off")

    fit_buildup.SetParameter(0,P_inf_init)
    fit_buildup.SetParameter(1,P_int_init)
    fit_buildup.SetParameter(2,tau_b_init)
    fit_buildup.SetParameter(3,time_off)
    
    # Perform fit
    result_buildup = graph_buildup.Fit(fit_buildup, "RS")
    
    # Extract parameters
    P_inf_fit = fit_buildup.GetParameter(0)
    P_int_fit = fit_buildup.GetParameter(1)
    tau_b_fit = fit_buildup.GetParameter(2)
    
    P_inf_err = fit_buildup.GetParError(0)
    P_int_err = fit_buildup.GetParError(1)
    tau_b_err = fit_buildup.GetParError(2)
    
    print(f"P_inf  = {P_inf_fit:.3f} ± {P_inf_err:.3f} % (plateau)")
    print(f"P_int  = {P_int_fit:.3f} ± {P_int_err:.3f} % (initial)")
    print(f"tau_b  = {tau_b_fit:.3f} ± {tau_b_err:.3f} s (buildup time)")
    print(f"Chi2/NDF = {result_buildup.Chi2():.2f} / {result_buildup.Ndf()} = {result_buildup.Chi2()/result_buildup.Ndf():.3f}")
else:
    print("⚠ No data in buildup range!")
    graph_buildup = None
    fit_buildup = None


print("="*70)
print("2. Constant Fit ")
# Initial parameter guess for constant fit
P_const_init = 70.0  # Constant polarization value (%)

#  constant range
mask_constant = (time_full >= constant_range_min) & (time_full <= constant_range_max)
time_constant = time_full[mask_constant]
pol_constant = pol_full[mask_constant]

if len(time_constant) > 0:
    n_constant = len(time_constant)
    x_constant = array('d', time_constant.tolist())
    y_constant = array('d', pol_constant.tolist())
    graph_constant = ROOT.TGraph(n_constant, x_constant, y_constant)
    graph_constant.SetMarkerStyle(22)
    graph_constant.SetMarkerSize(0.8)
    graph_constant.SetMarkerColor(ROOT.kOrange+7)
    
    #constant function with manual range
    fit_constant = ROOT.TF1("fit_constant", "[0]", constant_range_min, constant_range_max)
    fit_constant.SetParameter(0, P_const_init)
    result_constant = graph_constant.Fit(fit_constant, "RS")
    
    # Extract parameters
    P_const_fit = fit_constant.GetParameter(0)
    P_const_err = fit_constant.GetParError(0)
    
    print(f"P_const = {P_const_fit:.3f} ± {P_const_err:.3f} %")
    print(f"Chi2/NDF = {result_constant.Chi2():.2f} / {result_constant.Ndf()} = {result_constant.Chi2()/result_constant.Ndf():.3f}")
else:
    print("No data in constant range!")
    graph_constant = None
    fit_constant = None

# ==================== 3. DECAY FIT ====================
print("\n" + "="*35+" 3. Decay Fit "+ "="*35 )
# Initial parameter guesses for decay fit
P0_init = 55.0   # Initial amplitude (%)
tau_init = 120.0  # Decay time constant (s)
c_init = 3.0     # Baseline offset (%)
mask_decay = (time_full >= decay_range_min) & (time_full <= decay_range_max)
# Filter data for decay range
mask_decay = (time_full >= decay_range_min) & (time_full <= decay_range_max)
time_decay = time_full[mask_decay]
pol_decay = pol_full[mask_decay]

if len(time_decay) > 0:
    # Create TGraph
    n_decay = len(time_decay)
    x_decay = array('d', time_decay.tolist())
    y_decay = array('d', pol_decay.tolist())
    graph_decay = ROOT.TGraph(n_decay, x_decay, y_decay)
    graph_decay.SetMarkerStyle(20)
    graph_decay.SetMarkerSize(0.8)
    graph_decay.SetMarkerColor(ROOT.kRed)
    
    # Define decay function with manual range
    fit_decay = ROOT.TF1("fit_decay", "[0]*exp(-x/[1]) + [2]",decay_range_min,decay_range_max)
    # Set initial parameter values
    fit_decay.SetParameter(0, P0_init)
    fit_decay.SetParameter(1, tau_init)
    fit_decay.SetParameter(2, c_init)
    
    # Perform fit
    result_decay = graph_decay.Fit(fit_decay, "RS")
    P0_fit = fit_decay.GetParameter(0)
    # Extract parameters
    P0_fit = fit_decay.GetParameter(0)
    tau_fit = fit_decay.GetParameter(1)
    c_fit = fit_decay.GetParameter(2)
    
    P0_err = fit_decay.GetParError(0)
    tau_err = fit_decay.GetParError(1)
    c_err = fit_decay.GetParError(2)
    
    print(f"P0   = {P0_fit:.3f} ± {P0_err:.3f} %")
    print(f"tau  = {tau_fit:.3f} ± {tau_err:.3f} s (decay time)")
    print(f"c    = {c_fit:.3f} ± {c_err:.3f} % (baseline)")
    print(f"Chi2/NDF = {result_decay.Chi2():.2f} / {result_decay.Ndf()} = {result_decay.Chi2()/result_decay.Ndf():.3f}")
else:
    print("⚠ No data in decay range!")
    graph_decay = None
    fit_decay = None
print("\n" + "="*70)
# ==================== COMBINED PLOT ====================
print("\n" + "="*70)
print("Creating combined plot...")
print("="*70)

c_combined = ROOT.TCanvas("c_combined", "Three-Phase Fit", 1200, 800)
c_combined.SetGrid()

# Plot full data
x_full = array('d', time_full.tolist())
y_full = array('d', pol_full.tolist())
graph_full = ROOT.TGraph(len(time_full), x_full, y_full)
graph_full.SetTitle("Three-Phase Fit: Buildup + Constant + Decay;Time (s);Polarization (%)")
graph_full.SetMarkerStyle(20)
graph_full.SetMarkerSize(0.5)
graph_full.SetMarkerColor(ROOT.kBlack)
graph_full.Draw("AP")

# Overlay fits
if fit_buildup:
    fit_buildup.SetLineColor(ROOT.kGreen+2)
    fit_buildup.SetLineWidth(3)
    fit_buildup.Draw("SAME")
    fit_constant.SetLineColor(ROOT.kOrange+7)
if fit_constant:
    fit_constant.SetLineColor(ROOT.kOrange+7)
    fit_constant.SetLineWidth(3)
    fit_constant.Draw("SAME")
    fit_decay.SetLineColor(ROOT.kRed)
if fit_decay:
    fit_decay.SetLineColor(ROOT.kRed)
    fit_decay.SetLineWidth(3)
    fit_decay.Draw("SAME")

# Add legend
legend = ROOT.TLegend(0.15, 0.70, 0.40, 0.88)
legend.SetTextSize(0.025)
legend.AddEntry(graph_full, "Data", "p")
if fit_buildup:
    legend.AddEntry(fit_buildup, "Buildup Fit", "l")
if fit_constant:
    legend.AddEntry(fit_constant, "Constant Fit", "l")
if fit_decay:
    legend.AddEntry(fit_decay, "Decay Fit", "l")
legend.Draw()

# Add fit parameters text box on the right side
param_box = ROOT.TPaveText(0.62, 0.15, 0.88, 0.88, "NDC")
param_box.SetFillColor(ROOT.kWhite)
param_box.SetBorderSize(1)
param_box.SetTextAlign(12)  # Left align
param_box.SetTextSize(0.022)
param_box.SetTextFont(42)

# Section 1: Buildup Parameters
if fit_buildup:
    param_box.AddText("#bf{1. Buildup Phase}")
    param_box.AddText("function: P_{#infty} - (P_{#infty} - P_{int})e^{-t/#tau_{b}}")
    param_box.AddText(" ")
    param_box.AddText(f"P_{{#infty}} = {P_inf_fit:.2f} #pm {P_inf_err:.2f} %")
    param_box.AddText(f"P_{{int}} = {P_int_fit:.2f} #pm {P_int_err:.2f} %")
    param_box.AddText(f"#tau_{{b}} = {tau_b_fit:.2f} #pm {tau_b_err:.2f} s")
    param_box.AddText(f"#chi^{{2}}/NDF = {result_buildup.Chi2():.2f}/{result_buildup.Ndf()}")
    param_box.AddText(" ")

# Section 2: Constant Parameters
if fit_constant:
    param_box.AddText("#bf{2. Constant Phase}")
    param_box.AddText("function: P = constant")
    param_box.AddText(" ")
    param_box.AddText(f"P_{{const}} = {P_const_fit:.2f} #pm {P_const_err:.2f} %")
    param_box.AddText(f"#chi^{{2}}/NDF = {result_constant.Chi2():.2f}/{result_constant.Ndf()}")
    param_box.AddText(" ")

# Section 3: Decay Parameters
if fit_decay:
    param_box.AddText("#bf{3. Decay Phase}")
    param_box.AddText("function: P_{0}e^{-t/#tau} + c")
    param_box.AddText(" ")
    param_box.AddText(f"P_{{0}} = {P0_fit:.2f} #pm {P0_err:.2f} %")
    param_box.AddText(f"#tau = {tau_fit:.2f} #pm {tau_err:.2f} s")
    param_box.AddText(f"c = {c_fit:.2f} #pm {c_err:.2f} %")
    param_box.AddText(f"#chi^{{2}}/NDF = {result_decay.Chi2():.2f}/{result_decay.Ndf()}")

param_box.Draw()

c_combined.Update()
c_combined.Draw()

# Save plot
output_path = f"/Users/snip/Documents/MEOP/relaxation_plots/three_phase_fit.png"
c_combined.SaveAs(output_path)
print(f"\n✓ Combined plot saved to: {output_path}")

# %%
relaxation_data_full


# %%
from array import array

time_data = relaxation_data_full["time"].values
pol_data = relaxation_data_full["smoothed"].values
time_shifted = time_data - time_data.min()
x_arr = array('d', time_shifted.tolist())
y_arr = array('d', pol_data.tolist())

linear_fit = ROOT.TF1("linear_fit", "log([0]) - x/[1]", 0, time_shifted.max())


log_y_arr = array('d', np.log(pol_data).tolist())
graph_log = ROOT.TGraph(n_points, x_arr, log_y_arr)

linear_fit.SetParameter(0, P0_fit)
linear_fit.SetParameter(1, tau_fit)
linear_fit_result = graph_log.Fit(linear_fit, "RS")

graph_log.SetTitle("Log-Scale Relaxation Fit;Time (s);ln(Polarization (%))")
graph_log.SetMarkerStyle(20)
graph_log.SetMarkerSize(0.8)
graph_log.SetMarkerColor(ROOT.kBlue)    

log_canvas= ROOT.TCanvas("log_canvas", "Log-Scale Relaxation Fit", 1000, 700)
log_canvas.SetGrid()
graph_log.Draw("AP")
graph_log.GetXaxis().SetTitle("Time since decay start (s)")
graph_log.GetYaxis().SetTitle("ln(Polarization (%))")   
log_canvas.Draw()

linear_fit.SetLineColor(ROOT.kRed)
linear_fit.Draw("SAME")

# %%


# %%
plt.figure(figsize=(10,6))
plt.plot(df["time"], df["polarization"], label="Raw Data", alpha=0.5)
plt.plot(df["time"], df["smoothed"], label="Smoothed Data", color='red')
plt.xlabel("Time (s)")
plt.ylabel("Polarization")
plt.title("Polarization vs Time with Smoothing")
plt.legend()
plt.show()

# %%


# %%
import os
from pathlib import Path

# Directory containing all the Excel files
data_dir = "/Users/snip/Documents/MEOP/Diameter runs for fitting"
output_dir = "/Users/snip/Documents/MEOP/relaxation_plots"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get all Excel files in the directory
excel_files = list(Path(data_dir).glob("*.xlsx"))
print(f"Found {len(excel_files)} Excel files to process\n")

# Process each file
results = []
for file_path in excel_files:
    try:
        print(f"Processing: {file_path.name}")
        
        # Read the data
        df = pd.read_excel(str(file_path))
        
        # Extract diameter from filename
        name_parts = file_path.stem.split("_")
        diameter = float(file_path.stem.replace("_trial", "").replace("_", "."))
        
        # Drop and rename columns
        df.drop(columns=["Time - Voltage (Formula Result)", "Polarization (%) - Voltage (Formula Result)"], inplace=True)
        df.rename(columns={"Time - Voltage (Formula Result).1": "time", 
                          "Polarization (%) - Voltage (Formula Result).1": "polarization"}, inplace=True)
        
        # Convert time to seconds
        def time_to_seconds(time_str):
            parts = str(time_str).split(':')
            if len(parts) == 2:
                return int(parts[0]) * 60 + float(parts[1])
            elif len(parts) == 3:
                return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
            else:
                return float(time_str)
        
        df["time"] = df["time"].apply(time_to_seconds)
        df["time"] = df["time"] - df["time"].min()  # Start from zero
        
        # Apply smoothing
        df["smoothed"] = savgol_filter(df["polarization"], 60, 3)
        df["derivative"] = np.gradient(df["smoothed"])
        
        # Find max and decay start
        global_max_index = df["smoothed"].idxmax()
        decay_start_index = global_max_index
        global_max_pol = df["smoothed"][global_max_index]
        threshold_10_percent = global_max_pol * 0.9  # 10% decrement threshold
        
        for i in range(global_max_index, len(df)):
            if i + 5 < len(df):
                if all(df["derivative"].iloc[i:i+5] < 0) and df["smoothed"].iloc[i] <= threshold_10_percent:
                    decay_start_index = i
                    break
        
        max_index = df.loc[:decay_start_index, "smoothed"].idxmax()
        max_pol = df["smoothed"][max_index]
        min_index = df["smoothed"].idxmin()
        min_pol = df["smoothed"][min_index]
        
        # Calculate relaxation time
        target_value = max_pol / np.e
        relaxation_data = df.iloc[max_index:min_index+1]
        relaxation_data = relaxation_data[relaxation_data["smoothed"] > target_value]
        max_time = df["time"][max_index]
        
        # Check if we have valid relaxation data
        if len(relaxation_data) == 0:
            print(f"  ⚠ Warning: No data points above target value (1/e point)")
            print(f"  Using full decay range for analysis")
            relaxation_data = df.iloc[max_index:min_index+1]
            if len(relaxation_data) < 2:
                raise ValueError("Insufficient data points for relaxation analysis")
        
        relaxation_time = relaxation_data["time"].iloc[-1] - max_time
        
        # Store results
        results.append({
            "filename": file_path.name,
            "diameter": diameter,
            "relaxation_time": relaxation_time,
            "max_pol": max_pol
        })
        
        # Create and save plot
        plt.figure(figsize=(12, 7))
        plt.plot(df["time"], df["polarization"], label="Raw Data", alpha=0.5)
        plt.plot(df["time"], df["smoothed"], label="Smoothed Data", color='red', linewidth=2)
        
        # Only plot relaxation data if we have valid points above target
        if len(relaxation_data[relaxation_data["smoothed"] > target_value]) > 0:
            plt.plot(relaxation_data["time"], relaxation_data["smoothed"], label="Relaxation Data", color='green', linewidth=2)
        
        # Add vertical lines
        plt.vlines(df["time"][max_index], ymin=0, ymax=max_pol, colors='orange', linestyles='--', 
                  label='Max Polarization Point', linewidth=2)
        target_end_time = relaxation_data["time"].iloc[-1]
        plt.vlines(target_end_time, ymin=0, ymax=target_value, colors='cyan', linestyles='--', 
                  label='Target Value End (1/e point)', linewidth=2)
        plt.hlines(target_value, xmin=df["time"][max_index], xmax=target_end_time, 
                  colors='cyan', linestyles=':', alpha=0.5)
        
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Polarization (%)", fontsize=12)
        plt.title(f"Relaxation Time Analysis\nDiameter: {diameter} inch | τ = {relaxation_time:.3f}s", fontsize=14)
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        
        # Save with original filename
        output_filename = file_path.stem + "_relaxation.png"
        output_path = os.path.join(output_dir, output_filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ Diameter: {diameter} inch | Relaxation time: {relaxation_time:.3f}s")
        print(f"  ✓ Saved: {output_filename}\n")
        
    except Exception as e:
        print(f"  ✗ Error processing {file_path.name}: {str(e)}\n")

# Summary
print(f"\n{'='*60}")
print(f"Processing complete!")
print(f"Processed {len(results)} files successfully")
print(f"Plots saved to: {output_dir}")
print(f"{'='*60}\n")

# Display summary table
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('diameter')
print("\nSummary of Results:")
print(results_df.to_string(index=False))
results_df


