# Time-Series Generation

Here we consider scenarios where the tracking system results in a noisy measurement. To do this, we introduced probabilistic noise into the tracking mechanism, where each frame had a probability \( p \) of experiencing an error—representing transient disruptions in tracking. As \( p \) increased, these disturbances became more frequent. When a failure occurred, a measurement error was sampled from a zero-mean Gaussian distribution with variance 25 and added to the agents’ observed positions.

We evaluated how varying levels of detection noise (parameterized by \( p \)) affected system performance. In particular, we measured each agent’s **tracking error** \( E \) relative to its reference trajectory at steady-state (denoted as the norm of the error of each agent with respect to its corresponding reference trajectory) and counted the number of **collision avoidance maneuvers** \( N_A \) that were triggered. 

For each value of \( p \), we ran **100 independent trials**, where each trial simulated the system from random initial conditions over **100 time steps**. This simulation length was chosen because the unperturbed system typically stabilizes within the first 10 steps.


# Time-Series Visualization

This MATLAB script visualizes the evolution of a variable over time by plotting the **mean trajectory** across multiple trials along with a shaded region representing the **standard deviation**. This type of visualization is commonly used in experimental and simulation-based research to show variability across repeated runs or subjects.

---

## 🧠 Methodology

### 1. **Data Loading**
The script loads precomputed simulation or experimental data from a `.mat` file, which contains:
- A time vector (e.g., `p`) spanning the analysis window.
- Multiple realizations of a signal (e.g., population count, error metric) across trials.

---

### 2. **Data Preparation**
- Signals are reshaped or transposed to ensure each trial corresponds to a row.
- The number of trials and timepoints are inferred automatically.
- For each timepoint, the **mean** and **standard deviation** across all trials are computed.

---

### 3. **Plotting Mean ± Standard Deviation**
- The script uses `fill` to create a **shaded region** showing the envelope of ±1 standard deviation around the mean.
- The **mean signal** is plotted as a bold blue curve.
- The visual style is clean and publication-ready, with white backgrounds, LaTeX-formatted labels, and consistent axis settings.

---

### 4. **Plot Customization**
- Axis labels and legends are rendered using LaTeX for professional mathematical typesetting.
- Ticks, font sizes, and figure dimensions are tuned for clarity and export.

---

## 📌 Output

The result is two high-quality plots:
1. **$N_{\mathrm{A}}$ vs. time ($p$)** — showing the average trajectory of some count variable.
2. **$E$ vs. time ($p$)** — showing the average trajectory of an error or energy variable.
