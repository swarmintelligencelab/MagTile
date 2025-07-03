# 🧭 Discrete Trajectory Mapping and Visualization

This project provides a MATLAB-based workflow to map a real-world continuous trajectory onto a discrete grid and visualize both the reference (idealized) trajectory and the actual experimental maneuver. The goal is to compare ideal grid-following behavior with real movement in a structured, interpretable way.

---

## 🧠 Methodology

### 1. **Data Extraction**
A segment of raw trajectory data is extracted from a larger dataset. This represents the real movement of an agent in a 2D space.

---

### 2. **Grid Construction**
A 15×15 spatial grid is defined over the domain of interest. To avoid edge effects, the usable "interior" portion of the grid spans indices 3 to 13. This grid discretizes the continuous space into a finite set of navigable cells.

---

### 3. **Trajectory Discretization**
Each point in the continuous trajectory is mapped to its corresponding grid cell based on position. The result is a coarse, symbolic trajectory represented by grid indices.

---

### 4. **Reference Path Generation**
To make the trajectory suitable for grid-constrained systems (like robotic planners), the path is resolved into a sequence of axis-aligned, single-cell steps. This eliminates diagonal or multi-cell jumps, ensuring feasibility in practice.

---

### 5. **Verification**
A validation step ensures that the resolved reference trajectory only contains valid moves: one grid cell at a time in either the row or column direction, with no diagonals.

---

### 6. **Visualization**
Both the idealized (resolved) and real (continuous) trajectories are plotted together:
- The **reference trajectory** is shown as an orange, grid-snapping path with square markers.
- The **real trajectory** is displayed as a smooth, light blue curve.
- Start and end positions are clearly marked using distinct colors and shapes.
- Grid points are plotted in the background for spatial reference.

Custom axis ticks and labels using LaTeX formatting allow intuitive interpretation of spatial extents like ±xₘₐₓ and ±yₘₐₓ.

---

## 📌 Applications

- Visual debugging of motion planners
- Analyzing deviations between planned and executed paths
- Preparing datasets for grid-based control or reinforcement learning
- Robotic navigation or multi-agent behavior benchmarking
