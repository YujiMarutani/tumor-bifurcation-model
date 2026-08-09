# Treatment Timing, Metabolic Modulation, and Bifurcation-Proximity Control in a Coupled Tumor–Immune–Drug Model (v1.0.1)

**Author:** Yuji Marutani (O'Valley Salon of Meta Design / OAOL Research)  
**Release Version:** v1.0.1 (Computational Consistency Correction)  
**DOI:** [10.5281/zenodo.21864631](https://doi.org/10.5281/zenodo.21864631)  
**License:** Creative Commons Attribution 4.0 International (CC-BY-4.0)

---

## Key Canonical Values (v1.0.1 Synchronization)

This repository provides the canonical code and manuscript sources for the v1.0.1 release, enforcing numerical consistency across implementations and documentation:

| Metric / Parameter | Value | Definition / Notes |
| :--- | :--- | :--- |
| `I_amp` | `0.45` | Fixed drug input amplitude |
| $\Lambda_{\text{basal}}$ | `0.010192` | Basal effective control ratio |
| **Frozen Basal $C_+^*$** | **`989.81`** | **1D frozen-state basal equilibrium** |
| **Untreated $C(40)$** | **`983.46`** | **Actual 5D trajectory value at $t=40\text{d}$** |
| $\Lambda_{\max}$ | `0.667661` | Peak modeled control ratio ($\Lambda < 1.0$ throughout) |
| **Peak-Control Frozen $C_+^*$** | **`332.34`** | **Instantaneous frozen equilibrium at peak control** |
| Dynamic $C(40)$ ($p_T=2$) | `466.43` | Actual 5D trajectory value at peak dynamic control |
| $\Lambda_{\max} / \Lambda_{\text{basal}}$ | `65.51x` | Dynamic control ratio amplification fold |
| Panel A Range | `206.72 – 628.84` | Final tumor density in $(\Delta t, f)$ landscape |
| Panel B Range | `347.20 – 637.93` | Final tumor density in $(\Delta t, B)$ landscape |
| Panel C (Untreated) | `983.46 -> 983.46` | Invariant persistence across $p_T \in [0.2, 4.0]$ |
| Panel C (Fixed) | `819.90 -> 177.11` | Monotherapy response across $p_T$ |
| Panel C (Dynamic) | `875.78 -> 58.46` | Full dynamic harmony response across $p_T$ |

---

## Epistemic Boundary Statements
1. **$\Lambda(t) < 1.0$ Throughout**: No bifurcation crossing occurs in the calibrated regime.
2. **Trajectory vs. Equilibrium Distinction**: $C(40) \approx 983.46$ is the trajectory output, whereas $C_+^* \approx 989.81$ is the 1D frozen-state equilibrium.
3. **No Eradication Claim**: The study demonstrates bifurcation-proximity state contraction, not complete model-mediated eradication.
