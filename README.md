## Key Files & Variables

### `ML_project.m`
Main **MATLAB** script that consolidates outputs from all paradigms and generates the final dataset.

> *Note:*  
> Preprocessing was conducted in MATLAB using paradigm-specific scripts developed with internal lab functions.  
> These scripts are **not publicly available** due to institutional restrictions.  
> However, the main integration script (`ML_project.m`), which merges and structures the final dataset, *is* included in this repository.

---

### `data_matrix.mat`
⚠️ Due to GitHub file size limitations, this file is hosted externally:

https://drive.google.com/file/d/13OkyKUDMG6gW2IvOWi_5T4VTVr3qVILM/view?usp=drive_link

EEG features extracted from **selected electrodes only**:  
- **7 Alpha electrodes:** `P1`, `P2`, `P4`, `PO3`, `PO4`, `POz`, `Pz`  
- **8 Beta electrodes:** `C2`, `Cz`, `F1`, `F2`, `FC1`, `FC2`, `FCz`, `Fz`  

Each trial includes:
- **Maximum Power**  
- **Mean Power**  
- **Frequency at Maximum Power**

Additionally, each trial includes a **FOOOF vector** capturing the spectral decomposition of power.  
These vectors were averaged across selected electrodes to produce a compact frequency-domain representation.

> *FOOOF vector shapes vary across paradigms due to trial length:*  
> - RS Open: `(1, 1519)`  
> - RS Closed: `(1, 1519)`  
> - Oddball: `(1, 3295)`  
> - Lecture: `(1, 3346)`
---

### `data_matrix_64.mat`
⚠️ Due to GitHub file size limitations, this file is hosted externally:

https://drive.google.com/file/d/1J1W4xjJmFJP3d5mE88BZFfzmUfDw-9Yy/view?usp=sharing

Same structure as `data_matrix`, but includes **all 64 electrodes**, resulting in:  
- Enhanced spatial resolution  
- Expanded feature space  

> *Note:* The **FOOOF vector** remains unchanged and is identical across both datasets.
---

### `ML final project.ipynb`
Main **Python** script for training models and comparing configurations.  
Fully compatible with **Google Colab** for easy execution and reproducibility.

---
