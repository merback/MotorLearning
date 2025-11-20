# Learning Electromagnetic Fields Based on Finite Element Basis Functions

This repository implements the methodology from the paper *“Learning Electromagnetic Fields Based on Finite Element Basis Functions”*. It combines **Isogeometric Analysis (IGA)**, **Proper Orthogonal Decomposition (POD)**, and **Deep Neural Networks (DNNs)** to efficiently predict electromagnetic fields in a parametric model of a Permanent Magnet Synchronous Machine (PMSM). Instead of learning field values directly, the network predicts coefficients in a reduced POD basis derived from high-fidelity IGA simulations, ensuring physical consistency and enabling rapid evaluation of quantities of interest such as torque. Two surrogate models are provided:
- **Airgap field surrogate** for torque-focused applications.
- **Full field surrogate** for detailed rotor, stator, and airgap predictions.

## Related Publication
If you use this code, please cite:
```
10.1109/TMAG.2025.3629546
```

## Requirements
- Python 3.11+
- MATLAB R2024b (tested until R2026a) with MATLAB Curve Fitting Toolbox

Install Python dependencies:
```bash
pip install -r requirements.txt
```

Python dependencies:
```
pandas 2.3.3
torch 2.8.0
numpy 2.3.4
matplotlib 3.10.6
```

No other packages need to be downloaded at this point. To provide a fully working snapshot, we already include here the code snapshots for the packages:

- GeoPDEs: [https://rafavzqz.github.io/geopdes/]
- NURBS toolbox: [https://de.mathworks.com/matlabcentral/fileexchange/26390-nurbs-toolbox-by-d-m-spink]
 
We acknowledge these toolboxes and also reused parts of the code published [here](https://github.com/MichaelWiesheu/MotorOptimization/tree/Paper).

---

## Usage
To reproduce the full pipeline:

### 1. Data Generation (MATLAB)
Run:
```bash
matlab -batch "computeSamples(isgap)"
```
Where `isgap` is:
- `true` for airgap field
- `false` for full field

---

### 2. POD Computation (Python)
```bash
python perform_pod.py true
```
Replace `true` with `false` for full field.

---

### 3. Train Neural Network
```bash
python train_net.py true
```
Replace `true` with `false` for full field.

---

### 4. Analyze Predictions
```bash
python analyze_net.py true false
```
Arguments:
- first: `true` for airgap field, `false` for full field
- second: `false` to use your trained network

---

### 5. Compute Errors (MATLAB)
Run:
```bash
matlab -batch "computeErrors(isgap, paper)"
```
Arguments:
- `isgap`: `true` for airgap field, `false` for full field
- `paper`: `false` to use your generated predictions

---

## Output Files
- Directories: `sol_XX_training`,  `sol_XX_testing`, `sol_XX_validation`
- Torch model: `model_XX.pt`
- Predictions: `Predictions_XX_XX`
- Error tables: `Error_XX_XX.csv`

**Important:** The mean, standard deviation, and maximum values of the computed error tables may differ slightly from those in the paper because of numerical errors during data generation, projection matrix computation, and training, as well as minor hyperparameter differences.

---

## Reproducing Paper Results
Download the training, testing, and validation data as well as the two files *MatrixK_full.csv* and *MatrixK_gap.csv* from Zenodo and put them into the main folder:

Run:
```bash
python analyze_net.py true true
```
Then:
```matlab
computeErrors(true, true)
```
