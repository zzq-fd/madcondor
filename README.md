# MadGraph Parameter Scan Setup

This setup is designed for running parameter scans for MadGraph-produced processes. It automates the generation of parameter combinations, submission to HTCondor, and extraction of cross-section results.

## Overview

The workflow scans over multiple physics parameters:
- **sinθ** (sintheta): Mixing angle values
- **tanβ** (tanbeta): Ratio of vacuum expectation values
- **mA**: Input mass parameter for particle 36 (i.e. `mass(36)`)
- **Δm** (deltam): Defined as `m35 - m36` (so `m35 = mA + deltam`)
- **m35/m37**: Derived (we set `m37 = m35`)
- **m55** (mass55): Mass parameter for particle 55
- **mχ** (mchi): Mass parameter for particle 52
- **cos(β-α)** (cosbma): User input; the param card uses `sin(β-α) = sqrt(1 - cosbma^2)` for `FRBLOCK 3`

## Files Description

### Core Scripts

- **`make_joblist.py`**: Generates a job list file (`joblist.txt`) with all parameter combinations to scan
- **`runMadscan.sh`**: Executes a single MadGraph job for a given parameter set
- **`subMadscan.sub`**: HTCondor submission file for batch job processing
- **`extract_cross_sections.py`**: Extracts cross-sections from completed MadGraph runs and creates a summary table
- **`launch_card_template.dat`**: Template for MadGraph launch cards with parameter placeholders

### Data Files

- **`bbdm_2HDMa_type1_case1_scan.tar.gz`**: Compressed MadGraph process directory containing the physics model and process definition
- **`joblist.txt`**: Generated file containing all parameter combinations (one per line: `sintheta tanbeta mA(mass36) m55 mchi cosbma deltam`)
- **`generate_card.txt`**: MadGraph process generation card for creating the initial process directory

## Initial Setup: Creating the MadGraph Process Directory

Before running parameter scans, you need to set up MadGraph and generate the process directory. Follow these steps:

### Step 1: Clone This Directory

Clone or download this repository to your local machine:

```bash
git clone <repository-url>
cd madcondor
```

Or if you already have the directory, navigate to it:

```bash
cd /path/to/madcondor
```

### Step 2: Download and Install MadGraph5_aMC@NLO

Download MadGraph5_aMC@NLO from the official website:

```bash
# Download MadGraph (adjust version as needed)
wget https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.0.tar.gz

# Extract the archive
tar -xzf MG5_aMC_v3.5.0.tar.gz

# Navigate to the MadGraph directory
cd MG5_aMC_v3_5_0
```

**Note:** Replace the version number with the appropriate version you need. You can also download from: https://launchpad.net/mg5amcnlo/

### Step 3: Install the Model

The model directory `Pseudoscalar_2HDMI` is included in this repository. Copy it to MadGraph's models directory (can also be found here https://github.com/LHC-DMWG/model-repository/):

```bash
# From the madcondor directory, copy the model to MadGraph's models directory
# Make sure you're in the madcondor directory first
cd /path/to/madcondor

# Copy the model directory to MadGraph's models folder
cp -r Pseudoscalar_2HDMI /path/to/MG5_aMC_v3_5_0/models/

# Verify the model was copied
ls /path/to/MG5_aMC_v3_5_0/models/Pseudoscalar_2HDMI
```

**Note:** The model directory is named `Pseudoscalar_2HDMI`, but in `generate_card.txt` it's referenced as `Pseudoscalar_2HDMI-bbMET_5FS`. This is because MadGraph automatically appends the restriction file name (`restrict_bbMET_5FS.dat`) to the model name. The model directory itself should be copied as `Pseudoscalar_2HDMI`.

### Step 4: Run MadGraph with generate_card.txt

Run MadGraph with the generation card using input redirection:

```bash
# From the MG5_aMC_v3_5_0 directory, run with input redirection
./bin/mg5_aMC /path/to/madcondor/generate_card.txt
```

This will create a process directory named `bbdm_2HDMa_type1_case1_scan` (as specified in the `output` command in `generate_card.txt`).

### Step 5: Modify SubProcesses/setcuts.f

After the process is generated, you need to modify the cuts file to add special handling for the xd xd~ particles (particle ID 52) in the missing Et block.

Navigate to the generated process directory:

```bash
cd bbdm_2HDMa_type1_case1_scan
```

Edit the file `SubProcesses/setcuts.f` and locate the section that handles c-neutrino's (missing Et block). Add the following line in that section:

```fortran
if (abs(idup(i,1,iproc)).eq.52) is_a_nu(i)=.true.  ! no cuts on xd xd~
```

**How to find the right location:**

1. Open `SubProcesses/setcuts.f` in a text editor
2. Search for the section that sets `is_a_nu(i)=.true.` for neutrinos
3. Add the line above in that section, typically near other neutrino identification lines

**Example location (the exact line numbers may vary):**

```fortran
      do i=1,nexternal
        if (abs(idup(i,1,iproc)).eq.12) is_a_nu(i)=.true.  ! electron neutrino
        if (abs(idup(i,1,iproc)).eq.14) is_a_nu(i)=.true.  ! muon neutrino
        if (abs(idup(i,1,iproc)).eq.16) is_a_nu(i)=.true.  ! tau neutrino
        if (abs(idup(i,1,iproc)).eq.52) is_a_nu(i)=.true.  ! no cuts on xd xd~
      enddo
```

### Step 6: Modify Cards/run_card.dat

Edit the run card to add a missing Et cut and turn off systematics:

```bash
# Edit the run card
nano Cards/run_card.dat
# or use your preferred editor
```

**A. Add Missing Et Cut:**

Find the section `# Minimum and maximum pt's (for max, -1 means no cut)` and add the following line in that block:

```
150.0  = misset    ! minimum missing Et (sum of neutrino's momenta)
```

**B. Turn Off Systematics:**

Find the section `# Store info for systematics studies` and set systematics to off (typically by setting the flag to `False` or `0`, depending on the card format).

**Example modifications:**

```bash
# In the "Minimum and maximum pt's" block, add:
150.0  = misset    ! minimum missing Et (sum of neutrino's momenta)

# In the "Store info for systematics studies" block, set:
False  = use_syst  ! Enable systematics studies
```

### Step 7: Create the Process Archive

After making all modifications, create a tar.gz archive of the process directory for use with the parameter scan:

```bash
# From the parent directory of bbdm_2HDMa_type1_case1_scan
cd ..
tar -czf bbdm_2HDMa_type1_case1_scan.tar.gz bbdm_2HDMa_type1_case1_scan/

# Copy it to your madcondor directory
cp bbdm_2HDMa_type1_case1_scan.tar.gz /path/to/madcondor/
```

**Note:** This archive will be used by the HTCondor submission script. Make sure the path in `subMadscan.sub` points to this file.

### Verification Checklist

Before proceeding to parameter scans, verify:

- [ ] Model directory is installed in `MG5_aMC_v3_5_0/models/`
- [ ] Process directory `bbdm_2HDMa_type1_case1_scan` was generated successfully
- [ ] `SubProcesses/setcuts.f` contains the line for particle ID 52
- [ ] `Cards/run_card.dat` has the missing Et cut (150.0)
- [ ] `Cards/run_card.dat` has systematics turned off
- [ ] Process archive `bbdm_2HDMa_type1_case1_scan.tar.gz` was created

## Setup Instructions

### 1. Prerequisites

- **MadGraph5_aMC@NLO**: The MadGraph installation must be available
- **HTCondor**: For batch job submission (if using cluster)
- **Python 3**: For running the job list generator and cross-section extraction script
- **Bash**: For running the scan script

### 2. Configure Parameters in `make_joblist.py`

`make_joblist.py` now has a user-editable block at the top:

```python
USER_SCAN_INPUT = {
    "SINTHETA": 0.9,
    "TANBETA": 3,
    "M36": 1500,
    "M55": 300,
    "MCHI": "scan:[140,140.5,141]",
    "COSBMA": 0,
    "DELTAM": 0,
}
```

Supported input syntax:

- Fixed value (no scan): `3` or `0.9`
- Independent scan: `scan:[1,2,3]`
- Synchronized scan: `scan1:[...]`, `scan2:[...]`, ...

Example synchronized scan:

```python
"M36": "scan1:[600,700,800]",
"DELTAM": "scan1:[0,50,100]",
```

This generates aligned pairs by index:
- `(M36, DELTAM) = (600,0)`
- `(M36, DELTAM) = (700,50)`
- `(M36, DELTAM) = (800,100)`

Core modes kept in `make_joblist.py`:

- Default mode: use `USER_SCAN_INPUT`
- `--from-results`: rebuild joblist from existing `Events_*`
- `--preset`: built-in grids (`type1case1scan1`, `tan0p01to50`, `legacy_default`)
- `--repeat N`: repeat each generated point `N` times

### 3. Generate the Job List

Default mode (use top config):

```bash
python3 make_joblist.py -o joblist.txt
```

Alternative examples:

```bash
python3 make_joblist.py --preset type1case1scan1 -o joblist.txt
python3 make_joblist.py --from-results bbdm_2HDMa_type1_case1_scan1 -o joblist.txt
python3 make_joblist.py --repeat 5 -o joblist.txt
```

### 4. Configure Submission Files

Minimal checks before submit:

- `launch_card_template.dat` has placeholders used by `runMadscan.sh`:
  - `__TANBETA__`, `__SINTHETA__`, `__SINBMA__`
  - `__M35__`, `__M36__`, `__M37__`, `__M55__`, `__MCHI__`
- `subMadscan.sub` points to the correct queue file (`joblist.txt`)
- Resource settings in `subMadscan.sub` are appropriate:
  - `request_memory`, `request_cpus`, `+JobFlavour`

### 5. Submit to HTCondor

Recommended submit method:

```bash
./submitMadscan.sh bbdm_2HDMa_type1_case1_scan.tar.gz
```

This helper script automatically:

- Derives `OUTROOT` from tarball name
- Creates log/result folders under `OUTROOT`
- Passes `PROCESS_TAR` and `OUTROOT` to `condor_submit`

### 6. Monitor Jobs

```bash
condor_q
```

Logs are written to:

- `OUTROOT/logs/output/`
- `OUTROOT/logs/error/`
- `OUTROOT/logs/log/`

### 7. Extract Cross Sections

After jobs finish:

```bash
python3 extract_cross_sections.py OUTROOT
```

Output:

- `OUTROOT/cross_section_table.txt`

### 8. Plotting Workflow

All plotting scripts below read existing run outputs under `OUTROOT/results/`.

Cross-section scatter:

```bash
python3 plot_crosssection_scatter.py OUTROOT --xmin 1 --xmax 10
```

MET overlays:

```bash
python3 plot_met55_overlay.py OUTROOT --sin 0.9 --tan 3 --mA 600 --ma 300 --mchi 1 --cosbma 0 --deltam 0
python3 plot_met52.py OUTROOT --sin 0.9 --tan 3 --mA 600 --ma 300 --mchi 1 --cosbma 0 --deltam 0
```

Single-point pT distributions:

```bash
python3 plot_pt55.py OUTROOT --xmin 0 --xmax 500
python3 plot_pt52.py OUTROOT --xmin 0 --xmax 500
```

Type1/Type2 comparisons:

```bash
python3 comparison_crosssection.py bbdm_2HDMa_type1_case1_scan bbdm_2HDMa_type2_case1_scan --xmax 10
python3 comparison_met.py case1_scan --sin 0.9 --tan 3 --mA 600 --ma 300 --mchi 1 --cosbma 0 --deltam 0
```

### 9. Minimal End-to-End Example

```bash
python3 make_joblist.py -o joblist.txt
./submitMadscan.sh bbdm_2HDMa_type1_case1_scan.tar.gz
python3 extract_cross_sections.py bbdm_2HDMa_type1_case1_scan
python3 plot_crosssection_scatter.py bbdm_2HDMa_type1_case1_scan
```
