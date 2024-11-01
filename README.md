# Qombucha

**Qombucha** is a novel combinatorial optimization-based framework that leverages known brain developmental history to simultaneously infer the methylation profiles of progenitor brain cells and cell type compositions from the methylation data of bulk glioblastoma samples. Given a patient-by-methylation features matrix matrix `B` of bulk methylation data, and an incomplete progenitor cell reference panel (containing `numpy.nan`s) `C_p` of cell types-by-methylation feature matrix, and the best estimated guess of the full reference panel matrix `C_0`, `Qombucha` uses a block coordinate descent-like approach to iteratively update cell type composition matrix `P` and unknown entries in the progenitor profiles. It outputs the final matrices `P` and `C` after convergence. 

## Installation from source

### Gurobi

Qombucha uses the Gurobi solver to iteratively solve quadratic programs. You need a valid license to use Gurobi. Licenses are [free for academic purposes](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).

### Conda

The easiest way to run Qombucha is to install all dependencies through the [conda](https://docs.conda.io/en/latest/) environment file provided in the repository. Once you have `conda` installed in your system, run the following commands, to create a conda environment named `qombucha`:

```
git clone https://github.com/algo-cancer/Qombucha.git
cd Qombucha
conda env create -f qombucha-env.yaml
conda activate qombucha
```

### Installing the Qombucha Python package
`TODO`.


## Running Qombucha:

### Inputs
The input matrices that are required are: `B`, `C`, and `C_0` (as described in our manuscript)
They should be `numpy arrays`, saved using the command: `numpy.save(<file>, <array>, allow_pickle=True)`

Once the installed `conda` environment `qombucha` is activated, navigate to `<project_base_dir>/src/Qombucha` and run `python Qombucha.py <OPTIONS>`. For a list of valid `OPTIONS`, run the following command:

```
$ python Qombucha.py -h
usage: Qombucha.py [-h] -C REFERENCE -Cp REFERENCE_PARTIAL -B BULK -op OUTPUT_P -oc OUTPUT_C [-I ITERS] [-c THREADS] [-t MAX_RUN_TIME_PER_GUROBI_INSTANCE] [-Delta DELTA_PERCENT_OBJ_ERROR_THRESHOLD] [-Lambda LAMBDA_DIFF_SUCCESSIVE_DELTA] [-W WEIGHT_PER_PATIENT]
                   [-Alpha ALPHA_WEIGHT_REGULARIZATION]

options:
  -h, --help            show this help message and exit
  -C REFERENCE, --reference REFERENCE
                        Complete reference panel
  -Cp REFERENCE_PARTIAL, --reference_partial REFERENCE_PARTIAL
                        Partial reference panel
  -B BULK, --bulk BULK  Bulk methylation value matrix
  -op OUTPUT_P, --output_P OUTPUT_P
                        Path to output P matrix (cell type fractions per patient)
  -oc OUTPUT_C, --output_C OUTPUT_C
                        Path to output C matrix (cell types by methylation features)
  -I ITERS, --iters ITERS
                        Max number of iterations (default = 200)
  -c THREADS, --threads THREADS
                        Number of threads to use (default = 2)
  -t MAX_RUN_TIME_PER_GUROBI_INSTANCE, --max_run_time_per_gurobi_instance MAX_RUN_TIME_PER_GUROBI_INSTANCE
                        Max runtime per Gurobi instance in seconds (default = 100000)
  -Delta DELTA_PERCENT_OBJ_ERROR_THRESHOLD, --delta_percent_obj_error_threshold DELTA_PERCENT_OBJ_ERROR_THRESHOLD
                        Convergence threshold value for relative objective error (default = 0.0125)
  -Lambda LAMBDA_DIFF_SUCCESSIVE_DELTA, --lambda_diff_successive_delta LAMBDA_DIFF_SUCCESSIVE_DELTA
                        Convergence threshold value for change in successive deltas (default = 1e-9)
  -W WEIGHT_PER_PATIENT, --weight_per_patient WEIGHT_PER_PATIENT
                        Experimental: Weight per patient to be used in the objective (Not used by default)
  -Alpha ALPHA_WEIGHT_REGULARIZATION, --alpha_weight_regularization ALPHA_WEIGHT_REGULARIZATION
                        Experimental: Alpha for regularizing imputed C values. (default = 0, not used)
```
