# Checking the bacterial phylogeny with IQ-TREE

This analysis tests whether the bacterial relationships reported from FastTree
persist under IQ-TREE. It reuses the existing 483-genome, 5,035-amino-acid bac120
alignment and its Acholeplasmataceae outgroup. It does not regenerate assemblies,
markers, read mappings, genome clusters, or geographical tests.

The original alignment and FastTree remain in
`data/results/mycoplasmatales_phylogeny/`. New results are written separately to
`data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/`.

## Run

From the analysis repository root, submit:

```bash
sbatch scripts/revision_iqtree.sbatch
```

The launcher uses the installed IQ-TREE 3.0.1 and Biopython environment; its
dependencies are specified in `envs/iqtree.yaml`. It requests 16 CPUs, 64 GB RAM,
and 24 hours on `ycga`, with all working files on GPFS. ModelFinder selects a
protein model by BIC. Five independent maximum-likelihood searches accompany
1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates. Bootstrap
trees receive additional NNI optimization (`--bnni`). The random seed is
20261025. The full executed command and software identity are recorded.

Input checks reconcile alignment and tree identifiers, measure missing data,
and record identical sequences and repeated source assemblies. All input taxa
are retained; IQ-TREE may collapse identical sequences internally and restore
them in its output. Output validation requires all 483 tips. These duplicates
do not constitute independent biological evidence.

Submitting the same unchanged launcher resumes IQ-TREE checkpoints after an
interruption. Hash checks reject changed inputs, implementation, or executable.
Execution records are retained separately for every SLURM job.

## Interpretation

Successful inference does not automatically replace the manuscript tree. Inspect
the selected model, sequence-composition diagnostics, search likelihoods,
bootstrap convergence, topology, and both support measures first. If only one
search reaches the best likelihood, additional searches may be needed.

Compare the group containing the seven focal genomes, clades A and B,
MAGSP0031's placement, DT-68 and the octocoral genomes, Mycoplasma_K, and the
weakly supported clade-B/DT-68 relationship. FastTree local supports and IQ-TREE
bootstrap or SH-aLRT values are different measures. High concatenated-marker
support alone does not establish model adequacy or agreement among markers.

Review whether composition diagnostics or changes in the important branches
warrant a more complex protein model. Update the figure, methods, citations,
and biological interpretation only after this comparison; report disagreements
and weakly supported relationships explicitly.
