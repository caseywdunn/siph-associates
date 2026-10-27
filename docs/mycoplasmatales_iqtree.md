# Inferring the bacterial phylogeny with IQ-TREE

The reported bacterial phylogeny uses IQ-TREE 3.0.1 on the existing 483-genome,
5,035-amino-acid bac120 alignment, with Acholeplasmataceae as the outgroup.
This inference reuses the prepared markers without regenerating assemblies,
read mappings, genome clusters or geographical tests.

The alignment is in `data/results/mycoplasmatales_phylogeny/`; the original
FastTree output remains there as a historical result. The accepted IQ-TREE
outputs are in `data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/`.
Job 13941115 completed successfully on 2026-10-26. The selected
`bac120.treefile` uses **LG+F+I+R10**; `bac120.iqtree` records model selection,
search likelihoods and support calculations, and `inference_complete.json`
records output integrity checks. Selected public outputs are in
[`tables/mycoplasmatales/iqtree/`](../tables/mycoplasmatales/iqtree/).

## Run

For a required reproduction, submit from the analysis repository root. The
accepted run is complete; viewing or rendering it does not require another run.

```bash
sbatch scripts/revision_iqtree.sbatch
```

The launcher uses the installed IQ-TREE 3.0.1 and Biopython environment; its
dependencies are specified in `envs/iqtree.yaml`. It requests 16 CPUs, 64 GB RAM,
and 24 hours on `ycga`, with all working files on GPFS. ModelFinder selects a
protein model by BIC. Each of five independent maximum-likelihood searches
includes 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates.
Bootstrap trees receive additional NNI optimization (`--bnni`). The best-scoring
run supplies the reported tree and support values; replicates are not pooled
across runs. The random seed is 20261025. The full executed command and
software identity are recorded.

Input checks reconcile alignment and tree identifiers, measure missing data,
and record identical sequences and repeated source assemblies. All input taxa
are retained; IQ-TREE may collapse identical sequences internally and restore
them in its output. Output validation requires all 483 tips. These duplicates
do not constitute independent biological evidence.

Submitting the same unchanged launcher resumes IQ-TREE checkpoints after an
interruption. Hash checks reject changed inputs, implementation, or executable.
Execution records are retained separately for every SLURM job.

## Interpretation and plotting

Node labels report **SH-aLRT / ultrafast bootstrap support**, both as percentages.
The group containing the seven focal representatives has support 99.9/100.
Clades A and B, MAGSP0031's sister relationship to clade B, DT-68 and the
octocoral genomes, and Mycoplasma_K retain their reported relationships.
The clade-B/MAGSP0031-to-DT-68 relationship remains weak (64.2/71).

The five searches agree on the focal topology. Only one reached the best
full-tree likelihood, and 348/483 taxa failed the composition test. These
diagnostics limit interpretation of deeper relationships; high support from
concatenated markers does not establish model adequacy or agreement among
markers. Weakly supported branches remain unresolved biological relationships.

Figure 3 displays the 21-tip focal group with branch lengths and support from
the accepted full tree. It represents 20 assemblies because one source
assembly also occurs among the GTDB references. The plotting step uses the
verified 94-tip outgroup split when orienting the tree; IQ-TREE's displayed
first outgroup tip alone is not the biological root. See the
[comparative-context results](mycoplasmatales_context_results.md) and
[phylogeny workflow guide](workflow/mycoplasmatales_phylogeny/README.md).
