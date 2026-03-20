# phylim: a phylogenetic limit evaluation library built on [cogent3](https://cogent3.org/)
[![Coverage Status](https://coveralls.io/repos/github/HuttleyLab/PhyLim/badge.svg?branch=main)](https://coveralls.io/github/HuttleyLab/PhyLim?branch=main)
[![Release](https://github.com/HuttleyLab/PhyLim/actions/workflows/release.yml/badge.svg)](https://github.com/HuttleyLab/PhyLim/actions/workflows/release.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15717906.svg)](https://doi.org/10.5281/zenodo.15717906)


phylim evaluates the identifiability when estimating the phylogenetic tree using the Markov model. The identifiability is the key condition of the Markov model used in phylogenetics to fulfil consistency. 

Establishing identifiability relies on the arrangement of specific types of transition probability matrices (e.g., DLC and sympathetic) while avoiding other types. A key concern arises when a tree does not meet the condition that, for each node, a path to a tip must exist where all matrices along the path are DLC. Such trees are not identifiable 🪚🎄! For instance, in the figure below, tree *T'* contains a node surrounded by a specific type of non-DLC matrix, rendering it non-identifiable. In contrast, compare *T'* with tree *T*.

phylim provides a quick, handy method to check the identifiability of a model fit, where we developed a main [cogent3 app](https://cogent3.org/doc/app/index.html), `phylim`. phylim is compatible with [piqtree](https://github.com/iqtree/piqtree), a python library that exposes features from iqtree2.

The following content will demonstrate how to set up phylim and give some tutorials on the main identifiability check app and other associated apps.

<p align="center">
<img src="https://ndownloader.figshare.com/files/50904159" alt="tree1" width="600" height="300" />
</p>


## Installation

```pip install phylim```

Let's see if it has been done successfully. In the package directory:

```pytest```

Hope all tests passed! :white_check_mark: :blush:

## Run the check of identifiability

If you fit a model to an alignment and get the model result:

```python
>>> from cogent3 import get_app, make_aligned_seqs

>>> aln = make_aligned_seqs(
...    {
...        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
...        "Gorilla": "ATGCGGCGCGCGGAGGCCGCGCTCGCGGAG",
...        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
...    },
...    info={"moltype": "dna", "source": "foo"},
... )

>>> app_fit = get_app("model", "GTR")
>>> result = app_fit(aln)
```

You can easily check the identifiability by:

```python
>>> checker = get_app("phylim")

>>> checked = checker(result)
>>> checked.is_identifiable

True
```

The `phylim` app wraps all information about phylogenetic limits.

```python
>>> checked
```


<div class="c3table">
  <table>
    <thead class="head_cell">
      <tr>
        <th>Source</th>
        <th>Model Name</th>
        <th>Identifiable</th>
        <th>Has Boundary Values</th>
        <th>Version</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>brca1.fasta</td>
        <td>GTR</td>
        <td>True</td>
        <td>True</td>
        <td>2025.1.12</td>
      </tr>
    </tbody>
  </table>
</div>


You can also use features like classifying all matrices or checking boundary values in a model fit.

<details>
<summary>Label all transition probability matrices in a model fit</summary>


You can call `classify_model_psubs` to give the category of all the matrices:

```python
>>> from phylim import classify_model_psubs

>>> labelled = classify_model_psubs(result)
>>> labelled
```


<div class="c3table">
<table>

<caption>
<span class="cell_title">Substitution Matrices Categories</span>
</caption>
<thead class="head_cell">
<th>edge name</th><th>matrix category</th>
</thead>
<tbody>
<tr><td><span class="c3col_left">Gorilla</span></td><td><span class="c3col_left">DLC</span></td></tr>
<tr><td><span class="c3col_left">Human</span></td><td><span class="c3col_left">DLC</span></td></tr>
<tr><td><span class="c3col_left">Mouse</span></td><td><span class="c3col_left">DLC</span></td></tr>
</tbody>
</table>

</div>

</details>


<details>
<summary>Check if all parameter fits are within the boundary</summary>


```python
>>> from phylim import check_fit_boundary

>>> violations = check_fit_boundary(result)
>>> violations
BoundsViolation(source='foo', vio=[{'par_name': 'C/T', 'init': np.float64(1.0000000147345554e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000625906854e-06), 'lower': 1e-06, 'upper': 50}])
```

</details>

❗For users who want to check identifiability on a model with multiple likelihood functions (e.g. **split codon model**), please check https://github.com/HuttleyLab/PhyLim/issues/23#issuecomment-3125670158


## Check identifiability for piqtree

phylim provides an app, `phylim_to_model_result`, which allows you to build the likelihood function from a piqtree output tree.

```python
>>> phylo = get_app("piq_build_tree", model="GTR")
>>> tree = phylo(aln)

>>> lf_from = get_app("phylim_to_model_result")
>>> result = lf_from(tree)

>>> checker = get_app("phylim")
>>> checked = checker(result)
>>> checked.is_identifiable

True
```




## Colour the edges for a phylogenetic tree based on matrix categories

If you obtain a model fit, phylim can visualise the tree with labelled matrices. 

phylim provides an app, `phylim_style_tree`, which takes an edge-matrix category map and colours the edges:

```python
>>> from phylim import classify_model_psubs

>>> edge_to_cat = classify_model_psubs(result)
>>> tree = result.tree

>>> tree_styler = get_app("phylim_style_tree", edge_to_cat)
>>> tree_styler(tree)
```

<img src="https://ndownloader.figshare.com/files/50903022" alt="tree1" width="400" />


You can also colour edges using a user-defined edge-matrix category map, applicable to any tree object! 

```python
>>> from cogent3 import make_tree
>>> from phylim import SYMPATHETIC, DLC

>>> tree = make_tree("(A, B, C);")
>>> edge_to_cat = {"A":SYMPATHETIC, "B":SYMPATHETIC, "C":DLC}

>>> tree_styler = get_app("phylim_style_tree", edge_to_cat)
>>> tree_styler(tree)
```

<img src="https://ndownloader.figshare.com/files/50903019" alt="tree1" width="400" />
