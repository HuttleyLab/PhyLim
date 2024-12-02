# phylim: a phylogenetic limit evaluation library built on [cogent3](https://cogent3.org/)

phylim evaluates the identifiability when estimating the phylogenetic tree using the Markov model. The identifiability is the key condition of the Markov model used in phylogenetics to fulfil consistency. 

Establishing identifiability relies on the organisation of five types of transition probability matrices on a phylogenetic tree. A key concern arises when a tree does not meet the condition that, for each node, a path to a tip must exist where all matrices along the path are DLC. Such trees are not identifiable 🪚🎄! For instance, in the figure below, tree *T'* contains a node surrounded by a specific type of non-DLC matrix, rendering it non-identifiable. In contrast, compare *T'* with tree *T*.

phylim provides a quick, handy method to check the identifiability of a model fit, where we developed a main [cogent3 app](https://cogent3.org/doc/app/index.html), `phylim`. phylim is compatible with [piqtree2](https://github.com/iqtree/piqtree2), a python library that exposes features from iqtree2.

The following content will demonstrate how to set up phylim and give some tutorials on the main identifiability check app and other associated apps.

<p align="center">
<img src="https://figshare.com/ndownloader/files/50904159" alt="tree1" width="600" height="300" />
</p>

## Installation

```pip install phylim```

Let's see if it has been done successfully. In the package directory:

```pytest```

Hope all tests passed! :blush:

## Run the check of identifiability

If you fit a model to an alignment and get the model result:

```python
>>> from cogent3 import get_app, make_aligned_seqs

>>> algn = make_aligned_seqs(
...    {
...        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
...        "Gorilla": "ATGCGGCGCGCGGAGGCCGCGCTCGCGGAG",
...        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
...    },
...    info={"moltype": "dna", "source": "foo"},
... )

>>> app_fit = get_app("model", "GTR")
>>> result = app_fit(algn)
```

You can easily check the identifiability by:

```python
>>> app_ident_check = get_app("phylim")

>>> record = app_ident_check(result)
>>> record.is_identifiable

True
```

The `phylim` app wraps all information about phylogenetic limits.

```python
>>> record
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
        <td>False</td>
        <td>2024.9.20</td>
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


```
>>> from phylim import check_fit_boundary

>>> violations = check_fit_boundary(result)
>>> violations
BoundsViolation(source='foo', vio=[{'par_name': 'C/T', 'init': np.float64(1.0000000147345554e-06), 'lower': 1e-06, 'upper': 50}, {'par_name': 'A/T', 'init': np.float64(1.0000000625906854e-06), 'lower': 1e-06, 'upper': 50}])
```

</details>


## Check identifiability for piqtree2

phylim provides an app, `phylim_tree_to_likelihoodfunction`, which allows you to build the likelihood function from a piqtree2 output tree.

```python
>>> from piqtree2 import build_tree

>>> tree = build_tree(algn, model="GTR")
>>> app_inverter = get_app("phylim_tree_to_likelihoodfunction")

>>> result = app_inverter(tree)
>>> record = app_ident_check(result)
>>> record.is_identifiable

True
```


## Colour the edges for a phylogenetic tree based on matrix categories

If you obtain a model fit, phylim can visualise the tree with labelled matrices. 

phylim provides an app, `phylim_colour_edges`, which takes an edge-matrix category map and colours the edges:

```python
>>> from phylim import classify_model_psubs

>>> edge_to_cat = classify_model_psubs(result)
>>> tree = result.tree

>>> app_colour_edge = get_app("phylim_colour_edges", edge_to_cat)
>>> app_colour_edge(tree)
```

<img src="https://figshare.com/ndownloader/files/50903022" alt="tree1" width="400" />


You can also colour edges using a user-defined edge-matrix category map, applicable to any tree object! 

```python
>>> from cogent3 import make_tree
>>> from phylim.classify_matrix import SYMPATHETIC, DLC

>>> tree = make_tree("(A, B, C);")
>>> edge_to_cat = {"A":SYMPATHETIC, "B":SYMPATHETIC, "C":DLC}

>>> app_colour_edge = get_app("phylim_colour_edges", edge_to_cat)
>>> app_colour_edge(tree)
```

<img src="https://figshare.com/ndownloader/files/50903019" alt="tree1" width="400" />

