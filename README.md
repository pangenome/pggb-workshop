# Building Pangenome Graphs

Erik Garrison, Julian Lucas, Giulio Formenti, Nadolina Brajuka

For the **_HPRC annual meeting workshop_**, October 11, 2022.

## Learning objectives

This tutorial builds interactive understanding of the [PanGenome Graph Builder (`pggb`)](https://github.com/pangenome/pggb) command line tool. We'll build some graphs and inspect them to understand how the method works and the effects of some of its key parameters.

In this exercise you learn how to

- build pangenome graphs using `pggb`,
- explore `pggb`'s results,
- understand how parameters affect the built pangenome graphs.

There are other methods to build these graphs, like the [minigraph-cactus pipeline](https://doi.org/10.1101/2022.10.06.511217). We're presenting `pggb` because of its easy interactive use and flexibility with diverse inputs of various scales.

### Intro slides

[Slides to get us ready.](https://docs.google.com/presentation/d/1aXwZywy0d3_2sjaTbXTr403ns--mVU1RKDu0NtEt5X4/edit#slide=id.p)

## Getting started

Make sure you have `pggb` and its tools installed.

The easiest way to set things up using `docker`.

    docker pull ghcr.io/pangenome/pggb:latest

Also make sure you have checked out `pggb` repository:

    git clone https://github.com/pangenome/pggb.git

Note that the Docker image is built for `x86_64` and if you're on an M1 Mac or other platform you will need to use `docker build --target binary -t ${USER}/pggb:latest .` in the `pggb` repository to run the build build.

Now create a directory to work on for this tutorial:

    mkdir hprc-workshop
    cd hprc-workshop
    cp -r ~/pggb/data .

Now we set up a docker interactive session, mounting this directory in our `/root` or `$HOME`.

    # run docker with pggb's latest image
    docker run -it -v $(pwd):/root \
        ghcr.io/pangenome/pggb:latest /bin/bash
    cd /root # change into root's $HOME
    ls data  # should show our pggb test data

We can look at the results from outside of the docker container. That tends to be easier as the image doesn't include things like image viewers.

## How does the `pggb` graph build work?

In short, pangenome graphs are multiple alignments. They differ from multiple sequence alignments in that they are nonlinear, and can support any kind of variation that arises in DNA evolution---including inversions and segmental duplications---which form loops and other complex structures in the graph.

### All-to-all alignment

We begin with an alignment, with `wfmash`. This compares all sequences to each other and finds the best `N` mappings for each. It produces base-level alignments.

### Inducing the graph

These base-level alignments are converted into a graph with `seqwish`. A filter is applied to remove short matches, which anchors the graph on confident longer exact matches.

### Normalizing the graph

To normalize the graph and harmonize the allele representation, we use `smoothxg` to apply a local MSA across all parts of the graph.

### Downstream

We can do many things with these graphs. First, we get a number of diagnostic images out of the pipeline, based on the graphs. These give a human interface to the graph models that can help us to understand the alignments at a high level. We're also able to produce variant calls (in `pggb`), using `vg deconstruct`. The graphs from `pggb` can be used as reference systems for short read alignment with `vg giraffe` or long read alignment with `GraphAligner`. Using `odgi` we can use the graphs as reference systems to describe homology relationships between whole genomes.

## Build HLA pangenome graphs

The [human leukocyte antigen (HLA)](https://en.wikipedia.org/wiki/Human_leukocyte_antigen) system is a complex of genes on chromosome 6 in humans which encode cell-surface proteins responsible for the regulation of the immune system.

Let's build a pangenome graph from a collection of sequences of the DRB1-3123 gene:

    pggb -i data/HLA/DRB1-3123.fa.gz -n 12 -t 4 -o DRB1_3123.1

Run `pggb` without parameters to get information on the meaning of each parameter:

    pggb

Take a look at the files in the `DRB1_3123.1` folder.

We get a graph in GFA (`*.gfa`) and odgi (`*.og`) formats. These can be used downstream in many methods, including those in `vg`, like `vg giraffe`. You can visualize the GFA format graph with [`BandageNG`](https://github.com/asl/BandageNG), and use `odgi` directly on the `*.gfa` or `*.og` output.

### Understanding `odgi` visualizations

We obtain a series of diagnostic images that represent the pangenome alignment. These are created with `odgi viz` (1D matrix) and `odgi layout` with `odgi draw` (2D graph drawings).

First, the 2D layout gives us a view of the total alignment. For small graphs, we can look at the version that shows where specific paths go (`*.draw_multiqc.png`):

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

For larger ones, the `*.draw.png` result is usually more legible, but it lacks path information:

![draw.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.lay.draw.png)

We also get some 1D visualizations. These present the graph as a kind of matrix. Across the x-axis we have nodes of the graph (scaled by length) and across the y-axis we have paths, or sequences, which have been embedded in the graph.

This layout is capable of representing several kinds of information using color.

The default associates a color with each path. This is stable across different runs of `odgi viz`:

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.viz_multiqc.png)

We also have a view that shows the "self depth" across the graph.
In this case there are no looping paths, so the color is always gray=1x.

![viz_depth_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.viz_depth_multiqc.png)

We can look at orientation of paths using two views.

One shows the "position" of each path relative to the graph. It runs light to dark from 0 to path length.

![viz_pos_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.viz_pos_multiqc.png)

A similar view shows inverted regions of paths relative to the graph in red, while the forward orientation in black.

![viz_inv_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.viz_inv_multiqc.png)

And finally, a compressed view shows coverage across the pangenome coordinate space of all paths. It's a kind of heatmap. This helps when we have a lot of paths to consider:

![viz_O_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.417fcdf.9c6ea4f.smooth.final.og.viz_O_multiqc.png)

### Looking at the alignments

How many alignments were executed during the pairwise alignment (take a look at the `PAF` output)? Visualize the alignments:

    pafplot -s 2000 DRB1_3123.1/*.paf

Now, from outside the container, use a file browser to open images produced by the process. (On ubuntu linux we can use `eog` to view the PNGs in a whole folder: `eog DRB1_3123.1`.)

![wfmash.paf.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.1/DRB1-3123.fa.gz.510a9ad.wfmash.paf.png)

### Graph statistics and build process

Use `odgi stats` to obtain the graph length, and the number of nodes, edges, and paths.

    odgi stats -i DRB1_3123.1/*.og -S

Do you think the resulting pangenome graph represents the input sequences well? Check the length and the number of the input sequences to answer this question.

### The effect of haplotype count `-n`

What happens if we set a lower `-n`? This parameter determines how many mappings we have.

    pggb -i data/HLA/DRB1-3123.fa.gz -n 4 -t 4 -o DRB1_3123.2

Each sequence is aligned against its `-n`-1 best matches. Setting `-n 4` causes clustering of sequences into groups that are more similar.

Look at the output diagnostic images. Does the graph look more or less compact? (hint: Lowering `-n` causes alignment clustering that breaks the graph into two components.)

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.2/DRB1-3123.fa.gz.cff9e6a.417fcdf.53439a3.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.2/DRB1-3123.fa.gz.cff9e6a.417fcdf.53439a3.smooth.final.og.viz_multiqc.png)

Check graph statistics. Does this pangenome graph represent better or worse the input sequences than the previously produced graph?

    odgi stats -i DRB1_3123.2/*.og -S

Using `pafplot` we can see how this happens.

    pafplot -s 2000 DRB1_3123.2/*.paf

![wfmash.paf.png](https://github.com/pangenome/hprc-workshop/blob/main/DRB1_3123.2/DRB1-3123.fa.gz.cff9e6a.wfmash.paf.png)

### The effect of the minimum match filter `-k`

Another key parameter is `-k`, which affects the behavior of `seqwish`. This filter removes exact matches from alignments that are shorter than `-k`. Short matches occur in regions of high diversity. In practice, these short matches contribute little to the overall structure of the graph, and we remove them to further simplify the base graph structure.

Try setting a much higher `-k` than the default (`-k 19`):

    pggb -i data/HLA/DRB1-3123.fa.gz -n 12 -k 47 -t 4 -o DRB1_3123.3

The graph starts to become "braided". We might say that it is underaligned.

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.3/DRB1-3123.fa.gz.510a9ad.e34d4cd.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.3/DRB1-3123.fa.gz.510a9ad.e34d4cd.9c6ea4f.smooth.final.og.viz_multiqc.png)

We can go lower (try `-k 7` or `-k 0`) or higher (try `-k 79`).

    pggb -i data/HLA/DRB1-3123.fa.gz -n 12 -k 0 -t 4 -o DRB1_3123.4

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.4/DRB1-3123.fa.gz.510a9ad.692a77d.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.4/DRB1-3123.fa.gz.510a9ad.692a77d.9c6ea4f.smooth.final.og.viz_multiqc.png)

### Decreasing mapping segment length `-s` increases sensitivity

At the left end of the graphs, we see a kind of "forked tail" motif in the DRB1-3123 graphs.

We can resolve this by setting a lower mapping segment length, which affects the behavior of the very first step in the pipeline, `wfmash`'s mapping step (itself based on a heavily modified version of MashMap). This defaults to `-s 5k`. We can use `-s 1k` to guarantee we pick up on smaller homology segments, leading to a more complete alignment.

    pggb -i data/HLA/DRB1-3123.fa.gz -s 1k -n 12 -k 0 -t 4 -o DRB1_3123.5

Looking at the alignment plot shows how this works:

![wfmash.paf.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.5/DRB1-3123.fa.gz.db08837.wfmash.paf.png)

The graphs look compact relative to other settings, which can be confirmed by counting base pairs in the graph with `odgi stats -S`.

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.5/DRB1-3123.fa.gz.db08837.692a77d.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.5/DRB1-3123.fa.gz.db08837.692a77d.9c6ea4f.smooth.final.viz_multiqc.png)

### Decreasing the minimum pairwise identity `-p` increases sensitivity

The `-p` setting affects the level of pairwise divergence that's accepted in the mapping step. By dropping this very low, we recover mappings that were missed with the default setting of `-p 90`.

    pggb -i data/HLA/DRB1-3123.fa.gz -p 70 -s 1k -n 12 -k 0 -t 4 -o DRB1_3123.6

We can see the added mappings in the pafplot.

![wfmash.paf.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.6/DRB1-3123.fa.gz.362269e.wfmash.paf.png)

But the effect on graph topology seems minimal, because the missing mappings are collected by the transitive relationships in the other mappings.

![draw_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.6/DRB1-3123.fa.gz.362269e.692a77d.e4629a5.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](https://raw.githubusercontent.com/pangenome/hprc-workshop/main/DRB1_3123.6/DRB1-3123.fa.gz.362269e.692a77d.e4629a5.smooth.final.og.viz_multiqc.png)

### A word of caution...

Note that DRB1-3123 represents a very extreme situation in the human genome---these gene sequences are diverged by up to 20% and lie in the MHC class II region, which is a site of ongoing diversifying selection and frequent incomplete lineage sorting in the primate clade. Parameter settings for whole genomes and chromosomes often are more stringent than those we've tested here (e.g. `-k 79` or even `-k 311` helps to reduce complexity in human satellites).

### Trying other HLA genes

Choose another HLA gene from the `data` folder and explore how the statistics of the resulting graph change as` s`, `p`,` n` change. Produce scatter plots where on the x-axis there are the tested values of one of the `pggb` parameters (`s`, `p`, or `n`) and on the y-axis one of the graph statistics (length, number of nodes, or number of edges). You can do that using the final graph and/or the intermediate ones.

For example:

    pggb -i data/HLA/B-3106.fa.gz -n 9 -t 8 -o B-3106.1

Or

    pggb -i data/HLA/TAP2-6891.fa.gz -n 11 -t 8 -o TAP2-6891.1

To set `-n`, count the lines in the `.fai` index files. This gives the number of sequences in the input:

    wc -l data/HLA/TAP2-6891.fa.gz.fai

## Bonus: LPA pangenome graphs

[Lipoprotein(a) (LPA)](https://en.wikipedia.org/wiki/Lipoprotein(a)) is a low-density lipoprotein variant containing a protein called apolipoprotein(a). Genetic and epidemiological studies have identified lipoprotein(a) as a risk factor for atherosclerosis and related diseases, such as coronary heart disease and stroke.

Try to make LPA pangenome graphs. The input sequences are in `data/LPA/LPA.fa.gz`. Sequences in this locus have a peculiarity: which one? Hint: visualize the alignments and take a look at the graph layout (with `Bandage` and/or in the `.draw_multiqc.png` files).

Here's a hint:

    pggb -i data/LPA/LPA.fa.gz -n 14 -t 8 -o LPA.1
