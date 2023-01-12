# Week 1

In this weeks tutorial you are tasked with installing `FastQC` & `MultiQC` quality control tools to assess the quality of sequencing data provided under `week1/data`.

Following quality control assessment, install `TrimGalore!` to perform adapter trimming and low quality read filtering. You will re-run `FastQC` and `MultiQC` to check if the adapter trimming was successful.

## Git housekeeping

Make a clone of the GitHub repository `BarryDigby/MA5112` locally on your laptop:

```console
git clone git@github.com:BarryDigby/MA5112.git
```

Checkout a branch with the following naming convention: first character of your first name + surname. (e.g Barry Digby = bdigby )

```console
cd MA5112/
git checkout -b <your name>
```

You will now have a copy of the `main` branch to work on for the tutorial. Each week, I will release new materials on the `main` branch. You can work on your own branch and incorporate these new additions using the following commands:

```console
git pull origin/main
```

This command essentially means "download the latest materials from the `main` branch to my branch". You will get the following error the first time you run this:

```console
From github.com:BarryDigby/MA5112
 * branch            main       -> FETCH_HEAD
hint: You have divergent branches and need to specify how to reconcile them.
hint: You can do so by running one of the following commands sometime before
hint: your next pull:
hint: 
hint:   git config pull.rebase false  # merge
hint:   git config pull.rebase true   # rebase
hint:   git config pull.ff only       # fast-forward only
hint: 
hint: You can replace "git config" with "git config --global" to set a default
hint: preference for all repositories. You can also pass --rebase, --no-rebase,
hint: or --ff-only on the command line to override the configured default per
hint: invocation.
fatal: Need to specify how to reconcile divergent branches.
```

What Git is saying is: The branch `main` has some new changes and you want to merge them, but your branch has changes too. How should I handle this?

We want to use `git config pull.rebase false` to accept the latest changes from `main` without overwriting your own work.

To recap: for week 2 when I release new materials, run the command `cd MA5112` and then `git pull origin/main` to synchronise the latest materials to your fork.

## Installing tools

Using `Anaconda` create an environment called `week1`:

```console
conda create -n week1
conda activate week1
```

Install the tools required for the practical:

```console
conda install -c bioconda fastqc multiqc trim-galore
```

## 

`fastqc --help`

`multiqc --help`

`trim_galore --help`

`for file in data/* ; do trim_galore --cores 2 --adapter TGGAATTCTCGGGTGCCAAGG --length 17 --clip_r1 4 --three_prime_clip_r1 4 --max_length 40 --gzip $file; done`

`for file in data/*; do fastqc $file --outdir qc/; done`





