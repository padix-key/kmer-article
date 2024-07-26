# Raw files for the article 'A fast and simple approach to k-mer decomposition'

## Running benchmarks

The k-mer decomposition method is implemented in *Rust* and can be executed with

```
cargo run --release
```

which generates the run times for each method.
To generate the plot from the raw run times, run

```
python plot.py
```

This script not only create the benchmark figure, but the linear fits as well.

## Compiling article

The article is written using [`Typst`](https://github.com/typst/typst).
To compile the document after the benchmark was created, run

```
typst compile article.typ article.pdf
```
