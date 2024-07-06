#import "@preview/charged-ieee:0.1.0": ieee


#let kmer = [$k$-mer]
#let kmers = [$k$-mers]

#let todo(msg) = {
  [#text(fill: red, weight: "bold", size: 10pt)[TODO #msg]]
}

#let example(content) = {
  box(inset: 10pt)[#text(fill: luma(100))[*Example:* \ #content]]
}

#let seq(content) = [$mono(#content)$]


#show: ieee.with(
  title: [A fast and simple approach to #kmer decomposition],
  abstract: [
    Alignment searches are fast heuristic methods to identify similar regions between two sequences.
    This group of algorithms is ubiquitously used in a myriad of software to find homologous
    sequences or to map sequence reads to genomes.
    Often the first step in alignment searches is #kmer decomposition: listing all overlapping
    subsequences of length $k$.
    This article presents a simple integer representation of #kmers and shows how a sequence can
    be quickly decomposed into #kmers in constant time with respect to $k$.
  ],
  authors: (
    (
      name: "Patrick Kunzmann",
      // department: [],
      // organization: [],
      location: [Göttingen, Germany],
      email: "patrick.kunzm@gmail.com"
    ),
  ),
  index-terms: ("k-mers", "Alignment"),
  bibliography: bibliography("refs.bib"),
)


= Introduction
A common task in bioinformatics is finding similar regions between sequences.
In one scenario finding homologs of a sequence in another genome may reveal common functionalities
between them #todo(cite).
In another scenario the reads obtained from sequencing need to be mapped to a position on a
reference genome to assemble a genome or quantify the number of transcripts #todo(cite).
A _similar region_ can be formalized as so called _alignment_:
It specifies which position in one sequence corresponds to a position the the other sequence.
The dynamic programming algorithms to obtain the guaranteed best alignment solution @Needleman1970
is not computationally feasible for most modern applications: the length and number of sequences is
simply too large.

To solve this problem heursitic approaches emerged. Many modern algorithms (see #todo(cite) as
examples) build upon the concept of finding exact matches of length $k$ between the sequences
@Altschul1990.
These subsequences of length $k$ are termed #kmers.
The process of finding all overlapping #kmers in a sequence is called _k-mer decomposition_.

#example[
  Take the sequences #seq[TATGC] and #seq[ATGG]:
  Decomposition into 3-mers yields the enumeration
  #footnote[
    The more common mathematical term would be '_sequence_'.
    For better distinction from the biological sequence definition, the term '_enumeration_' is used
    in this article.
  ]
  $(mono("TAT"), mono("ATG"), mono("TGC"))$ or $(mono("ATG"), mono("TGG"))$, respectively.
  When comparing the enumerations, we find a match in the common #seq[ATG] sequence.
]

To accelerate the identification of #kmer matches, the obtained #kmers are often stored in a
_hash table_ for fast lookup of the positions where a given #kmer appears.
This requires #kmer _hashing_, i.e. mapping a #kmer to an integer, the _hash_.
If this mapping is unambiguous, i.e. two different #kmers are guaranteed to get different
integers assigned, it is called _perfect hashing_.

Although #kmer decomposition is the fundament of many modern alignment tools, most literature
about their underlying algorithms focus on how the #kmers are used to find the alignment.
This article puts the spotlight on #kmer decomposition itself: It presents an intuitive perfect
hash of a #kmer and describes a fast algorithm that decomposes a sequence into these hash values
in constant time with respect to $k$.
Finally this paper proposes a simple way to give this #kmer hashes a pseudorandom ordering, a
desirable property for certain #kmer based methods, such as _minimizers_ #todo(cite) and
_syncmers_ #todo(cite).

= Methods

== Sequence representation
Each sequence is subject to an _alphabet_ $Omega$, that enumerates the symbols that are allowed in
a type of sequence.
Let the _symbol code_ be the 0-based position of a symbol in $Omega$.
Let the _sequence code_ be the symbol codes for all symbols in the sequence.

#example[
  Take the DNA sequence $s = mono("TATGC")$.
  The underlying alphabet comprises the nucleotides, hence
  $Omega = (mono("A"), mono("C"), mono("G"), mono("T"))$.
  The symbol code for the first symbol in the sequence (#seq[T]), is the 0-based position of it
  in $Omega$ (3).
  Doing this encoding for all symbols in $s$ yields the sequence code $(3, 0, 3, 2, 1)$
]

Note that this integer mapping of a symbol generalizes the definition of a sequence beyond simple
text: If $Omega$ does contain other objects than characters, each enumeration of these objects
can be considered a sequence.

Specifically, in case that the sequence is represented as ASCII text
#footnote[This is true for almost every sequence encountered in biology.]
mapping a sequence to sequence code can be implemented using fast array access
#todo[explain conversion with pseudocode].

== #kmer representation
The aim of the method presented in this article is to represent each #kmer unambigously as single
integer.
Analogous to the symbol code, this integer will be called _#kmer code_.
First, the #kmer is converted into its sequence code as explained above.
Next, the length of $Omega$, written as $|Omega|$, is used as radix to compute the #kmer code $c_k$
from the sequence code $c$ via #footnote[0-based indexing]

$ c_k = sum_(i=1)^k c(i-1) times |Omega|^(k-i). $ <equation_kmer_code>

#example[
  Take the $3$-mer #seq[ATG] that uses again $Omega = (mono("A"), mono("C"), mono("G"), mono("T"))$
  as base alphabet:
  The sequence code of this #kmer is $(0, 3, 2)$.
  The corresponding #kmer code calculates as $c_k = 0 times 4^2 + 3 times 4^1 + 2 times 4^0 = 14$.
]

Note that $c_k$ can be again envisioned as element of an alphabet $Omega_k$ that enumerates all
possible #kmers.
As $Omega_k$ contains every combination of $|Omega|$ symbols in each of its $k$ positions,
the length of such alphabet is

$ |Omega_k| = |Omega|^k. $

== #kmer decomposition
Performing #kmer decomposition of a seqeunce into #kmer codes requires application of
@equation_kmer_code for each overlapping #kmer. Thus,

$ c_k (j) = sum_(i=1)^k c(i+j-1) times |Omega|^(k-i). $ <equation_naive>

A naive implementation of this formula has a time complexity of $O(n k)$, where $n$ is the length of
the sequence.
However, it ignores the relation between two consecutive #kmer codes:

$ c_k (j+1)
  &= sum_(i=1)^k c(i+j) times |Omega|^(k-i) \
  &= lr([ sum_(i=1)^(k-1) c(i+j) times |Omega|^(k-i)  ]) + c(k+j) |Omega|^(k-k) \
  &= lr([ sum_(i=2)^(k) c(i+j-1) times |Omega|^(k-i+1) ]) + c(k+j) \
  &= |Omega| lr([ sum_(i=2)^(k) c(i+j-1) times |Omega|^(k-i) ]) + c(k+j) \
  // &= |Omega| lr([ lr([ sum_(i=1)^(k) c(i+j-1) times |Omega|^(k-i) ]) - c(j) |Omega|^(k-1) ]) \+ c(k+j) \
  &= |Omega| lr([ c_k (j) - c(j) |Omega|^(k-1) ]) + c(k+j). $ <equation_decomposition>

Intuitively, the #kmer code of the previous #kmer is taken, the symbol code of its first symbol
is removed, the remainder is shifted to the left and the symbol code of the entering symbol is
added.

As @equation_decomposition has no sum anymore, the time complexity is reduced to $O(n)$.
Only $c_k (0)$ needs to be computed according to @equation_naive.
In the implementation of @equation_decomposition further speedup can be achieved if $|Omega|$ is a
power of two.
This is true e.g. for unambigous nucleotide sequences with $|Omega| = 4$.
In this case the multiplication with $|Omega|$ can be substituted with a fast bit shift operation
#todo[cite].

== Pseudorandom ordering
In some scenarios the order of #kmer codes is relevant.
For example minimizers #todo(cite) select only the smallest #kmer from a window of #kmers.
This decreases the number of considered #kmers and in consequence improves the speed of finding
#kmer matches between sequences.
However, using the #kmer code directly to determine the ordering is not ideal:
Especially, if the symbols in $Omega$ are ordered alphabetically, the order in $Omega_k$ is
lexicographic.
This means that #kmers with low complexity such as #seq[AAAAA...] would always be the smallest #kmer
leading to more spurious than significant #kmer matches downstream #todo(cite).
A simple way to remedy this behavior is to apply a pseudorandom ordering to the #kmers #todo(cite).

This can for example be achieved by choosing a appropriate #kmer hashing function.
However, in the presented case a #kmer is already represented as integer, the #kmer code.
Therefore, only an injective #footnote[one-to-one] function $f$ is required to obtain an integer
defining the ordering for a #kmer code.
Also I argue, that quality of randomness #todo[find better term] is less important than the speed
of computation for the use case of alignments.
A _linear congruential generator_ (LCG) with _full period_ is suitable in this scenario.
#todo[
  More elaboration.
  Usually used for generate next random value in series, here for random mapping.
  Fast implementation using bit truncation.
]

$ f(x) = (a x + c) mod m $

= Results and Discussion

== Decomposition performance
#figure(
  image("benchmark/benchmark.svg", width: 100%),
  caption: [
    This is a caption.
  ],
) <figure_benchmark>

@figure_benchmark shows how the mentioned decomposition methods compare to each other.

This makes the decomposition especially fast for methods that use long #kmers #todo[cite, kallisto].
#todo[Name speedup at k=31]

Note that the shown benchmark also include sequence encoding itself:
If the implementing library already store sequences in their sequence code form, #kmer decomposition
becomes faster that shown in the benchmark.

= Conclusion
This code representation of sequences and #kmers as well as the fast decomposition method is
implemented in the _Python_ bioinformatics library _Biotite_.