#import "ieee_template.typ": ieee
#import "@preview/algo:0.3.3": algo, i, d, comment


#let kmer = box[$k$-mer]
#let kmers = box[$k$-mers]

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

#show figure.caption: set align(left)
// Forbid line breaks in equations
#show math.equation: it => box[#it]


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

To solve this problem heuristic approaches emerged. Many modern algorithms
(see @Bray2016 @Steinegger2017 as examples) build upon the concept of finding exact matches of
length $k$ between the sequences @Altschul1990.
These subsequences of length $k$ are termed #kmers.
The process of finding all overlapping #kmers in a sequence is commonly titled
_k-mer decomposition_.

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
desirable property for certain #kmer based methods, such as _minimizers_ @Roberts2004 and
_syncmers_ @Edgar2021.

The techniques presented in this article are also implemented in the _Python_ bioinformatics library
_Biotite_ @Kunzmann2023.

= Algorithm

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

Specifically, in case that the sequence is represented as ASCII text
#footnote[This is true for almost every sequence type encountered in biology.]
mapping a sequence $S$ to sequence code can be implemented using fast array access as described in
@figure_encode, in contrast to employing a slower associative array
#footnote[This container type is also termed _map_ or _dictionary_ depending on the programming
language.].

#figure(
  algo(
    title: [#smallcaps("Encode")],
    parameters: ($S$,$Omega$),
    comment-prefix: [#sym.triangle.stroked.r ],
    indent-size: 15pt,
    indent-guides: 1pt + gray,
    row-gutter: 10pt,
    column-gutter: 5pt,
    inset: 5pt,
  )[
    $s_("illegal") <- |Omega|$ \
    $m_("symbol"->"code") <- "repeat"(s_("illegal"), 256)$\

    for $i$, $s$ in enumerate($Omega$): #i \
      $m_("symbol"->"code")["as_ascii"(s)] <- i$ #d \

    $C <- "repeat"(0, |S|)$ \
    for $i$, $s$ in enumerate($S$): #i \
      $C[i] <- m_("symbol"->"code")["as_ascii"(s)]$ #d \

    return $C$
  ],
  caption: [
    Sequence encoding into sequence code.
    The input sequence $S$ is subject to alphabet $Omega$, which contains only ASCII-encoded
    symbols.
    As symbol codes are 0-based, the greatest symbol code is $|Omega| - 1$.
    Hence, $s_("illegal") <- |Omega|$ can be used as marker value to check for symbols that are
    included ib not in $Omega$.
    The array $m_("symbol"->"code")$ can be indexed with the ASCII code of a symbol to obtain the
    corresponding symbol code.
    For symbols that are not part of $Omega$, $m_("symbol"->"code")$ would give $s_("illegal")$.
  ],
) <figure_encode>

== #kmer representation
The aim of the method presented in this article is to represent each #kmer unambigously as single
integer.
Analogous to the symbol code, this integer will be called _#kmer code_.
First, the #kmer is converted into its sequence code as explained above.
Next, the length of $Omega$, written as $|Omega|$, is used as radix to compute the #kmer code $c_k$
from the sequence code $c$ via

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

The sequence code $c$ can also be restored from the #kmer code $c_k$ via

$ c(i) = (c_k div |Omega|^(k-i)) mod |Omega|, $ <equation_kmer_decode>

where '$div$' denotes integer division.

== #kmer decomposition
Performing #kmer decomposition of a seqeunce into #kmer codes requires application of
@equation_kmer_code for each overlapping #kmer. Thus,

$ c_k (j) = sum_(i=1)^k c(i+j-1) times |Omega|^(k-i), $ <equation_naive>

where $j$ defines the 0-based sequence position where the #kmer starts.
A naive implementation of this formula has a time complexity of $O(n k)$, where $n$ is the length of
the sequence.
However, it ignores the relation between two consecutive #kmer codes:

$ c_k (j+1)
  &= sum_(i=1)^k c(i+j) times |Omega|^(k-i) \
  &= lr([ sum_(i=1)^(k-1) c(i+j) times |Omega|^(k-i)  ]) + c(k+j) |Omega|^(k-k) \
  &= lr([ sum_(i=2)^(k) c(i+j-1) times |Omega|^(k-i+1) ]) + c(k+j) \
  &= |Omega| lr([ sum_(i=2)^(k) c(i+j-1) times |Omega|^(k-i) ]) + c(k+j) \
  //&= |Omega| lr([ lr([ sum_(i=1)^(k) c(i+j-1) times |Omega|^(k-i) ]) - c(j) |Omega|^(k-1) ]) \+ c(k+j) \
  &= |Omega| lr([ c_k (j) - c(j) |Omega|^(k-1) ]) + c(k+j).
$ <equation_decomposition>

Intuitively, the #kmer code of the previous #kmer is taken, the symbol code of its first symbol
is removed, the remainder is shifted to the left and the symbol code of the entering symbol is
added.

As @equation_decomposition contains no sum anymore, the time complexity is reduced to $O(n)$.
Only $c_k (0)$ needs to be computed according to @equation_naive.
In the implementation of @equation_decomposition potentially further speedup can be achieved if
$|Omega|$ is a power of two.
This is true e.g. for unambigous nucleotide sequences with $|Omega| = 4$.
In this case the compiler may substitute this multiplication with a fast bit shift operation
depending on the hardware architecture.

Note that @equation_decomposition only works for linear #kmers.
Some algorithms in homology search use _spaced_ #kmers @Ma2002, which contains ignored positions.
In this case @equation_naive can be still used.

== Pseudorandom ordering
In some scenarios the order of #kmer codes is relevant.
For example minimizers @Roberts2004 select only the smallest #kmer from a window of #kmers.
This decreases the number of considered #kmers and in consequence improves the speed of finding
#kmer matches between sequences.
However, using the #kmer code directly to determine the ordering is not ideal:
Especially, if the symbols in $Omega$ are ordered alphabetically, the order in $Omega_k$ is
lexicographic.
This means that #kmers with low complexity such as #seq[AAAAA...] would always be the smallest #kmer
leading to more spurious than significant #kmer matches downstream @Roberts2004.
A simple way to remedy this behavior is to apply a pseudorandom ordering to the #kmers #todo(cite).
This can for example be achieved by choosing a appropriate #kmer hashing function
#todo[cite review with hash function].

However, in the presented case a #kmer is already represented as integer: the #kmer code.
Therefore, only an injective #footnote[one-to-one] function $sigma(c_k)$ is required to obtain an
integer defining the ordering for a #kmer $c_k$, i.e. the #kmer code $p$ is defined to be smaller
than $q$ if $sigma(p) lt sigma(q)$.
Furthermore, for the use case of computing sequence alignments, the quality of randomness
is arguably less important than the speed of computation.
Hence, a _linear congruential generator_ (LCG) is appealing in this case.
It predicts the next random number in a sequence via

$ x_(n+1) = (a x_n + b) mod m $ <equation_lcg>

In a LGC with _full period_ the sequence does only repeat after $m$ elements.
To achieve the full period attention has to be paid to the choice of $a$ and $m$.
Furthermore, $b$ and $M$ need to be coprime, which can be trivially achieved by setting $b=1$
@Tezuka1995.

For the purpose of #kmer ordering, the LCG should be used to map a #kmer code $c_k$ to a unique
pseudorandom value $sigma(c_k)$ that defines the #kmer ordering.
Thus one can employ @equation_lcg to define

$ sigma(c_k) = (a c_k + b) mod m. $ <equation_lcg_kmer>

$sigma(c_k)$ is only injective, if each $c_k$ is mapped to a unique value.
To ensure this property, an LGC with full period is used.
If one represents #kmer codes as 64-bit integers
#footnote[$|Omega|^k$ quickly leads to an combinatorial explosion of $|Omega_k|$, making 64-bit integers necessary.]
$m = 2^(64)$ is required.
When carefully implemented, the modulo computation in @equation_lcg is free due to automatic bit
truncation.
For $a$ one can resort to published multipliers @Steele2021 that ensure both, full periodicity and
good randomness for the chosen $m$.
As example the following combination fulfills the requirements:

$ a = "d1342543de82ef95"_16 \
  b = 1 \
  m = 2^64
$

#example[
  Take the nucleotide $3$-mers #seq[ATG] and #seq[TGG] with corresponding #kmer codes $p=14$ and
  $q=58$, respectively.
  $sigma(p) = 8131822755183663655$ and $sigma(q) = 7336488451890104259$, i.e $sigma(p) gt sigma(q)$.
  Thus, by the newly defined order, the $3$-mer #seq[TGG] is smaller than #seq[ATG] in contrast
  to their lexicographic order.
]

= Results and Discussion

The presented decomposition methods were benchmarked for different #kmer lengths on a 1000~bp
nucleotide sequence as shown in @figure_benchmark.
The benchmark was run on a system with _Apple M3 Max_ processor \@ 4.05 GHz using an implementation
written in _Rust_ (Supplementary File 1).
As expected, the naive method scales linearly with $k$
($T approx (0.97 + 0.36 k) thin mu s$, $R^2=0.9997$).
In contrast, the fast decomposition method based on @equation_decomposition runs in constant time
($T approx 1.55 thin mu s$).
As such potential optimization is hardware-related, this result may depend on the actual
architecture and programming language.

#figure(
  image("benchmark/benchmark.svg", width: 100%),
  caption: [
    Run time of #kmer decomposition using different methods.
    Decomposition was run on a sequence with length 1000.
    The displayed run time includes also the consersion into sequence code.
    *naive*: Naive application of @equation_naive for each sequence position.
    *fast*: Application of @equation_decomposition.
  ]
) <figure_benchmark>

In summary the fast decomposition method is already faster than the naive method for
$k gt.eq 2$, i.e. for any practically relevant #kmer length.
The fast method is especially advantageous for algorithms that utilize long #kmers.
For example, by default _Minimap~2_ @Li2018 uses $k=15$ and _Kallisto_ @Bray2016 uses $k=31$.
For these examples, the fast decomposition method is $tilde.op #h(0cm) 4 #h(0cm) times$ and
$tilde.op #h(0cm) 8 #h(0cm) times$ faster than the naive method, respectively.

Note that the implementation used for the benchmark also includes sequence encoding itself:
If the implementing library already stores sequences in their sequence code form, #kmer
decomposition becomes faster than shown in the benchmark.

= Conclusion

This article advocates representing #kmers as integer in memory for mainly two reasons:
First, it reduces the time complexity of #kmer decomposition to $O(n)$.
Since modern sequence alignment algorithms strive to be as fast as possible, this performance gain
may pose a crucial advantage.
Second, many current application of #kmers already implicitly rely on conversion of #kmers
into an integer by means of hashing.
Among other applications, this includes

- comparison of two sets of #kmers to approximate sequence identity (e.g. @Edgar2004),
- #todo[another example] and
- efficiently finding match positions between two sequences (e.g. @Steinegger2017).

Already having #kmers as unique integers at hand removes the need for hashing them and thus may
further speeds up such applications.

In addition, representing a sequence as array of integers, has the advantage of generalizing the
definition of a sequence beyond simple text:
If the alphabet $Omega$ does contain other objects than single characters as symbols, e.g. arbitrary
strings or integers, each enumeration of these objects can be considered a sequence.
This allows alphabets to contain more symbols than the 95 printable ASCII characters, which would
enable, for example, creating and representing more fine-grained structural alphabets
@Brevern2000 @VanKempen2024 @Wang2008 in the future.