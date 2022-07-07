# Configuration and topology files

The current state of a system, as specified by oxDNA, is described by two files: a configuration file and a topology file. The configuration file contains all the general information (timestep, energy and box size) and the orientations and positions of each nucleotide. The topology file, on the other hand, keeps track of the backbone-backbone bonds between nucleotides in the same strand. Working configuration and topology files can be found in the `examples` directory.

## Configuration file

The first three rows of a configuration file contain the timestep T at which the configuration has been printed, the length of the box sides Lx, Ly and Lz and the total, potential and kinetic energies, Etot, U and K, respectively:

```text
t = T
b = Lz Ly Lz
E = Etot U K
```

After this header, each row contains position of the centre of mass, orientation, velocity and angular velocity of a single nucleotide in the following order: 

$$
\overbrace{r_x r_y r_z}^{\rm centre-of-mass\,\, position\,\, r} \underbrace{b_x b_y b_z}_{{\rm base\,\, vector\,\,} \vec{a}_1} \overbrace{n_x n_y n_z}^{{\rm base\,\, normal\,\, vector\,\,} \vec{a}_3} \underbrace{v_x v_y v_z}_{\rm Velocity} \overbrace{L_x L_y L_z}^{\rm Angular\,\, velocity}
$$

$\vec{a}_1$, $\vec{a}_2 = \vec{a}_3 \times \vec{a}_1$ and $\vec{a}_3$ define the local reference frame through which the position of all interaction sites relative to the centre of mass are calculated. The position of the interaction sites can be recovered as follows:

* for oxDNA1:
  * hydrogen-bonding/repulsion site $ = \vec{r} + 0.4 \, \vec{a}_1$;
  * stacking site $ = \vec{r} + 0.34 \, \vec{a}_1$;
  * backbone repulsion site $ = \vec{r} - 0.4 \, \vec{a}_1$.
* For oxDNA2:
  * hydrogen-bonding site $ = \vec{r} + 0.4 \, \vec{a}_1$;
  * stacking site $ = \vec{r} + 0.34 \, \vec{a}_1$;
  * backbone repulsion site $ = \vec{r} - 0.34 \, \vec{a}_1 + 0.3408 \, \vec{a}_2$.
  
```{note}
When simulating very large structures the size of the trajectory files stored on disk may become very large. In this case it may be convenient to avoid printing the last six columns by setting `trajectory_print_momenta = false` in the input file, thus decreasing the size of the trajectory by $\approx 40\%$. 
```

## Topology file

The topology file stores the intra-strand, fixed bonding topology (*i.e.* which nucleotides share backbone links). The first row contains the total number of nucleotides N and the number of strands Ns: 

```text
N Ns
```

After this header, the *i*-th row specifies strand, base and 3' and 5' neighbours of the *i*-th nucleotide in this way:

```text
S B 3' 5'
```

where S is the index of the strand (starting from 1) which the nucleotide belongs to, B is the base (*A*, *C*, *G*, and *T* for DNA or *A*, *C*, *G*, and *U* for RNA, but see below for more options) and 3' and 5' specify the index of the nucleotides with which the *i*-th nucleotide is bonded in the 3' and 5' direction, respectively. A -1 signals that the nucleotide terminates the strand in either 3' or 5' direction. 

```{warning}
OxDNA's convention is to list nucleotides in the 3' {math}`\to` 5' order. Note that this is the opposite of how most of the other DNA-related tools behave.
```

The topology file of a strand of sequence GCGTTG would be:

```text
6 1
1 G -1 1
1 C 0 2
1 G 1 3
1 T 2 4
1 T 3 5
1 G 4 -1
```

Specifying the topology in this way can simplify the process of simulating, for example, circular DNA.

### Special nucleotides

Internally, oxDNA encodes the base types with integers as follows:

```text
A = 0
G = 1
C = 2
T = U = 3
```

According to these values, two nucleotides interact through the Watson-Crick mechanism (*i.e.* can be hydrogen bonded) if the sum of their base types is 3. 

This property can be leveraged to extend the canonical base pairing and create very specific topologies. Indeed, in oxDNA a nucleotide of type B which is larger than 9 or smaller than 0 behaves as a nucleotide of type $B \bmod 4$ if $B$ is positive or $3 - ((3 - B) \bmod 4)$ if $B$ is negative, but can be hydrogen-bonded **only** with a nucleotide of type $B'$ for which $B + B' = 3$.

For instance, $B = 13$, for which $B \bmod 4 = 1$, would correspond to a nucleotide with the same property of a Guanine. However, such a nucleotide would bond only to a nucleotide with base type $B' = -10$.

```{warning}
The CUDA backend supports base types whose absolute values do not exceed $2^{10} - 1 = 511$. In other words, base types larger than $511$ or smaller than $-511$ are not allowed.
```

## Converting to and from oxDNA configurations

Many DNA-nanotechnology tools support exporting to oxDNA configurations, as well as converting between different formats (see *e.g.* [tacoxDNA](https://doi.org/10.1002/jcc.26029), [oxView](https://academic.oup.com/nar/article/48/12/e72/5843822), [Adenita](https://doi.org/10.1093/nar/gkaa593), [MrDNA](https://doi.org/10.1093/nar/gkaa200)).
