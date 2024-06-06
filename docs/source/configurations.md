# Configuration and topology files

The current state of a system, as specified by oxDNA, is described by two files: a configuration file and a topology file. The configuration file contains all the general information (timestep, energy and box size) and the orientations and positions of each nucleotide. The topology file, on the other hand, keeps track of the backbone-backbone bonds between nucleotides in the same strand. Working configuration and topology files can be found in the `examples` directory.

## Configuration file

```{warning}
Nucleotides are listed in the 3' {math}`\to` 5' order if using the "classic" topology format or in the 5' {math}`\to` 3' order if using the "new" format (see [Topology file](#topology-file)).
```

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
  
```{warning}
The position of the centre of mass of the nucleotides in oxDNA1 (0.4 length units away from the backbone site) is different from what the PhD thesis of T. E. Ouldridge specifies (0.24 length units away from the backbone site). This change has no effect on the thermodynamics, and the extent to which it changes the dynamics is arguably very small.
```
  
```{note}
When simulating very large structures the size of the trajectory files stored on disk may become very large. In this case it may be convenient to avoid printing the last six columns by setting `trajectory_print_momenta = false` in the input file, thus decreasing the size of the trajectory by $\approx 40\%$. 
```

## Topology file

The topology file stores information about the fixed-bonding topology (*i.e.* which nucleotides share backbone links), as well as the strand sequences. Along with the original topology format (the "classic" one which is in the 3'-5' direction), starting with version 3.6, oxDNA also supports a new format that is simpler and more flexible (and in 5'-3' direction). Note that this "new" format is not necessarily supported by the other tools of the oxDNA ecosystem. 

```{note}
You can interconvert between the classic and new formats by using the `utils/convert.py` script.
```

### Classic format (3' {math}`\to` 5')

The first row contains the total number of nucleotides N and the number of strands Ns: 

```text
N Ns
```

After this header, the *i*-th row specifies strand, base and 3' and 5' neighbours of the *i*-th nucleotide in this way:

```text
S B 3' 5'
```

where S is the index of the strand (starting from 1) which the nucleotide belongs to, B is the base (*A*, *C*, *G*, and *T* for DNA or *A*, *C*, *G*, and *U* for RNA, but see below for more options) and 3' and 5' specify the index of the nucleotides with which the *i*-th nucleotide is bonded in the 3' and 5' direction, respectively. A -1 signals that the nucleotide terminates the strand in either 3' or 5' direction. 

```{warning}
OxDNA's "classic topology" convention is to list nucleotides in the 3' {math}`\to` 5' order. Note that this is the opposite of how most of the other DNA-related tools behave.
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

### New format (5' {math}`\to` 3')

The new format (introduced in oxDNA 3.6) lists nucleotides in the more common order (5' {math}`\to` 3'), is (arguably) simpler, and is more flexible. However, it is not fully supported by all oxDNA-related tools yet.

The first row of the topology should contain the total number of nucleotides N, the number of strands Ns and the "5->3" string to signal that the file has the "new" format:

```text
N Ns 5->3
```

Each of the Ns rows that follows contains the details of a single strand as a space-separated list of elements. The first element is the sequence of the strand (in 5' {math}`\to` 3' order), while additional elements can be specified with the key=value syntax. The topology of the same system used as an example of "classic" topology (*i.e.* a system composed of a single strand of sequence GTTGCG) is

```text
6 1 5->3
GTTGCG
```

As another example, the topology of a system composed by two complementary strands would be

```
12 2 5->3
GTTGCG
CGCAAC
```

In addition, DNA and RNA interactions also support the `type=DNA|RNA` (which defaults to `DNA`) and `circular=true|false` (which defaults to `false`) specifiers. The former sets the type of strand, while the latter, if set to `true`, indicates that the strand is circular. For example, if you wanted to explicitally note that the above strands are DNA, the file would read

```
12 2 5->3
GTTGCG type=DNA
CGCAAC type=DNA
```

```{note}
Setting `type` only affects the force field when using the DNA/RNA hybrid model (`interaction_type=NA`). For normal DNA or RNA, the interactions will be determined by the `interaction_type=DNA|DNA2|RNA|RNA2|NA|LJ...` parameter in the [input file](input.md)
```

```{note}
The `examples/PERSISTENCE_LENGTH/NEW_TOPOLOGY` directory contains an example that uses this new format.
```

### Special nucleotides

Internally, oxDNA encodes the base types with integers as follows:

```text
A = 0
G = 1
C = 2
T = U = 3
```

According to these values, two nucleotides interact through the Watson-Crick mechanism (*i.e.* can be hydrogen bonded) if the sum of their base types is 3. 

This property can be leveraged to extend the canonical base pairing and create very specific topologies. Indeed, in oxDNA a nucleotide of custom type X which is larger than 9 or smaller than 0 behaves as a nucleotide of type $X \bmod 4$ if $X$ is positive or $3 - ((3 - X) \bmod 4)$ if $X$ is negative, but can be hydrogen-bonded **only** with a nucleotide of type $X'$ for which $X + X' = 3$.

```{note}
In the classic topology format, custom nucleotide types are set by using numbers instead of letters. For instance, the following topology line specifies that the corresponding nucleotide (which is part of strand 1 and bonded to nucleotides 2 and 4) has a custom type `-10`: `1 -10 2 4`.

In the new topology format, custom nucleotide types can be set by enclosing them between brackets. As an example, the following line sets the sequence of a DNA strand made of 6 nucleotides, with the third one having a custom type `-10`: `AA(-10)GCT type=DNA`.
```

For instance, $B = 13$, for which $B \bmod 4 = 1$, would correspond to a nucleotide with the same property of a Guanine. However, such a nucleotide would bond only to a nucleotide with base type $B' = -10$.

```{warning}
The CUDA backend supports base types whose absolute values do not exceed $2^{9} - 1 = 511$. In other words, base types larger than $511$ or smaller than $-511$ are not allowed.
```

## Converting to and from oxDNA configurations

Many DNA-nanotechnology tools support exporting to oxDNA configurations, as well as converting between different formats (see *e.g.* [tacoxDNA](https://doi.org/10.1002/jcc.26029), [oxView](https://academic.oup.com/nar/article/48/12/e72/5843822), [Adenita](https://doi.org/10.1093/nar/gkaa593), [MrDNA](https://doi.org/10.1093/nar/gkaa200)).
