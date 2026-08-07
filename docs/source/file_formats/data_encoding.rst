.. _zero_one_two_etc:

Allele and genotype encoding
============================

|Software| works with bi-allelic loci;
*reference (major)* allele ``a`` and
*alternative (minor)* allele ``A``.
The terms major and minor indicates their frequency in a population,
but this is not a requirement for |Software|.
These alleles are respectively numerically encoded as ``0`` and ``1``.

Combining two alleles in a diploid individual gives three possible *genotypes*:

* reference homozygote ``a/a`` encoded as ``0`` (from ``0+0``),
* heterozygote ``a/A`` or ``A/a`` encoded as ``1`` (from ``0+1`` or ``1+0``), and
* alternative homozygote ``A/A`` encoded as ``2`` (from ``1+1``).

The codes ``0`` and ``1`` for alleles, and ``0``, ``1``, and ``2`` for genotypes
represent the number (dose) of alternative alleles and
are hence called *allele dosage* or just *dosage*.

Missing or uncertain alleles and genotypes are numerically encoded as ``9``.

When the origin of alleles that an individual inherited is known,
we have four *phased genotypes*: ``aa``, ``aA``, ``Aa``, and ``AA``,
where the *paternal allele* is listed first and
the *maternal allele* is listed second.

Related to phased genotypes is the concept of a *haplotype*,
which represents a sequence of alleles an individual inherited together from one parent.
For example, individual with genotypes at three loci of ``a/a``, ``a/a``, and ``A/A``
has haplotypes ``[a,a,A]`` and ``[a,a,A]``,
which we numerically encode as ``[0,0,1]`` and ``[0,0,1]``.

When working with the X chromosome:

(1) *heterogametic genotypes* (for males in mammals (XY) and for females in birds (ZW))
should be coded as:
``0`` (reference allele ``a`` on the X chromosome of the XY genotype),
``1`` (alternative allele ``A`` on the X chromosome of the XY genotype), or
``9`` (missing) and

(2) *Homogametic genotypes* (females in mammals (XX) and males in birds (ZZ))
should be coded as described above for autosomes
(since they have the XX genotype).
