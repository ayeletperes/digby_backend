"""Schema for gene-usage QTL (guQTL) results.

One database per species and locus, matching the layout the genomic and AIRR-seq
databases already use, so the same discovery walks it. Each database holds a
single analysis run; `Run` records which one, and no other table carries a run id.
Re-running the analysis means rebuilding the file.

The shape follows the output of lib_manuscript/analysis/igqtl.R, documented in
husa_manuscript/docs/analyses/igqtl.md. Two things about that data drive the
choices here:

- p-values reach 1e-53, so they need double precision, and `neglog10_p` is stored
  alongside because every plot and every ordering wants it rather than the raw
  value.
- positions are locus-relative. IGH sits on a contig literally named `igh` while
  IGK and IGL use chr2 and chr22, so a position is only meaningful within its
  locus - which is another reason one database per locus.
"""

from sqlalchemy import (
    Boolean, Column, Float, ForeignKey, Index, Integer, String, Text,
    UniqueConstraint,
)
from sqlalchemy.ext.declarative import declarative_base

from db.db_propertymixin import Details_Mixin

Base = declarative_base()
metadata = Base.metadata


class Details(Details_Mixin, Base):
    __tablename__ = "details"


class Run(Base):
    """Provenance for one analysis run: one study's scan at this locus.

    A database holds every study built at its locus, one run row each. Nothing
    is ever compared across them - the scans have different cohorts, different
    thresholds and therefore p-values on different scales - so `run_id` is a
    filter every view applies, not a dimension anything aggregates over.
    """
    __tablename__ = 'qtl_run'

    id = Column(Integer, primary_key=True)
    label = Column(String(100))
    generated_at = Column(String(40))
    script = Column(String(500))
    # Which study's cohort was scanned. A scan is computed within a cohort and
    # never pooled across them, so this qualifies every query rather than
    # selecting a file, the same way study_id does in the repertoire schema.
    #
    # Stated by whoever builds the database, not derived. The run's own config
    # names a metadata file whose stem happens to carry the project, and reading
    # a name for an identity is how a wrong answer gets in - the same mistake as
    # taking a gene's segment from the fourth letter of its name. Null where a
    # database predates the column; the API reports that as unknown rather than
    # guessing.
    project = Column(String(60))
    # the analysis config verbatim, as JSON: every filter and its value, which is
    # what makes a figure reproducible
    config = Column(Text)


class Threshold(Base):
    """Significance threshold and cohort size, per analysis.

    The usage scan and the two pairing scans each have their own threshold, set
    from the number of independent variants after LD collapse rather than the
    number tested.
    """
    __tablename__ = 'qtl_threshold'

    run_id = Column(Integer, ForeignKey('qtl_run.id'), nullable=False, index=True)
    id = Column(Integer, primary_key=True)
    analysis = Column(String(20), nullable=False)        # usage | pairing
    # pairing only: which side the scan was anchored on, and P(J|D) or P(D|J)
    grouped_by = Column(String(20))
    conditional = Column(String(20))

    n_subjects = Column(Integer)
    n_excluded = Column(Integer)
    n_variants = Column(Integer)
    n_independent = Column(Integer)
    n_asc = Column(Integer)
    threshold = Column(Float, nullable=False)
    n_significant_variants = Column(Integer)
    n_independent_significant = Column(Integer)

    __table_args__ = (
        UniqueConstraint('run_id', 'analysis', 'conditional', name='uq_threshold_analysis'),
    )


class Subject(Base):
    """A genotyped subject.

    The identifier is opaque on purpose: this cohort mixes numeric Watson ids
    (110040277) with VDJbase-style names (P28_I61_S1) in the same column.
    """
    __tablename__ = 'qtl_subject'

    run_id = Column(Integer, ForeignKey('qtl_run.id'), nullable=False, index=True)
    id = Column(Integer, primary_key=True)
    subject = Column(String(60), nullable=False)
    ancestry = Column(String(40))


class Variant(Base):
    """A tested variant, with the genomic feature it falls in.

    `contig` and `pos` are locus-relative; see the module docstring.
    """
    __tablename__ = 'qtl_variant'

    run_id = Column(Integer, ForeignKey('qtl_run.id'), nullable=False, index=True)
    id = Column(Integer, primary_key=True)
    variant = Column(String(60), nullable=False)
    contig = Column(String(20))
    pos = Column(Integer)
    maf = Column(Float)
    missing_rate = Column(Float)
    # the scan's subject count and smallest genotype class for this variant, and
    # whether that class is big enough to trust the fit. Per variant, not per
    # association: they were on every one of the 667,542 IGH association rows
    n = Column(Integer)
    min_genotype_group = Column(Integer)
    well_powered = Column(Boolean)

    # from variant_features.tsv; populated for IGH, where the association rows
    # themselves carry no gene context
    gene = Column(String(60))
    feature = Column(String(40))
    sub_feature = Column(String(40))
    distance_to_gene = Column(Float)

    # The strongest association this variant showed against any ASC, and which
    # ASC that was. Derived from qtl_usage_association and stored beside what it
    # summarises: the whole-locus Manhattan is one point per variant, and finding
    # 9,402 maxima in 658,140 rows took the endpoint several seconds on every
    # request for a database that never changes after it is built. Computed once
    # at build time instead, which is 50 ms there and 23 ms to read back.
    #
    # Null on a database built before these existed; the API notices and computes
    # the aggregate live rather than reporting a variant as untested.
    best_neglog10_p = Column(Float)
    best_asc_id = Column(Integer, ForeignKey('qtl_asc.id'))
    best_significant = Column(Boolean)

    __table_args__ = (
        Index('ix_qtl_variant_pos', 'pos'),
    )


class Asc(Base):
    """An allele similarity cluster: the feature whose usage is being explained.

    Names are not consistent across segments in the source data - D clusters are
    written in full (IGHD5-12) where V and J are abbreviated (V1-18, J4) - so the
    name is stored as given and `segment` carries the distinction.
    """
    __tablename__ = 'qtl_asc'

    run_id = Column(Integer, ForeignKey('qtl_run.id'), nullable=False, index=True)
    id = Column(Integer, primary_key=True)
    asc = Column(String(60), nullable=False)
    segment = Column(String(2))                          # V | D | J

    # empty for IGH, where gene coordinates were not resolved
    asc_position = Column(Float)
    asc_span = Column(Float)
    n_member = Column(Integer)


class UsageAssociation(Base):
    """One variant tested against one ASC's usage. The Manhattan plot reads this.

    Effect direction is relative to the genotype matrix coding, which is not
    always the minor allele, so `beta` must not be described in terms of a named
    allele without resolving polarity first.
    """
    __tablename__ = 'qtl_usage_association'

    id = Column(Integer, primary_key=True)
    variant_id = Column(Integer, ForeignKey('qtl_variant.id'), nullable=False)
    asc_id = Column(Integer, ForeignKey('qtl_asc.id'), nullable=False)

    beta = Column(Float)
    se = Column(Float)
    t_stat = Column(Float)
    p_value = Column(Float, nullable=False)
    # -log10(p), stored because every plot and sort wants it and the raw value
    # reaches 1e-53
    neglog10_p = Column(Float, nullable=False)

    significant = Column(Boolean)
    distance_to_asc = Column(Float)

    # is_cis is per variant AND cluster - it is the distance to that cluster -
    # so unlike the power columns, which moved to qtl_variant, it belongs here
    is_cis = Column(Boolean)
    is_lead = Column(Boolean, default=False)

    __table_args__ = (
        # variant_id needs no index of its own: it leads the unique constraint
        # above, so SQLite uses that index for a lookup on it
        UniqueConstraint('variant_id', 'asc_id', name='uq_usage_variant_asc'),
        Index('ix_usage_asc_p', 'asc_id', 'neglog10_p'),
        Index('ix_usage_significant', 'significant', 'neglog10_p'),
    )


class AscUsage(Base):
    """One subject's usage of one ASC: the phenotype the scan explains.

    `usage` is the share of that subject's repertoire for the ASC's own segment,
    with a pseudocount; `logit_usage` is what was actually regressed. Zero counts
    are kept deliberately - a deleting variant driving usage to zero is the
    signal, not missing data.
    """
    __tablename__ = 'qtl_asc_usage'

    id = Column(Integer, primary_key=True)
    subject_id = Column(Integer, ForeignKey('qtl_subject.id'), nullable=False)
    asc_id = Column(Integer, ForeignKey('qtl_asc.id'), nullable=False)

    count = Column(Integer)
    total = Column(Integer)
    # how many ASCs of this segment the share was taken over: part of the
    # pseudocount, so `usage` cannot be recomputed without it
    n_asc = Column(Integer)
    usage = Column(Float)
    logit_usage = Column(Float)

    __table_args__ = (
        UniqueConstraint('subject_id', 'asc_id', name='uq_asc_usage_subject_asc'),
        Index('ix_asc_usage_asc', 'asc_id'),
    )


class Dosage(Base):
    """A subject's genotype at a variant.

    `dosage` may be non-integer where genotypes were imputed; `genotype` is it
    rounded to 0, 1 or 2, which is what the boxplot groups on.
    """
    __tablename__ = 'qtl_dosage'

    id = Column(Integer, primary_key=True)
    variant_id = Column(Integer, ForeignKey('qtl_variant.id'), nullable=False)
    subject_id = Column(Integer, ForeignKey('qtl_subject.id'), nullable=False)

    dosage = Column(Float)
    genotype = Column(Integer)

    __table_args__ = (
        # as above: variant_id leads the unique constraint, so it is covered
        UniqueConstraint('variant_id', 'subject_id', name='uq_dosage_variant_subject'),
    )


class PairingAssociation(Base):
    """A variant tested against a gene's partner distribution, by MANOVA.

    `conditional` is P(J|D) or P(D|J). Both scans cover the same variants, so it
    has to be a filter everywhere or the counts double-report.
    """
    __tablename__ = 'qtl_pairing_association'

    id = Column(Integer, primary_key=True)
    conditional = Column(String(20), nullable=False)
    variant_id = Column(Integer, ForeignKey('qtl_variant.id'), nullable=False)
    anchor_gene = Column(String(40), nullable=False)

    n = Column(Integer)
    pillai = Column(Float)
    f_stat = Column(Float)
    p_value = Column(Float, nullable=False)
    neglog10_p = Column(Float, nullable=False)
    min_genotype_group = Column(Integer)
    significant = Column(Boolean)

    __table_args__ = (
        UniqueConstraint('conditional', 'variant_id', 'anchor_gene',
                         name='uq_pairing_conditional_variant_anchor'),
        Index('ix_pairing_anchor_p', 'conditional', 'anchor_gene', 'neglog10_p'),
    )


class CellTest(Base):
    """Post-hoc test of one D/J cell under a pairing hit.

    A cell is only meaningful inside a row the omnibus already called
    significant, which `omnibus_significant` records.
    """
    __tablename__ = 'qtl_cell_test'

    id = Column(Integer, primary_key=True)
    conditional = Column(String(20), nullable=False)
    variant_id = Column(Integer, ForeignKey('qtl_variant.id'), nullable=False)
    d_gene = Column(String(40), nullable=False)
    j_gene = Column(String(40), nullable=False)

    n = Column(Integer)
    beta = Column(Float)
    p_value = Column(Float)
    mean_low = Column(Float)
    mean_high = Column(Float)
    delta_mean = Column(Float)
    # how many subjects sat either side of the split, which is what says whether
    # a difference in means is worth anything
    n_low = Column(Integer)
    n_high = Column(Integer)
    omnibus_p_value = Column(Float)
    omnibus_significant = Column(Boolean)
    min_genotype_group = Column(Integer)
    # the pipeline's own call on whether to star this cell; `marked_strict`
    # applies the tighter rule. Kept rather than re-derived so the dashboard and
    # the manuscript figures agree.
    marked = Column(Boolean)
    marked_strict = Column(Boolean)

    __table_args__ = (
        UniqueConstraint('conditional', 'variant_id', 'd_gene', 'j_gene',
                         name='uq_cell_conditional_variant_pair'),
        Index('ix_cell_variant', 'variant_id'),
    )


class DjEnrichment(Base):
    """One subject's observed D/J pairing against what independence would give."""
    __tablename__ = 'qtl_dj_enrichment'

    id = Column(Integer, primary_key=True)
    subject_id = Column(Integer, ForeignKey('qtl_subject.id'), nullable=False)
    d_gene = Column(String(40), nullable=False)
    j_gene = Column(String(40), nullable=False)

    count = Column(Integer)
    depth = Column(Integer)
    # the marginals the expectation is built from; kept so a ratio can be checked
    d_total = Column(Integer)
    j_total = Column(Integer)
    expected = Column(Float)
    enrichment = Column(Float)
    p_d = Column(Float)
    p_j = Column(Float)
    p_j_given_d = Column(Float)
    p_d_given_j = Column(Float)

    __table_args__ = (
        UniqueConstraint('subject_id', 'd_gene', 'j_gene',
                         name='uq_dj_subject_pair'),
    )
