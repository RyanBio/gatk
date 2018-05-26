package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public abstract class BreakEndVariantType extends SvType {

    /**
     * Technically, a BND-formatted variant should have two VCF records, for mates, hence we also have this field.
     * Upstream mate is defined as the location in a mate pair that has a lower coordinate according to
     * the reference sequence dictionary.
     */
    private final boolean isTheUpstreamMate;

    private BreakEndVariantType(final String id, final Allele altAllele, final boolean isTheUpstreamMate,
                                final Map<String, String> typeSpecificExtraAttributes) {
        super(id, altAllele, SVContext.NO_LENGTH, typeSpecificExtraAttributes);
        this.isTheUpstreamMate = isTheUpstreamMate;
    }

    public final boolean isTheUpstreamMate() {
        return isTheUpstreamMate;
    }

    @Override
    public final String toString() {
        return BREAKEND_STR;
    }

    //==================================================================================================================

    private static String getIDString(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
        // if no strand switch or different contig, "", otherwise append INV55/33
        final String bndtype = narl.getStrandSwitch().equals(StrandSwitch.NO_SWITCH) || !narl.getLeftJustifiedLeftRefLoc().getContig().equals(narl.getLeftJustifiedRightRefLoc().getContig())? ""
                : (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE) ? INV55 : INV33);
        String locationPartOfString = makeLocationPartOfID(narl.getLeftJustifiedLeftRefLoc().getContig(),
                narl.getLeftJustifiedLeftRefLoc().getStart(), narl.getLeftJustifiedRightRefLoc().getContig(),
                narl.getLeftJustifiedRightRefLoc().getEnd());
        return BREAKEND_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + bndtype + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
               locationPartOfString + (forUpstreamLoc ? "1" : "2");
    }

    private static String getRefBaseString(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc,
                                           final ReferenceMultiSource reference) {
        try {
            byte[] refBases = reference.getReferenceBases(forUpstreamLoc ? narl.getLeftJustifiedLeftRefLoc() :
                    narl.getLeftJustifiedRightRefLoc())
                    .getBases();
            return new String(refBases);
        } catch (final IOException ioex) {
            throw new GATKException("Could not read reference for extracting reference bases.", ioex);
        }
    }

    public enum SupportedType {
        INTRA_CHR_STRAND_SWITCH_55,// intra-chromosome strand-switch novel adjacency, alignments left-flanking the novel adjacency
        INTRA_CHR_STRAND_SWITCH_33,// intra-chromosome strand-switch novel adjacency, alignments right-flanking the novel adjacency

        INTRA_CHR_REF_ORDER_SWAP,// intra-chromosome reference-order swap, but NO strand-switch, novel adjacency

        INTER_CHR_STRAND_SWITCH_55,// pair WY in Fig.1 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_STRAND_SWITCH_33,// pair XZ in Fig.1 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER, // the green pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
        INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER; // the red pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
    }

    /**
     * Breakend variant type for inversion suspects: those with novel adjacency between two reference locations
     * on the same chromosome but the novel adjacency brings them together in a strand-switch fashion.
     * This is to be distinguished from the more general "translocation" breakends, which are novel adjacency between
     * reference locations without strand switch if the reference bases are from the same chromosome.
     *
     * Note that dispersed duplication with some copies inverted could also lead to breakpoints with strand switch.
     */
    abstract private static class IntraChromosomalStrandSwitchBreakEnd extends BreakEndVariantType {
        static final Map<String, String> INV55_FLAG = Collections.singletonMap(INV55, "");
        static final Map<String, String> INV33_FLAG = Collections.singletonMap(INV33, "");

        private IntraChromosomalStrandSwitchBreakEnd(final String id, final Allele altAllele, final boolean isTheUpstreamMate,
                                                     final Map<String, String> typeSpecificExtraAttributes) {
            super(id, altAllele, isTheUpstreamMate, typeSpecificExtraAttributes);
        }

        @VisibleForTesting
        static String extractInsertedSequence(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
            final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
            if (ins.isEmpty()) {
                return ins;
            } else {
                return forUpstreamLoc ? ins : SequenceUtil.reverseComplement(ins);
            }
        }
    }

    public static final class IntraChromosomalStrandSwitch55BreakEnd extends IntraChromosomalStrandSwitchBreakEnd {

        @VisibleForTesting
        public IntraChromosomalStrandSwitch55BreakEnd(final String id, final Allele altAllele, final boolean isTheUpstreamMate) {
            super(id, altAllele, isTheUpstreamMate, INV55_FLAG);
        }

        private IntraChromosomalStrandSwitch55BreakEnd(final NovelAdjacencyAndAltHaplotype narl,
                                                       final ReferenceMultiSource reference,
                                                       final boolean isTheUpstreamMate) {
            super(BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                                       extractInsertedSequence(narl, isTheUpstreamMate),
                                        isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc()),
                    isTheUpstreamMate, INV55_FLAG);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference) {
            return new Tuple2<>(new IntraChromosomalStrandSwitch55BreakEnd(narl, reference, true),
                                new IntraChromosomalStrandSwitch55BreakEnd(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc) {
            return Allele.create(refBase + insertedSequence + "]" + novelAdjRefLoc.toString() + "]");
        }
    }

    public static final class IntraChromosomalStrandSwitch33BreakEnd extends IntraChromosomalStrandSwitchBreakEnd {

        @VisibleForTesting
        public IntraChromosomalStrandSwitch33BreakEnd(final String id, final Allele altAllele, final boolean isTheUpstreamMate) {
            super(id, altAllele, isTheUpstreamMate, INV33_FLAG);
        }

        private IntraChromosomalStrandSwitch33BreakEnd(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSource reference,
                                                       final boolean isTheUpstreamMate) {
            super(BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                            extractInsertedSequence(narl, isTheUpstreamMate),
                            isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc()),
                    isTheUpstreamMate, INV33_FLAG);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference) {
            return new Tuple2<>(new IntraChromosomalStrandSwitch33BreakEnd(narl, reference, true),
                                new IntraChromosomalStrandSwitch33BreakEnd(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc) {
            return Allele.create("[" + novelAdjRefLoc.toString() + "[" + insertedSequence + refBase);
        }

    }

    public static final class IntraChromosomeRefOrderSwap extends BreakEndVariantType {

        @VisibleForTesting
        public IntraChromosomeRefOrderSwap(final String id, final Allele altAllele, final boolean isTheUpstreamMate) {
            super(id, altAllele, isTheUpstreamMate, noExtraAttributes);
        }

        private IntraChromosomeRefOrderSwap(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSource reference,
                                            final boolean isTheUpstreamMate) {
            super(BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    constructAltAllele(BreakEndVariantType.getRefBaseString(narl, isTheUpstreamMate, reference),
                            narl.getComplication().getInsertedSequenceForwardStrandRep(),
                            isTheUpstreamMate ? narl.getLeftJustifiedRightRefLoc(): narl.getLeftJustifiedLeftRefLoc(),
                            isTheUpstreamMate),
                    isTheUpstreamMate, noExtraAttributes);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference) {
            return new Tuple2<>(new IntraChromosomeRefOrderSwap(narl, reference, true),
                                new IntraChromosomeRefOrderSwap(narl, reference, false));
        }

        private static Allele constructAltAllele(final String refBase, final String insertedSequence, final SimpleInterval novelAdjRefLoc,
                                                 final boolean forUpstreamLoc) {
            if (forUpstreamLoc) {
                return Allele.create(refBase + insertedSequence + "[" + novelAdjRefLoc.toString() + "[");
            } else {
                return Allele.create("]" + novelAdjRefLoc + "]" + insertedSequence + refBase);
            }
        }
    }

    public static final class InterChromosomeBreakend extends BreakEndVariantType {

        @VisibleForTesting
        public InterChromosomeBreakend(final String id, final Allele altAllele, final boolean isTheUpstreamMate) {
            super(id, altAllele, isTheUpstreamMate, noExtraAttributes);
        }

        private InterChromosomeBreakend(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSource reference,
                                        final boolean isTheUpstreamMate) {
            super(BreakEndVariantType.getIDString(narl, isTheUpstreamMate),
                    constructAltAllele(narl, reference, isTheUpstreamMate),
                    isTheUpstreamMate, noExtraAttributes);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyAndAltHaplotype narl,
                                                                                       final ReferenceMultiSource reference) {

            return new Tuple2<>(new InterChromosomeBreakend(narl, reference, true),
                                new InterChromosomeBreakend(narl, reference, false));
        }

        // see VCF spec 4.2 for BND format ALT allele field for SV, in particular the examples shown in Fig.1, Fig.2 and Fig.5 of Section 5.4
        private static Allele constructAltAllele(final NovelAdjacencyAndAltHaplotype narl, final ReferenceMultiSource reference,
                                                 final boolean forUpstreamLoc) {
            final String refBase = BreakEndVariantType.getRefBaseString(narl, forUpstreamLoc, reference);
            final String insertedSequence = extractInsertedSequence(narl, forUpstreamLoc);
            final SimpleInterval novelAdjRefLoc = forUpstreamLoc ? narl.getLeftJustifiedRightRefLoc() : narl.getLeftJustifiedLeftRefLoc();

            // see Fig.5 of Section 5.4 of spec Version 4.2 (the green pairs)
            final boolean upstreamLocIsFirstInPartner =
                    narl.getTypeInferredFromSimpleChimera().equals(TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER);
            if (narl.getStrandSwitch().equals(StrandSwitch.NO_SWITCH)) {
                if (forUpstreamLoc == upstreamLocIsFirstInPartner) {
                    return Allele.create(refBase + insertedSequence + "[" + novelAdjRefLoc.toString() + "[");
                } else {
                    return Allele.create("]" + novelAdjRefLoc + "]" + insertedSequence + refBase);
                }
            } else if (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE)){
                return Allele.create(refBase + insertedSequence + "]" + novelAdjRefLoc.toString() + "]");
            } else {
                return Allele.create("[" + novelAdjRefLoc.toString() + "[" + insertedSequence + refBase);
            }
        }

        private static String extractInsertedSequence(final NovelAdjacencyAndAltHaplotype narl, final boolean forUpstreamLoc) {
            final String ins = narl.getComplication().getInsertedSequenceForwardStrandRep();
            if (ins.isEmpty()) {
                return ins;
            } else {
                if (narl.getStrandSwitch() == StrandSwitch.NO_SWITCH) {
                    return narl.getComplication().getInsertedSequenceForwardStrandRep();
                } else {
                    return forUpstreamLoc == (narl.getStrandSwitch().equals(StrandSwitch.FORWARD_TO_REVERSE) ) ? ins: SequenceUtil.reverseComplement(ins);
                }
            }
        }
    }
}
