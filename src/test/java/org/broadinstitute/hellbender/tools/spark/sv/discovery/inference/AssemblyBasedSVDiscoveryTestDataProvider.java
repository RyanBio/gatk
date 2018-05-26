package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public abstract class AssemblyBasedSVDiscoveryTestDataProvider {

    public abstract static class AssemblyBasedSVDiscoveryTestDataForSimpleChimera {
        // these are used for building
        public final AlignmentInterval firstAlignment;
        public final AlignmentInterval secondAlignment;
        public final String evidenceAssemblyContigName;
        public final byte[] evidenceContigSeq;

        // these are used for testing
        public final SimpleChimera manuallyCuratedSimpleChimera;
        public final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble;

        public final List<SvType> manuallyCuratedSVTypes;
        public final List<VariantContext> manuallyCuratedVariants;

        public final Class<? extends BreakpointsInference> inferencer;

        public AssemblyBasedSVDiscoveryTestDataForSimpleChimera(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                                                final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                                                final SimpleChimera manuallyCuratedSimpleChimera,
                                                                final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble,
                                                                final List<SvType> manuallyCuratedSVTypes,
                                                                final List<VariantContext> manuallyCuratedVariants,
                                                                final Class<? extends BreakpointsInference> inferencer) {
            this.firstAlignment = firstAlignment;
            this.secondAlignment = secondAlignment;
            this.evidenceAssemblyContigName = evidenceAssemblyContigName;
            this.evidenceContigSeq = evidenceContigSeq;
            this.manuallyCuratedSimpleChimera = manuallyCuratedSimpleChimera;
            this.manuallyCuratedBiPathBubble = manuallyCuratedBiPathBubble;
            this.manuallyCuratedSVTypes = manuallyCuratedSVTypes;
            this.manuallyCuratedVariants = manuallyCuratedVariants;
            this.inferencer = inferencer;
        }

        abstract public SAMSequenceDictionary getAppropriateDictionary();

        abstract public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer();
    }

    // data block ======================================================================================================

    private static String makeID(final String prefix, final String chr1, final int start, final String chr2, final int stop,
                                 final String postfix) {
        return prefix + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                SvType.makeLocationPartOfID(chr1, start, chr2, stop) + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + postfix;
    }

    static final byte[] EMPTY_BYTE_ARRAY = new byte[]{};

    private static final Allele INV_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INV.name()), false);
    private static final Allele DEL_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.DEL.name()), false);
    private static final Allele INS_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INS.name()), false);
    private static final Allele DUP_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.DUP.name()), false);

    // utils block =====================================================================================================

    /**
     * Note that {@code delRange} is expected to be pre-process to VCF spec compatible,
     * e.g. if chr1:101-200 is deleted, then {@code delRange} should be chr1:100-200
     */
    static final VariantContextBuilder makeDeletion(final SimpleInterval delRange, final Allele refAllele) {

        return new VariantContextBuilder()
                .chr(delRange.getContig()).start(delRange.getStart()).stop(delRange.getEnd())
                .alleles(Arrays.asList(refAllele, DEL_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.DEL.name(), delRange.getContig(), delRange.getStart(), delRange.getContig(), delRange.getEnd(), ""))
                .attribute(VCFConstants.END_KEY, delRange.getEnd())
                .attribute(SVLEN, - delRange.size() + 1)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.name());
    }

    static final SvType makeDeletionType(final SimpleInterval delRange, final boolean isFromDupContraction) {
        return new SimpleSVType.Deletion(
                makeID((isFromDupContraction? DUP_TAN_CONTRACTION_STRING : SimpleSVType.SupportedType.DEL.name()),
                        delRange.getContig(), delRange.getStart(), delRange.getContig(), delRange.getEnd(), ""),
                DEL_SYMB_ALLELE, -delRange.size()+ 1,
                isFromDupContraction ? Collections.singletonMap(DUP_TAN_CONTRACTION_STRING, "") :Collections.emptyMap());
    }

    static final VariantContextBuilder makeInversion(final SimpleInterval invertedRegion, final Allele refAllele) {
        return new VariantContextBuilder()
                .chr(invertedRegion.getContig()).start(invertedRegion.getStart() - 1).stop(invertedRegion.getEnd())     // TODO: 5/2/18 VCF spec doesn't requst left shift by 1 for inversion POS
                .alleles(Arrays.asList(refAllele, INV_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.INV.name(), invertedRegion.getContig(), invertedRegion.getStart() - 1, invertedRegion.getContig(), invertedRegion.getEnd(), ""))
                .attribute(VCFConstants.END_KEY, invertedRegion.getEnd())
                .attribute(SVLEN, 0)                                                                 // TODO: 5/2/18 this is following VCF spec,
                .attribute(SVTYPE, SimpleSVType.SupportedType.INV.name());
    }

    static final SvType makeInversionType(final SimpleInterval invRange, final boolean isInv55) {
        return new SimpleSVType.Inversion(
                makeID((isInv55? INV55 : INV33),
                        invRange.getContig(), invRange.getStart() - 1, invRange.getContig(), invRange.getEnd(), ""),
                INV_SYMB_ALLELE, invRange.size(),
                Collections.singletonMap((isInv55) ? INV55 : INV33, ""));
    }

    static final VariantContextBuilder makeInsertion(final String chr, final int pos, final int end, final int svLen,
                                                     final Allele refAllele) {

        return new VariantContextBuilder().chr(chr).start(pos).stop(end)
                .alleles(Arrays.asList(refAllele, INS_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.INS.name(), chr, pos, chr, end, ""))
                .attribute(VCFConstants.END_KEY, end)
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, SimpleSVType.SupportedType.INS.name());
    }

    static final SvType makeInsertionType(final SimpleInterval insertionPos, final int insLen) {
        return new SimpleSVType.Insertion(
                makeID(SimpleSVType.SupportedType.INS.name(),
                        insertionPos.getContig(), insertionPos.getStart(), insertionPos.getContig(), insertionPos.getEnd(), ""),
                INS_SYMB_ALLELE, insLen,
                Collections.emptyMap());
    }

    static final VariantContextBuilder makeTandemDuplication(final SimpleInterval duplicatedRange, final Allele refAllele,
                                                             final int svLen, final Map<String, Object> attributes) {
        return new VariantContextBuilder().chr(duplicatedRange.getContig()).start(duplicatedRange.getStart() - 1).stop(duplicatedRange.getStart() - 1) // TODO: 5/24/18 by left shifting 1, we are treating it as insertions
                .id(makeID(DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING,
                        duplicatedRange.getContig(), duplicatedRange.getStart() - 1 , duplicatedRange.getContig(), duplicatedRange.getStart() - 1, ""))
                .alleles(Arrays.asList(refAllele, DUP_SYMB_ALLELE))
                .attribute(VCFConstants.END_KEY, duplicatedRange.getStart() - 1) // TODO: 5/24/18 see todo above
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DUP.name());
    }

    static final SvType makeTandemDuplicationType(final SimpleInterval duplicatedRange, final int svLen) {
        return new SimpleSVType.DuplicationTandem(
                makeID(DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING,
                        duplicatedRange.getContig(), duplicatedRange.getStart(), duplicatedRange.getContig(), duplicatedRange.getEnd(), ""),
                DUP_SYMB_ALLELE, svLen,
                Collections.singletonMap(DUP_TAN_EXPANSION_STRING, ""));
    }

    static final VariantContextBuilder makeBND(final SimpleInterval upstreamLoc, final SimpleInterval dnstreamLoc,
                                               final Allele refAllele, final String insertedSeq, final String bndSubtypeString,
                                               final boolean forUpstreamLoc, final boolean refBaseFirst, final boolean bracketPointsLeft) {

        final Allele altAllele;
        if (refBaseFirst) {
            if (bracketPointsLeft)
                altAllele = Allele.create(refAllele.getBaseString() + insertedSeq + "]" + (forUpstreamLoc ? dnstreamLoc.toString() : upstreamLoc.toString()) + "]");
            else
                altAllele = Allele.create(refAllele.getBaseString() + insertedSeq + "[" + (forUpstreamLoc ? dnstreamLoc.toString() : upstreamLoc.toString()) + "[");
        } else {
            if (bracketPointsLeft)
                altAllele = Allele.create("]" + (forUpstreamLoc ? dnstreamLoc.toString() : upstreamLoc.toString()) + "]" + insertedSeq + refAllele.getBaseString());
            else
                altAllele = Allele.create("[" + (forUpstreamLoc ? dnstreamLoc.toString() : upstreamLoc.toString()) + "[" + insertedSeq + refAllele.getBaseString());
        }
        if (forUpstreamLoc) {
            return new VariantContextBuilder().chr(upstreamLoc.getContig()).start(upstreamLoc.getStart()).stop(upstreamLoc.getEnd())
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .id(makeID(BREAKEND_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + bndSubtypeString,
                            upstreamLoc.getContig(), upstreamLoc.getStart(), dnstreamLoc.getContig(), dnstreamLoc.getEnd(), "_1"))
                    .attribute(SVTYPE, BREAKEND_STR);
        } else {
            return new VariantContextBuilder().chr(upstreamLoc.getContig()).start(upstreamLoc.getStart()).stop(upstreamLoc.getEnd())
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .id(makeID(BREAKEND_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + bndSubtypeString,
                            upstreamLoc.getContig(), upstreamLoc.getStart(), dnstreamLoc.getContig(), dnstreamLoc.getEnd(), "_2"))
                    .attribute(SVTYPE, BREAKEND_STR);
        }
    }

    static final SvType makeBNDType(final String id, final Allele altAllele, final boolean isTheUpstreamMate,
                                    final BreakEndVariantType.SupportedType type) {
        switch (type) {
            case INTRA_CHR_STRAND_SWITCH_55:
                return new BreakEndVariantType.IntraChromosomalStrandSwitch55BreakEnd(id, altAllele, isTheUpstreamMate);
            case INTRA_CHR_STRAND_SWITCH_33:
                return new BreakEndVariantType.IntraChromosomalStrandSwitch33BreakEnd(id, altAllele, isTheUpstreamMate);
            case INTRA_CHR_REF_ORDER_SWAP:
                return new BreakEndVariantType.IntraChromosomeRefOrderSwap(id, altAllele, isTheUpstreamMate);
            case INTER_CHR_STRAND_SWITCH_55:
            case INTER_CHR_STRAND_SWITCH_33:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER:
                return new BreakEndVariantType.InterChromosomeBreakend(id, altAllele, isTheUpstreamMate);
            default:
                throw new GATKException("Unrecognized type: " + type.name());
        }
    }
}
