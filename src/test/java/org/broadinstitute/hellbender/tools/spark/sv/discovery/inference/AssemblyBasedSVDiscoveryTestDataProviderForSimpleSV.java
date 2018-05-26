package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch.NO_SWITCH;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.CONTIG_NAMES;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.MAX_ALIGN_LENGTH;

/**
 * Provides test data for testing several methods involved in the SV variant caller,
 * but specifically on simple types.
 * NO TESTS ARE RUN IN THIS PARTICULAR CLASS.
 */
public final class AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV {

    public static final class TestDataForSimpleSV extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera {
        public final DistancesBetweenAlignmentsOnRefAndOnRead distances;


        private TestDataForSimpleSV(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                    final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                    final SimpleChimera manuallyCuratedSimpleChimera,
                                    final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble,
                                    final List<SvType> manuallyCuratedSVTypes,
                                    final List<VariantContext> manuallyCuratedVariants,
                                    final DistancesBetweenAlignmentsOnRefAndOnRead distances,
                                    final Class<? extends BreakpointsInference> inferencer) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, manuallyCuratedSimpleChimera, manuallyCuratedBiPathBubble, manuallyCuratedSVTypes, manuallyCuratedVariants, inferencer);
            this.distances = distances;
        }

        @Override
        public SAMSequenceDictionary getAppropriateDictionary() {
            return TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return inferencer;
        }
    }

    public static final boolean testDataInitialized;
    public static final TestDataForSimpleSV forSimpleDeletion_plus;
    public static final TestDataForSimpleSV forSimpleDeletion_minus;
    public static final TestDataForSimpleSV forSimpleInsertion_plus;
    public static final TestDataForSimpleSV forSimpleInsertion_minus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_fudgedDel_plus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_fudgedDel_minus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_fatIns_plus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_fatIns_minus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_DelAndIns_plus;
    public static final TestDataForSimpleSV forLongRangeSubstitution_DelAndIns_minus;
    public static final TestDataForSimpleSV forDeletionWithHomology_plus;
    public static final TestDataForSimpleSV forDeletionWithHomology_minus;
    public static final TestDataForSimpleSV forSimpleTanDupContraction_plus;
    public static final TestDataForSimpleSV forSimpleTanDupContraction_minus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansion_ins_plus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansion_ins_minus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansion_dup_plus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansion_dup_minus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_ins_plus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_ins_minus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_dup_plus;
    public static final TestDataForSimpleSV forSimpleTanDupExpansionWithNovelIns_dup_minus;
    public static final TestDataForSimpleSV forComplexTanDup_1to2_pseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_1to2_pseudoHom_minus;
    public static final TestDataForSimpleSV forComplexTanDup_2to1_pseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_2to1_pseudoHom_minus;
    public static final TestDataForSimpleSV forComplexTanDup_3to2_noPseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_3to2_noPseudoHom_minus;
    public static final TestDataForSimpleSV forComplexTanDup_2to3_noPseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_2to3_noPseudoHom_minus;
    public static final TestDataForSimpleSV forComplexTanDup_1to2_short_pseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_1to2_short_pseudoHom_minus;
    public static final TestDataForSimpleSV forComplexTanDup_2to3_short_noPseudoHom_plus;
    public static final TestDataForSimpleSV forComplexTanDup_2to3_short_noPseudoHom_minus;

    static {
        try ( final ByteArrayOutputStream outputStream = new ByteArrayOutputStream() ){

            final List<TestDataForSimpleSV> simpleDeletion = forSimpleDeletion(outputStream);
            forSimpleDeletion_plus = simpleDeletion.get(0);
            forSimpleDeletion_minus = simpleDeletion.get(1);

            final List<TestDataForSimpleSV> simpleInsertion = forSimpleInsertion(outputStream);
            forSimpleInsertion_plus = simpleInsertion.get(0);
            forSimpleInsertion_minus = simpleInsertion.get(1);

            final List<TestDataForSimpleSV> longRangeSubstitution = forLongRangeSubstitution(outputStream);
            forLongRangeSubstitution_fudgedDel_plus = longRangeSubstitution.get(0);
            forLongRangeSubstitution_fudgedDel_minus = longRangeSubstitution.get(1);
            forLongRangeSubstitution_fatIns_plus = longRangeSubstitution.get(2);
            forLongRangeSubstitution_fatIns_minus = longRangeSubstitution.get(3);
            forLongRangeSubstitution_DelAndIns_plus = longRangeSubstitution.get(4);
            forLongRangeSubstitution_DelAndIns_minus = longRangeSubstitution.get(5);

            final List<TestDataForSimpleSV> deletionWithHomology = forDeletionWithHomology(outputStream);
            forDeletionWithHomology_plus = deletionWithHomology.get(0);
            forDeletionWithHomology_minus = deletionWithHomology.get(1);

            final List<TestDataForSimpleSV> simpleTandemDuplicationContraction = forSimpleTandemDuplicationContraction(outputStream);
            forSimpleTanDupContraction_plus = simpleTandemDuplicationContraction.get(0);
            forSimpleTanDupContraction_minus = simpleTandemDuplicationContraction.get(1);

            final List<TestDataForSimpleSV> simpleTandemDuplicationExpansion = forSimpleTandemDuplicationExpansion(outputStream);
            forSimpleTanDupExpansion_ins_plus = simpleTandemDuplicationExpansion.get(0);
            forSimpleTanDupExpansion_ins_minus = simpleTandemDuplicationExpansion.get(1);
            forSimpleTanDupExpansion_dup_plus = simpleTandemDuplicationExpansion.get(2);
            forSimpleTanDupExpansion_dup_minus = simpleTandemDuplicationExpansion.get(3);

            final List<TestDataForSimpleSV> simpleTandemDuplicationExpansionWithNovelInsertion = forSimpleTandemDuplicationExpansionWithNovelInsertion(outputStream);
            forSimpleTanDupExpansionWithNovelIns_ins_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(0);
            forSimpleTanDupExpansionWithNovelIns_ins_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(1);
            forSimpleTanDupExpansionWithNovelIns_dup_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(2);
            forSimpleTanDupExpansionWithNovelIns_dup_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(3);

            final List<TestDataForSimpleSV> complexTandemDuplication = forComplexTandemDuplication(outputStream);
            forComplexTanDup_1to2_pseudoHom_plus = complexTandemDuplication.get(0);
            forComplexTanDup_1to2_pseudoHom_minus = complexTandemDuplication.get(1);
            forComplexTanDup_2to1_pseudoHom_plus = complexTandemDuplication.get(2);
            forComplexTanDup_2to1_pseudoHom_minus = complexTandemDuplication.get(3);
            forComplexTanDup_3to2_noPseudoHom_plus = complexTandemDuplication.get(4);
            forComplexTanDup_3to2_noPseudoHom_minus = complexTandemDuplication.get(5);
            forComplexTanDup_2to3_noPseudoHom_plus = complexTandemDuplication.get(6);
            forComplexTanDup_2to3_noPseudoHom_minus = complexTandemDuplication.get(7);

            final List<TestDataForSimpleSV> shortComplexTandemDuplication = forComplexTandemDuplicationIns(outputStream);
            forComplexTanDup_1to2_short_pseudoHom_plus = shortComplexTandemDuplication.get(0);
            forComplexTanDup_1to2_short_pseudoHom_minus = shortComplexTandemDuplication.get(1);
            forComplexTanDup_2to3_short_noPseudoHom_plus = shortComplexTandemDuplication.get(2);
            forComplexTanDup_2to3_short_noPseudoHom_minus = shortComplexTandemDuplication.get(3);

            testDataInitialized = true;
        } catch (final Exception ioex) {
            throw new GATKException("Failed to create test data ", ioex);
        }
    }

    public static List<TestDataForSimpleSV> getAllTestData() {
        final List<TestDataForSimpleSV> testDataForSimpleSVs = Arrays.asList(
                forSimpleDeletion_plus,
                forSimpleDeletion_minus,
                forSimpleInsertion_plus,
                forSimpleInsertion_minus,
                forLongRangeSubstitution_fudgedDel_plus,
                forLongRangeSubstitution_fudgedDel_minus,
                forLongRangeSubstitution_fatIns_plus,
                forLongRangeSubstitution_fatIns_minus,
                forLongRangeSubstitution_DelAndIns_plus,
                forLongRangeSubstitution_DelAndIns_minus,
                forDeletionWithHomology_plus,
                forDeletionWithHomology_minus,
                forSimpleTanDupContraction_plus,
                forSimpleTanDupContraction_minus,
                forSimpleTanDupExpansion_ins_plus,
                forSimpleTanDupExpansion_ins_minus,
                forSimpleTanDupExpansion_dup_plus,
                forSimpleTanDupExpansion_dup_minus,
                forSimpleTanDupExpansionWithNovelIns_ins_plus,
                forSimpleTanDupExpansionWithNovelIns_ins_minus,
                forSimpleTanDupExpansionWithNovelIns_dup_plus,
                forSimpleTanDupExpansionWithNovelIns_dup_minus,
                forComplexTanDup_1to2_pseudoHom_plus,
                forComplexTanDup_1to2_pseudoHom_minus,
                forComplexTanDup_2to1_pseudoHom_plus,
                forComplexTanDup_2to1_pseudoHom_minus,
                forComplexTanDup_3to2_noPseudoHom_plus,
                forComplexTanDup_3to2_noPseudoHom_minus,
                forComplexTanDup_2to3_noPseudoHom_plus,
                forComplexTanDup_2to3_noPseudoHom_minus,
                forComplexTanDup_1to2_short_pseudoHom_plus,
                forComplexTanDup_1to2_short_pseudoHom_minus,
                forComplexTanDup_2to3_short_noPseudoHom_plus,
                forComplexTanDup_2to3_short_noPseudoHom_minus);
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    // same event, two representations from opposite strands
    public static List<Tuple2<TestDataForSimpleSV, TestDataForSimpleSV>> getAllTestDataPaired() {
        final List<TestDataForSimpleSV> allTestData = getAllTestData();
        final List<Tuple2<TestDataForSimpleSV, TestDataForSimpleSV>> testDataForSimpleSVs
                = new ArrayList<>(allTestData.size()/2);
        for (int i = 0; i < allTestData.size() - 1; i += 2) {
            testDataForSimpleSVs.add(new Tuple2<>(allTestData.get(i), allTestData.get(i+1)));
        }
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    private static VariantContextBuilder addStandardAttributes(final VariantContextBuilder inputBuilder,
                                                               final int alignLength, final String contigName,
                                                               final String variantType, final int end,
                                                               final int svLen, final String altSeq) {
        return inputBuilder
                .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1)
                .attribute(ALIGN_LENGTHS, alignLength).attribute(MAX_ALIGN_LENGTH, alignLength)
                .attribute(CONTIG_NAMES, contigName)
                .attribute(VCFConstants.END_KEY, end)
                .attribute(SVTYPE, variantType)
                .attribute(SVLEN, svLen)
                .attribute(SEQ_ALT_HAPLOTYPE, altSeq);
    }

    /**
     * 40-'A' + 10-'C'+10-'T' + 40-'G' where the segment 10-'C'+10-'T' is deleted (forward strand representation description).
     *
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleDeletion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        // simple deletion '+' strand representation
        final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'G');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        String contigName = "simple_del_+";

        final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000060-17000060");
        final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", "");
        final byte[] altSeq = EMPTY_BYTE_ARRAY;
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SIMPLE_DEL, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100040), 1 ,40, TextCigarCodec.decode("40M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100061, 100100), 41 ,80, TextCigarCodec.decode("40S40M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 0, 17000040, 17000061, 40, 41);
        final List<SvType> type = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000060"), false));
        final List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000060"), Allele.create("A", true)),
                        40, contigName, SimpleSVType.SupportedType.DEL.name(), 17000060, -20, "").make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple deletion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        contigName = "simple_del_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 1 ,40, TextCigarCodec.decode("40M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 41 ,80, TextCigarCodec.decode("40S40M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 0, 17000040, 17000061, 40, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * 100-'A' + 100-'T' and a 50 bases of 'C' is inserted at the A->T junction point (forward strand description)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleInsertion(final ByteArrayOutputStream outputStream) throws IOException {
        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple insertion '+' strand representation
        final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'A');
        final byte[] insertedSeq  = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(50, (byte)'C');
        final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'T');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(insertedSeq);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        String contigName = "simple_ins_+";
        final String insSeqString = new String(insertedSeq);

        final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000100-17000100");
        final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000100-17000100");
        final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
        final byte[] altSeq = Arrays.copyOf(insertedSeq, insertedSeq.length);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SIMPLE_INS, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000100), 1 ,100, TextCigarCodec.decode("100M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000101, 17000200), 151 ,250, TextCigarCodec.decode("100S100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, 50, 17000100, 17000101, 100, 151);
        final List<SvType> type = Collections.singletonList(makeInsertionType(new SimpleInterval("21:17000100-17000100"), 50));
        final List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeInsertion("21", 17000100, 17000100, 50, Allele.create("G", true)),
                        100, contigName, SimpleSVType.SupportedType.INS.name(), 17000100, 50, insSeqString).make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple insertion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(insertedSeq);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(insertedSeq);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        contigName = "simple_ins_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000101, 17000200), 1 ,100, TextCigarCodec.decode("100M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000100), 151 ,250, TextCigarCodec.decode("100S100M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, 50, 17000100, 17000101, 100, 151);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * fudged deletion case:
     * 100-'A' + 100-'G' where the middle 30-'A'+30-'G' is substituted with 10-'C' (forward strand representation)
     * fat insertion case:
     * 50-'A' + 50-'G' where the middle 10-'A'+10-'G' is substituted with 60-'C' (forward strand representation)
     * Two linked variants case:
     * 100-'A' + 100-'G' where the middle 30-'A'+30-'G' is substituted with 55-'C' (forward strand representation)
     */
    private static List<TestDataForSimpleSV>
    forLongRangeSubstitution(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {//fudged deletion case
            // '+' strand representation
            final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'A');
            final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'G');
            final byte[] substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(10, (byte)'C');
            outputStream.reset();
            outputStream.write(leftRefFlank, 0, 70);outputStream.write(substitution);outputStream.write(rightRefFlank, 0, 70);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "rpl_fudged_del_+";
            final String insSeqString = StringUtil.bytesToString(substitution);

            final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000070-17000070");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000130-17000130");
            final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] altSeq = Arrays.copyOf(substitution, substitution.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, RPL, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 1 ,70, TextCigarCodec.decode("70M80S"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 81 ,150, TextCigarCodec.decode("80S70M"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 10, 17000070, 17000131, 70, 81);
            final List<SvType> type = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000070-17000130"), false));
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeDeletion(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true)),
                    70, contigName, SimpleSVType.SupportedType.DEL.name(), 17000130, -60, insSeqString)
                    .attribute(INSERTED_SEQUENCE, insSeqString).attribute(INSERTED_SEQUENCE_LENGTH, insSeqString.length()).make());
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(substitution);
            outputStream.reset();
            outputStream.write(rightRefFlank, 0, 70);outputStream.write(substitution);outputStream.write(leftRefFlank, 0, 70);
            contigSeq = outputStream.toByteArray();
            contigName = "rpl_fudged_del_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 1 ,70, TextCigarCodec.decode("70M80S"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 81 ,150, TextCigarCodec.decode("80S70M"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 10, 17000070, 17000131, 70, 81);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        {//fat insertion case
            // '+' strand representation
            final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(50, (byte)'A');
            final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(50, (byte)'G');
            final byte[] substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(60, (byte)'C');
            outputStream.reset();
            outputStream.write(leftRefFlank, 0, 40);outputStream.write(substitution);outputStream.write(rightRefFlank, 0, 40);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "rpl_fat_ins_+";
            final String insSeqString = StringUtil.bytesToString(substitution);

            final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000040-17000040");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000060-17000060");
            final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] altSeq = Arrays.copyOf(substitution, substitution.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, RPL, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 1 ,40, TextCigarCodec.decode("40M100S"), true, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 101 ,140, TextCigarCodec.decode("100S40M"), true, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 60, 17000040, 17000061, 40, 101);
            final List<SvType> type = Collections.singletonList(makeInsertionType(new SimpleInterval("21:17000040-17000040"), 60));
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeInsertion("21", 17000040, 17000060, 60, Allele.create("G", true)),
                            40, contigName, SimpleSVType.SupportedType.INS.name(), 17000060, 60, insSeqString).make());
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(substitution);
            outputStream.reset();
            outputStream.write(rightRefFlank, 0, 40);outputStream.write(substitution);outputStream.write(leftRefFlank, 0, 40);
            contigSeq = outputStream.toByteArray();
            contigName = "rpl_fat_ins_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000061, 17000100), 1 ,40, TextCigarCodec.decode("40M100S"), false, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000040), 101 ,140, TextCigarCodec.decode("100S40M"), false, 60, 0, 60, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(20, 60, 17000040, 17000061, 40, 101);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        {//two linked variants case
            // '+' strand representation
            final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'A');
            final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'G');
            final byte[] substitution = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(55, (byte)'C');
            outputStream.reset();
            outputStream.write(leftRefFlank, 0, 70);outputStream.write(substitution);outputStream.write(rightRefFlank, 0 ,70);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "rpl_linked_del_ins_+";
            final String insSeqString = StringUtil.bytesToString(substitution);

            final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000070-17000070");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000130-17000130");
            final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications("", insSeqString);
            final byte[] altSeq = Arrays.copyOf(substitution, substitution.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, RPL, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 1 ,70, TextCigarCodec.decode("70M125S"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 126 ,195, TextCigarCodec.decode("125S70M"), true, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 55, 17000070, 17000131, 70, 126);
            final List<SvType> type = Arrays.asList(
                    makeDeletionType(new SimpleInterval("21:17000070-17000130"), false),
                    makeInsertionType(new SimpleInterval("21:17000070-17000070"), 55)
            );
            final List<VariantContext> variant = Arrays.asList(
                    addStandardAttributes(makeDeletion(new SimpleInterval("21:17000070-17000130"), Allele.create("T", true)), 70, contigName, SimpleSVType.SupportedType.DEL.name(), 17000130, -60, "")
                            .attribute(LINK, SimpleSVType.SupportedType.INS.name() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + SvType.makeLocationPartOfID("21", 17000070, "21", 17000070)).make(),
                    addStandardAttributes(makeInsertion("21", 17000070, 17000070, 55, Allele.create("T", true)), 70, contigName, SimpleSVType.SupportedType.INS.name(), 17000070, 55, "")
                            .attribute(LINK, SimpleSVType.SupportedType.DEL.name() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + SvType.makeLocationPartOfID("21", 17000070, "21", 17000130)).make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

            // '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(substitution);
            outputStream.reset();
            outputStream.write(rightRefFlank, 0, 70);outputStream.write(substitution);outputStream.write(leftRefFlank, 0 ,70);
            contigSeq = outputStream.toByteArray();
            contigName = "rpl_linked_del_ins_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 17000131, 17000200), 1 ,70, TextCigarCodec.decode("70M125S"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000070), 126 ,195, TextCigarCodec.decode("125S70M"), false, 60, 0, 70, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(60, 55, 17000070, 17000131, 70, 126);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));
        }

        return result;
    }

    /**
     * 40-'C' + 'ATCG' + 34 bases of unique sequence + 'ATCG' + 40-'T' is shrunk to 40-'C' + 'ATCG' + 40-'T' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forDeletionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple deletion with homology '+' strand representation
        final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'C');
        final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'T');
        final byte[] homology = new byte[]{'A', 'T', 'C', 'G'};
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(homology);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        String contigName = "simple_del_with_hom_+";

        final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000078-17000078");
        final BreakpointComplications complications = new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications(new String(homology), "");
        final byte[] altSeq = EMPTY_BYTE_ARRAY;
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SIMPLE_DEL, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000044), 1 ,44, TextCigarCodec.decode("44M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000079, 17000122), 41 ,84, TextCigarCodec.decode("40S44M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(34, -4, 17000044, 17000044, 44, 41);
        final List<SvType> type = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000078"), false));
        final List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000078"), Allele.create("A", true)),
                        40, contigName, SimpleSVType.SupportedType.DEL.name(), 17000078, -38, "").make());
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        // simple deletion with homology '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(homology);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(homology);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        contigName = "simple_del_with_hom_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000079, 17000122), 1 ,44, TextCigarCodec.decode("44M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000044), 41 ,84, TextCigarCodec.decode("40S44M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(34, -4, 17000044, 17000044, 44, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances,  BreakpointsInference.SimpleInsertionDeletionBreakpointsInference.class));

        return result;
    }

    /**
     * 40-'A' + 20-'C' + 40-'G' is shrunk to 40-'A' + 10-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationContraction(final ByteArrayOutputStream outputStream) {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        // simple tandem duplication contraction '+' strand representation
        final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'G');
        final byte[] doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(20, (byte)'C');
        outputStream.reset();
        outputStream.write(leftRefFlank, 0, 40);outputStream.write(doubleDup, 0, 10);outputStream.write(rightRefFlank, 0 ,40);
        byte[] contigSeq = outputStream.toByteArray();
        String contigName = "simple_del_dup_contraction_+";

        final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000040-17000040");
        final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000050-17000050");
        final BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("CCCCCCCCCC", "",
                new SimpleInterval("21:17000041-17000050"), 2, 1, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Collections.singletonList(Strand.POSITIVE), Collections.emptyList());
        byte[] altSeq = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(10, (byte)'C');
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, DEL_DUP_CONTRACTION, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 1 ,50, TextCigarCodec.decode("50M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 17000051, 17000100), 41 ,90, TextCigarCodec.decode("40S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, -10, 100050, 100051, 50, 41);
        final List<SvType> type = Collections.singletonList(makeDeletionType(new SimpleInterval("21:17000040-17000050"), true));
        final List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeDeletion(new SimpleInterval("21:17000040-17000050"), Allele.create("A", true)), 50, contigName, SimpleSVType.SupportedType.DEL.name(), 17000050, -10, "")
                        .attribute(DUP_REPEAT_UNIT_REF_SPAN, "21:17000040-17000050").attribute(DUPLICATION_NUMBERS, "2,1").attribute(DUP_TAN_CONTRACTION_STRING, "").make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

        // simple tandem duplication contraction '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(doubleDup);
        outputStream.reset();
        outputStream.write(rightRefFlank, 0, 40);outputStream.write(doubleDup, 0, 10);outputStream.write(leftRefFlank, 0 ,40);
        contigSeq = outputStream.toByteArray();
        contigName = "simple_del_dup_contraction_-";
        region1 = new AlignmentInterval(new SimpleInterval("21", 17000051, 17000100), 1 ,50, TextCigarCodec.decode("50M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 17000001, 17000050), 41 ,90, TextCigarCodec.decode("40S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(0, -10, 17000050, 17000051, 50, 41);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

        return result;
    }

    /**
     * case that will be called as insertion
     * 40-'A' + 10-'C' + 40-'G' is expanded to 40-'A' + 20-'C' + 40-'G' (forward strand representation)
     *
     * case that will be called as duplication
     * 40-'A' + 55-'C' + 40-'G' is expanded to 40-'A' + 110-'C' + 40-'G' (forward strand representation)
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationExpansion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {// insertion case
            // '+' strand representation
            final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'A');
            final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'G');
            final byte[] doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(20, (byte)'C');
            outputStream.reset();
            outputStream.write(leftRefFlank);outputStream.write(doubleDup);outputStream.write(rightRefFlank);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "simple_dup_exp_too_small_+";

            final SimpleInterval leftBreakpoint = new SimpleInterval("21:17000040-17000040");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:17000040-17000040");
            final BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", "",
                    new SimpleInterval("21:100041-100050"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("10M", "10M"));
            final byte[] altSeq = Arrays.copyOf(doubleDup, doubleDup.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_EXPANSION, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 1 ,50, TextCigarCodec.decode("50M50S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100041, 100090), 51 ,100, TextCigarCodec.decode("50S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-10, 0, 100050, 100041, 50, 51);
            final List<SvType> type = Collections.singletonList(makeInsertionType());
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeInsertion(), )
                            .attribute().attribute().attribute().attribute().make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(doubleDup);
            outputStream.reset();
            outputStream.write(rightRefFlank);outputStream.write(doubleDup);outputStream.write(leftRefFlank);
            contigSeq = outputStream.toByteArray();
            contigName = "simple_dup_exp_too_small_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 100041, 100090), 1 ,50, TextCigarCodec.decode("50M50S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 51 ,100, TextCigarCodec.decode("50S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-10, 0, 100050, 100041, 50, 51);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        {// duplication case
            // '+' strand representation
            final byte[] leftRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'A');
            final byte[] rightRefFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(40, (byte)'G');
            final byte[] doubleDup = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(110, (byte)'C');
            outputStream.reset();
            outputStream.write(leftRefFlank);outputStream.write(doubleDup);outputStream.write(rightRefFlank);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "simple_dup_exp_+";

            final SimpleInterval leftBreakpoint = new SimpleInterval("21:100040-100040");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:100040-100040");
            final BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", "",
                    new SimpleInterval("21:100041-100095"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("55M", "55M"));
            final byte[] altSeq = Arrays.copyOf(doubleDup, doubleDup.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_EXPANSION, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100095), 1 ,95, TextCigarCodec.decode("95M95S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100041, 100135), 96 ,190, TextCigarCodec.decode("95S95M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, 0, 100095, 100041, 95, 96);
            final List<SvType> type = Collections.singletonList(makeTandemDuplicationType());
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeTandemDuplication(), )
                            .attribute().attribute().attribute().attribute().make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(doubleDup);
            outputStream.reset();
            outputStream.write(rightRefFlank);outputStream.write(doubleDup);outputStream.write(leftRefFlank);
            contigSeq = outputStream.toByteArray();
            contigName = "simple_dup_exp_-";
            region1 = new AlignmentInterval(new SimpleInterval("21", 100041, 100135), 1 ,95, TextCigarCodec.decode("95M95S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100095), 96 ,190, TextCigarCodec.decode("95S95M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, 0, 100095, 100041, 95, 96);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        return result;
    }

    /**
     * Real event, which will be output as INS (but the event was actually from a hg38 sample, but doesn't matter)
     * repeat:     chr21:26849022-26849037
     * repeat sequence: CCGGGAAATGCTTTTT
     * insertedSequenceForwardStrandRep: TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC
     *
     * Real event, which will be output as DUP
     * leftFlank:  chr21:25297101-25297163
     * repeat:     chr21:25297164-25297252
     * rightFlank: chr21:25297253-25297300
     * GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT
     * AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG
     * CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT
     *
     * insertedSequenceForwardStrandRep: CTCTCTCTCT
     */
    private static List<TestDataForSimpleSV>
    forSimpleTandemDuplicationExpansionWithNovelInsertion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSV> result = new ArrayList<>();

        {// complications slightly different due to small difference in assembly contig sequence
            AlignmentInterval region1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm029081:tig00000\t0\t21\t26847644\t60\t1394M1675S\t*\t0\t0\tTATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26849022,+,1704S657M2I706M,60,2;chr10,97348533,+,1388S317M1364S,0,0;\tMD:Z:1204A189\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1389\tXS:i:0", true);
            AlignmentInterval region2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm029081:tig00000\t2048\t21\t26849022\t60\t1704H657M2I706M\t*\t0\t0\tCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26847644,+,1394M1675S,60,1;chr10,97348533,+,1388S317M1364S,0,0;\tMD:Z:1363\tRG:Z:GATKSVContigAlignments\tNM:i:2\tAS:i:1345\tXS:i:0", true);

            String contigName = "simple_dup_exp_too_small_1_2_with_ins_+";
            byte[] contigSeq = "TATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT".getBytes();
            String insertedSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC";
            final SimpleInterval leftBreakpoint = new SimpleInterval("21:26849021-26849021");
            final SimpleInterval rightBreakpoint = new SimpleInterval("21:26849021-26849021");
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", insertedSeq,
                    new SimpleInterval("21:26849022-26849037"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("16M", "16M"));
            final byte[] altSeq = "CCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTT".getBytes();

            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_EXPANSION, altSeq);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(14, 310, 26849037, 26849022, 1394, 1705);
            final List<SvType> type = Collections.singletonList(makeInsertionType());
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeInsertion(), )
                            .attribute().attribute().attribute().attribute().make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            contigName = "simple_dup_exp_too_small_1_2_with_ins_-";
            contigSeq = "AGAGATAATATGTTAATATTAATTTGCTTTAATTAATAAAATATTCCTAGAATATGATTTTAAATTCACATGACATTCGGAAGTCACTTATATTTGCAGTATTCAAGAATTACATGAAGTATGTGCTCATTAAGAGCGATTTTAAATTTTATCTTTTGTCACCACCTGGTGGTAGAACATTTTCTCACTTTTCTTTTTTCATTTTAAGCCTTACGTTGTTGCTAAGGAAAATTTCTAAACCAGGGAAACCAGTTAACCAGTATTGATTGTCAAGTATTAATGTATACTAAACATACTGTATGTACTTAATATCTATTGAATAAATGAATGAGTTTTTAATCTGTTTTAGTGGATTGCAAATGTACTTAATCTTTTCATTTGAATGACAAGAAGTTCTCAAATGTAATGAATAAAAGCGGTCCACAGTATTATTTTGTCAGAGGCAGTTTTCTTTCTTTAAAGTTTAATTGGATTTTTGTTTATGGTAATACAAACACATAATTAAAAAGTAGTTGATTTAAAAACGCATTTTATTATGCCCTCAAGAATAAATTCATTTAACATTGTAATATTTTTAAAATATTTTCCATAATACTAGGAGGTAGGAGAATATTTAATTAGAATTTATTAAACATGATTGAGAAGATGTGAAACTTTTGCAAATTCTTAACTTGAACTATTAAGAATAAAAAAGGAAAGGAAAGAACAAAAAAAAAAAAAAAACACCAAGTCGCACATCATTTGTTGTAGACCATCGCCACCTTGTGTTTGGAGAAAGATGTACTAAAGTCTGAGAATCTTTTAAGGTGGAAATTTTTTATTGCTGGAAAATTCCTATTTATAATTAAAACAAAACTTTACTGGAAACAACATCATTCCAAAGTCTGACAATATAATTACCTCCTACACAATTAGATAGGTTAAAATGTCAAAACAGTTGATTAAAAGTTATGAGTCAGAGTCATTTCTAAGCCCCCAGTTTGACTCTCGCTATGACATCTTGGGAAAACTACTTAACCCTTTTGAGTCTCATTTCTGCATCAGTTACACAAAGGGAGTCGATCTGGCAGATCATACTCCTCTGAGTTCTGTGAGGGCATGATGTGGAGGGAATATCAATTTTGAGTCTAAATGCAAGCTACAATAATTTTTCTGTCTCCAATATTTGCAATACTTCTGGGTATATAACCTGTTTAAACTAGAAGTTTTCTAACTACTTTCTCATTTGTGATCAATTAAAATATCCATTTGTCCAATATTACTTTAGACCATATGTTCCCAAAAGTAGGACTAATACCAAACAGCTTCCCTATTTCTAATTAGTCAATCTTATGGTACTAAAATGTGAATAAAAAGCATTTCCCGGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCCAGGGGGCGGAGCCTGCAGTGAGCCGAGATTGCGCCACTGCACTCCAGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCATTTCCCGGTTATGAAATATGCACATTAGTTGTTTCTCTAGACAATGTTTACATTAAAACAATTCTAATAGCTCATGTATATTACCAGATTTCCTCATTTCATGGCTACTTGTTGAATACTTACTTTGGTATAGACGTCGTGAAAATTGTGAATAGAATATTGTGTCTATAACATATGAATCAGTGAAACAAATAAGAACACATAAGGTGTAAATAGAAATTGACACCTTGATGCTTATCTCTCACTTTAAAATGCAGACAACTTGTTCTTTTTTTCAGACAGTCCCTCTGTTCATAGAGGACAATAATGTGTTTAATAGGTCAAACTGTTAATTGTTCCAATAATTTTATGTCAAATGTTTTGAGTTTTGTTCTGTTATCTCTAGAAATAAAAATTAGCTTACATATGTAAATTAGAACATGTTGAGATTAATGTTTATATAAGCCCATATGTATTCTATAAAAATAAAACCCAAAATATGGGGAGTAGGGATTAGGATAGCATGTGTTATTTTGGGTGCCCATTCAAATCAAGGAGCAAAAATTGGATTATGAAACTCACTACAAGAACCCCATATTACATGTCACTTTTACTAGTGCCTTGTACAAATTTACTCCTGTGGTCTTGTTTCTACATTCAAATGTTTTTCCAATTCACATTTTGCTTATAAATATACATATATGTCCTATGTATATATTTGCATATTACCACAAATGCCTGATGAGTGCTAATCTTTAAATAAAAGAGAAAAAGTAAGTGGATTAAATGTTCTTCTGAGCATAGATGGAGTTTTTCAGCTTAGTAGTTGGCAACTCCAACATCCATGTCTCTGTTGAAGCGATTAAAAAAGAATCATTTCCCAGAACTACAGGTGCACAACAGATTGGAGTCATTATGCAAAAGAACAGAAATCTGGATAGAAACCCAGAGAATAAAAAAGCCATGTTCCCGGCCATCCTTTTGCCTTTAACTCACTGTTTTTTCAACCCTACTTGTATCTGCAACATGAAATATGAAACAATTATAATGGATTTGTAGTAACCCATTCAAGATCAGAAATGTATGCTTTTAAACAAACCTCCAATAAGAATACCCTGCAGGTAGAAACACACTAAAGGTCAAACTTACGTTAACATATTTTTAAAAGCTAAAGTTTATGTAGAGAAGTTGATAGACTGAAATACATTACATATGTCTTTAAAGATACAGATTTTCTGCTGTTTTTGTCCCTTTTTGCCCCAATTTAAGGGAATCACATCTACCTTCACTCATCTTCATTCAGCCTCTGGGCCACCTCTCCTTGAGATGGAGATTATGAGTGTCCAAAATCTCCATCTCAAATCTCCATCTCAAAGCCACTCTGAGGCTGTAACTGTTGTCACCATA".getBytes();
            insertedSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC";
            complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", insertedSeq,
                    new SimpleInterval("21:26849022-26849037"), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("16M", "16M"));
            final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_EXPANSION, altSeq);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(14, 310, 26849037, 26849022, 1366, 1677);
            region1 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t2064\t21\t26849022\t60\t1704H657M3I706M\t*\t0\t0\tCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26847644,-,1394M1676S,60,1;chr10,97348533,-,1388S317M1365S,0,0;\tMD:Z:1363\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:1344\tXS:i:0", true);
            region2 = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00001\t16\t21\t26847644\t60\t1394M1676S\t*\t0\t0\tTATGGTGACAACAGTTACAGCCTCAGAGTGGCTTTGAGATGGAGATTTGAGATGGAGATTTTGGACACTCATAATCTCCATCTCAAGGAGAGGTGGCCCAGAGGCTGAATGAAGATGAGTGAAGGTAGATGTGATTCCCTTAAATTGGGGCAAAAAGGGACAAAAACAGCAGAAAATCTGTATCTTTAAAGACATATGTAATGTATTTCAGTCTATCAACTTCTCTACATAAACTTTAGCTTTTAAAAATATGTTAACGTAAGTTTGACCTTTAGTGTGTTTCTACCTGCAGGGTATTCTTATTGGAGGTTTGTTTAAAAGCATACATTTCTGATCTTGAATGGGTTACTACAAATCCATTATAATTGTTTCATATTTCATGTTGCAGATACAAGTAGGGTTGAAAAAACAGTGAGTTAAAGGCAAAAGGATGGCCGGGAACATGGCTTTTTTATTCTCTGGGTTTCTATCCAGATTTCTGTTCTTTTGCATAATGACTCCAATCTGTTGTGCACCTGTAGTTCTGGGAAATGATTCTTTTTTAATCGCTTCAACAGAGACATGGATGTTGGAGTTGCCAACTACTAAGCTGAAAAACTCCATCTATGCTCAGAAGAACATTTAATCCACTTACTTTTTCTCTTTTATTTAAAGATTAGCACTCATCAGGCATTTGTGGTAATATGCAAATATATACATAGGACATATATGTATATTTATAAGCAAAATGTGAATTGGAAAAACATTTGAATGTAGAAACAAGACCACAGGAGTAAATTTGTACAAGGCACTAGTAAAAGTGACATGTAATATGGGGTTCTTGTAGTGAGTTTCATAATCCAATTTTTGCTCCTTGATTTGAATGGGCACCCAAAATAACACATGCTATCCTAATCCCTACTCCCCATATTTTGGGTTTTATTTTTATAGAATACATATGGGCTTATATAAACATTAATCTCAACATGTTCTAATTTACATATGTAAGCTAATTTTTATTTCTAGAGATAACAGAACAAAACTCAAAACATTTGACATAAAATTATTGGAACAATTAACAGTTTGACCTATTAAACACATTATTGTCCTCTATGAACAGAGGGACTGTCTGAAAAAAAGAACAAGTTGTCTGCATTTTAAAGTGAGAGATAAGCATCAAGGTGTCAATTTCTATTTACACCTTATGTGTTCTTATTTGTTTCACTGATTCATATGTTATAGACACAATATTCTATTCACAATTTTCACGACGTCTATACCAAAGTAAGTATTCAACAAGTAGCCATGAAATGAGGAAATCTGGTAATATACATGAGCTATTAGAATTGTTTTAATGTAAACATTGTCTAGAGAAACAACTAATGTGCATATTTCATAACCGGGAAATGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAGGCTCCGCCCCCTGGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCTCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCGGGAAATGCTTTTTATTCACATTTTAGTACCATAAGATTGACTAATTAGAAATAGGGAAGCTGTTTGGTATTAGTCCTACTTTTGGGAACATATGGTCTAAAGTAATATTGGACAAATGGATATTTTAATTGATCACAAATGAGAAAGTAGTTAGAAAACTTCTAGTTTAAACAGGTTATATACCCAGAAGTATTGCAAATATTGGAGACAGAAAAATTATTGTAGCTTGCATTTAGACTCAAAATTGATATTCCCTCCACATCATGCCCTCACAGAACTCAGAGGAGTATGATCTGCCAGATCGACTCCCTTTGTGTAACTGATGCAGAAATGAGACTCAAAAGGGTTAAGTAGTTTTCCCAAGATGTCATAGCGAGAGTCAAACTGGGGGCTTAGAAATGACTCTGACTCATAACTTTTAATCAACTGTTTTGACATTTTAACCTATCTAATTGTGTAGGAGGTAATTATATTGTCAGACTTTGGAATGATGTTGTTTCCAGTAAAGTTTTGTTTTAATTATAAATAGGAATTTTCCAGCAATAAAAAATTTCCACCTTAAAAGATTCTCAGACTTTAGTACATCTTTCTCCAAACACAAGGTGGCGATGGTCTACAACAAATGATGTGCGACTTGGTGTTTTTTTTTTTTTTTTGTTCTTTCCTTTCCTTTTTTATTCTTAATAGTTCAAGTTAAGAATTTGCAAAAGTTTCACATCTTCTCAATCATGTTTAATAAATTCTAATTAAATATTCTCCTACCTCCTAGTATTATGGAAAATATTTTAAAAATATTACAATGTTAAATGAATTTATTCTTGAGGGCATAATAAAATGCGTTTTTAAATCAACTACTTTTTAATTATGTGTTTGTATTACCATAAACAAAAATCCAATTAAACTTTAAAGAAAGAAAACTGCCTCTGACAAAATAATACTGTGGACCGCTTTTATTCATTACATTTGAGAACTTCTTGTCATTCAAATGAAAAGATTAAGTACATTTGCAATCCACTAAAACAGATTAAAAACTCATTCATTTATTCAATAGATATTAAGTACATACAGTATGTTTAGTATACATTAATACTTGACAATCAATACTGGTTAACTGGTTTCCCTGGTTTAGAAATTTTCCTTAGCAACAACGTAAGGCTTAAAATGAAAAAAGAAAAGTGAGAAAATGTTCTACCACCAGGTGGTGACAAAAGATAAAATTTAAAATCGCTCTTAATGAGCACATACTTCATGTAATTCTTGAATACTGCAAATATAAGTGACTTCCGAATGTCATGTGAATTTAAAATCATATTCTAGGAATATTTTATTAATTAAAGCAAATTAATATTAACATATTATCTCT\t*\tSA:Z:21,26849022,-,1704S657M3I706M,60,3;chr10,97348533,-,1388S317M1365S,0,0;\tMD:Z:1204A189\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1384\tXS:i:0", true);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpointsDetectedFromReverseStrand, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        {
            // simple tandem duplication expansion with novel insertion '+' strand representation
            final byte[] leftRefFlank = "GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT".getBytes();                     //63
            final byte[] rightRefFlank = "CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT".getBytes();                                   //48
            final byte[] insertedSeq = "CTCTCTCTCT".getBytes();                                                                           //10
            final byte[] dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG".getBytes();    //89
            outputStream.reset();
            outputStream.write(leftRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(rightRefFlank);
            byte[] contigSeq = outputStream.toByteArray();
            String contigName = "simple_dup_exp_1_2_with_ins_+";

            final SimpleInterval leftBreakpoint = new SimpleInterval("21", 25297163, 25297163);
            final SimpleInterval rightBreakpoint = new SimpleInterval("21", 25297163, 25297163);
            final BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications("", new String(insertedSeq),
                    new SimpleInterval("21", 25297164,25297252), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList("89M", "89M"));
            final byte[] altSeq = Arrays.copyOfRange(contigSeq, leftRefFlank.length, contigSeq.length - rightRefFlank.length);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_EXPANSION, altSeq);
            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 1 ,152, TextCigarCodec.decode("152M147S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 163 ,299, TextCigarCodec.decode("162S137M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-89, 10, 25297252, 25297164, 152, 163);
            final List<SvType> type = Collections.singletonList(makeTandemDuplicationType());
            final List<VariantContext> variant = Collections.singletonList(
                    addStandardAttributes(makeTandemDuplication(), )
                            .attribute().attribute().attribute().attribute().make()
            );
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));

            // simple tandem duplication expansion with novel insertion '-' strand representation
            SequenceUtil.reverseComplement(leftRefFlank);
            SequenceUtil.reverseComplement(rightRefFlank);
            SequenceUtil.reverseComplement(insertedSeq);
            SequenceUtil.reverseComplement(dup);
            outputStream.reset();
            outputStream.write(rightRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(leftRefFlank);
            contigSeq = outputStream.toByteArray();
            contigName = "simple_dup_exp_1_2_with_ins_-";

            region1 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 1 ,137, TextCigarCodec.decode("137M162S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 148 ,299, TextCigarCodec.decode("147S152M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-89, 10, 25297252, 25297164, 137, 148);
            result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithPreciseDupRangeBreakpointsInference.class));
        }

        return result;
    }

    /**
     * These test data was based on a real observation on a locally-assembled contig
     * "TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC"
     * with two alignment records chr18:312579-312718 140M135S
     *                            chr18:312610-312757 127S148M
     * for a tandem repeat expansion event from 1 copy to 2 copies with also a pseudo-homology

     * Return a list of eight entries for positive and reverse strand representations for:
     * 1. expansion from 1 unit to 2 units with pseudo-homology
     * 2. contraction from 2 units to 1 unit with pseudo-homology
     * 3. contraction from 3 units to 2 units without pseudo-homology
     * 4. expansion from 2 units to 3 units without pseudo-homology
     */
    private static List<TestDataForSimpleSV>
    forComplexTandemDuplication(final ByteArrayOutputStream outputStream) {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        String contigName = "cpx_dup_exp_1_2_pseudo_+";
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        SimpleInterval leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        SimpleInterval rightBreakpoint = new SimpleInterval("20", 312609, 312609);
        BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(pseudoHomology, "",
                new SimpleInterval("20", 312610, 312705), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312718));
        byte[] altSeq = (firstRepeat+secondRepeat+pseudoHomology).getBytes();
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1 ,140, TextCigarCodec.decode("140M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 128 ,275, TextCigarCodec.decode("127S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-109, -13, 312718, 312610, 140, 128);
        List<SvType> type = Collections.singletonList(makeTandemDuplicationType());
        List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeTandemDuplication(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_1_2_pseudo_-";
        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 1 ,148, TextCigarCodec.decode("148M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 136 ,275, TextCigarCodec.decode("135S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-109, -13, 312718, 312610, 148, 136);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        contigName = "cpx_dup_contract_2_1_pseudo_+";
        final byte[] contigSeqForComplexContractionWithPseudoHomology = fakeRefSeqForComplexExpansionWithPseudoHomology;
        leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        rightBreakpoint = new SimpleInterval("20", 312705, 312705);
        complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat+pseudoHomology, "",
                new SimpleInterval("20", 312610, 312705), 2, 1, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Collections.singletonList(Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312814));
        altSeq = (firstRepeat+pseudoHomology).getBytes();
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1, 140, TextCigarCodec.decode("140M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 32, 179, TextCigarCodec.decode("31S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-13, -109, 312718, 312706, 140, 32);
        type = Collections.singletonList(makeDeletionType());
        variant = Collections.singletonList(
                addStandardAttributes(makeDeletion(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionWithPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_contract_2_1_pseudo_-";
        final byte[] contigSeqForComplexContractionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionWithPseudoHomology, contigSeqForComplexContractionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionWithPseudoHomology_reverseStrand);
        breakpoints =  new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 1, 148, TextCigarCodec.decode("148M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 40, 179, TextCigarCodec.decode("39S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-13, -109, 312718, 312706, 148, 40);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionWithPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // third test: contraction from 3 units to 2 units without pseudo-homology
        contigName = "cpx_dup_contract_3_2_+";
        final byte[] fakeRefSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, rightRefFlank).getBytes();
        leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        rightBreakpoint = new SimpleInterval("20", 312705, 312705);
        complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat+secondRepeat, "",
                new SimpleInterval("20", 312610, 312705), 3, 2, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312897));
        altSeq = (firstRepeat+secondRepeat).getBytes();
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 32, 262, TextCigarCodec.decode("31S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-96, -192, 312801, 312706, 223, 32);
        type = Collections.singletonList(makeDeletionType());
        variant = Collections.singletonList(
                addStandardAttributes(makeDeletion(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionNoPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_contract_3_2_-";
        final byte[] contigSeqForComplexContractionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionNoPseudoHomology, contigSeqForComplexContractionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionNoPseudoHomology_reverseStrand);
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 1, 231, TextCigarCodec.decode("231M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 40, 262, TextCigarCodec.decode("39S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-96, -192, 312801, 312706, 231, 40);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexContractionNoPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        contigName = "cpx_dup_exp_2_3_+";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = fakeRefSeqForComplexContractionNoPseudoHomology;
        leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        rightBreakpoint = new SimpleInterval("20", 312609, 312609);
        complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat, "",
                new SimpleInterval("20", 312610, 312705), 2, 3, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312801));
        altSeq = (firstRepeat+secondRepeat+firstRepeat).getBytes();
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 128, 358, TextCigarCodec.decode("127S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-192, -96, 312801, 312610, 223, 128);
        type = Collections.singletonList(makeTandemDuplicationType());
        variant = Collections.singletonList(
                addStandardAttributes(makeTandemDuplication(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_2_3_-";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 1, 231, TextCigarCodec.decode("231M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 136, 358, TextCigarCodec.decode("135S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-192, -96, 312801, 312610, 231, 136);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        return result;
    }

    /**
     * See {@link #forComplexTandemDuplication(ByteArrayOutputStream)} .
     * Here we are simply making the repeat sequence shorter than 50 bases so that the end type will be INS instead of DUP
     */
    private static List<TestDataForSimpleSV>
    forComplexTandemDuplicationIns(final ByteArrayOutputStream outputStream) {

        final List<TestDataForSimpleSV> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";              // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";      // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAA";   // 42
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAA";   // 42
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                // 13


        // first test : expansion from 1 unit to 2 units with pseudo-homology
        String contigName = "cpx_dup_exp_small_1_2_pseudo_+";
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        SimpleInterval leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        SimpleInterval rightBreakpoint = new SimpleInterval("20", 312609, 312609);
        BreakpointComplications complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(pseudoHomology, "",
                new SimpleInterval("20", 312610, 312651), 1, 2, Collections.singletonList(Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312664));
        byte[] altSeq = String.format("%s%s%s", firstRepeat, secondRepeat, pseudoHomology).getBytes();
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312664), 1 ,86, TextCigarCodec.decode("86M81S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312703), 74 ,167, TextCigarCodec.decode("73S94M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        DistancesBetweenAlignmentsOnRefAndOnRead distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, -13, 312664, 312610, 86, 74);
        final List<SvType> type = Collections.singletonList(makeInsertionType());
        final List<VariantContext> variant = Collections.singletonList(
                addStandardAttributes(makeInsertion(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_small_1_2_pseudo_-";
        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312703), 1 ,94, TextCigarCodec.decode("94M73S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312664), 82 ,167, TextCigarCodec.decode("81S86M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-55, -13, 312664, 312610, 94, 82);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        // second test: expansion from 2 units to 3 units without pseudo-homology
        contigName = "cpx_dup_exp_too_small_2_3_+";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        leftBreakpoint = new SimpleInterval("20", 312609, 312609);
        rightBreakpoint = new SimpleInterval("20", 312609, 312609);
        complications = new BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications(firstRepeat, "",
                new SimpleInterval("20", 312610, 312651), 2, 3, Arrays.asList(Strand.POSITIVE, Strand.POSITIVE), Arrays.asList(Strand.POSITIVE, Strand.POSITIVE, Strand.POSITIVE),
                new SimpleInterval("20", 312610, 312693));
        altSeq = String.format("%s%s%s", firstRepeat, secondRepeat, firstRepeat).getBytes();
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312693), 1, 115, TextCigarCodec.decode("115M81S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312732), 74, 196, TextCigarCodec.decode("73S123M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-84, -42, 312693, 312610, 115, 74);
        type = Collections.singletonList(makeInsertionType());
        variant = Collections.singletonList(
                addStandardAttributes(makeInsertion(), )
                        .attribute().attribute().attribute().attribute().make()
        );
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        contigName = "cpx_dup_exp_too_small_2_3_-";
        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        breakpoints = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, NO_SWITCH, complications, SMALL_DUP_CPX, altSeq);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312732), 1, 123, TextCigarCodec.decode("123M73S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312693), 82, 196, TextCigarCodec.decode("81S115M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, region1, region2, NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        distances = new DistancesBetweenAlignmentsOnRefAndOnRead(-84, -42, 312693, 312610, 123, 82);
        result.add(new TestDataForSimpleSV(region1, region2, contigName, contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, simpleChimera, breakpoints, type, variant, distances, BreakpointsInference.SmallDuplicationWithImpreciseDupRangeBreakpointsInference.class));

        return result;
    }
}
