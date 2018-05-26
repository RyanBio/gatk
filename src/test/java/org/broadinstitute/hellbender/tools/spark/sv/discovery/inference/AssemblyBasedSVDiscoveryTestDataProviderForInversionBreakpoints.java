package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_55;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * Extracted here because inversions are treated differently.
 */
public class AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints extends AssemblyBasedSVDiscoveryTestDataProvider {

    public static final class TestDataForInversion extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera {

        private TestDataForInversion(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                     final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                     final SimpleChimera manuallyCuratedSimpleChimera,
                                     final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble,
                                     final List<SvType> manuallyCuratedSVTypes,
                                     final List<VariantContext> manuallyCuratedVariants) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, manuallyCuratedSimpleChimera, manuallyCuratedBiPathBubble, manuallyCuratedSVTypes, manuallyCuratedVariants,
                    BreakpointsInference.IntraChrStrandSwitchBreakpointInference.class);
        }

        @Override
        public SAMSequenceDictionary getAppropriateDictionary() {
            return TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return BreakpointsInference.IntraChrStrandSwitchBreakpointInference.class;
        }
    }

    // the chromosome that the long contig1 is supposed to be mapped to is actually chr19, but to make tests runnable, we could only use "20" or "21"
    // todo: this should be fixed, but since the exact mapped to chromosome is not important now, we push it to later
    public static final String chrForLongContig1 = "20";
    public static final String homologyForLongContig1 = "TTTCTTTTTTCTTTTTTTTTTTTAGTTGATCATTCTTGGGTGTTTCTCGCAGAGGGGGATTTGGCAGGGTCATAGGACAACAGTGGAGGGAAGGTCAGCAGATAAACAAGTGAACAAAGGTCTCTGGTTTTCCTAGGCAGAGGACCCTGCGGCCTTCCGCAGTGTTTGTGTCCCCGGGTACTTGAGATTAGGGAGTGGTGATGACTCTTAACGAGCATGCTGCCTTCAAGCATCTGTTTAACAAAGCACATCTTGCACCGCCCTTAATCCATTTAACCCTGAGTGGACACAGCACATGTTTCAGAGAGCACAGGGTTGGGGGTAAGGTCATAGGTCAACAGGATCCCAAGGCAGAAGAATTTTTCTTAGTACAGAACAAAATGAAGTCTCCCATGTCTACCTCTTTCTACACAGACACAGCAACCATCCGATTTCTCAATCTTTTCCCCACCTTTCCCCCGTTTCTATTCCACAAAACTGCCATTGTCATCATGGCCCGTTCTCAATGAGCTGTTGGGTACACCTCCCAGACGGGGTGGTGGCCGGGCAGAGGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCACCCCCCACCTCCCGGACGGGGCGGCTGGCTGGGCGGGGGGCTGACCCCCCCCACCTCCCTCC";
    public static final String LONG_CONTIG1 =
                    "TTTTTTTTTTTTTTTCTGAGACCGAATCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTTACTGCAAGCTCTGCCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCCCACCCCCCCACCCCCCCAGGTAGCTG" +
                    "GGACTACAGGTGTCTGCCACCACACCTGGCTAAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAAGATGGTTTCGCTCTCCTGACCTCGCGATCCGCCCACCTCGGCCTCTCAAAGTGCTGGGATTACAGGCCTGAGCCACTGCGCCC" +
                    "TGCCTGACCCCTTCTTTTAAAACAAATCTTTTGTCTTTGTCTTCATTTCTGCATTCGTCCCCTTCGTTCAATCCTGTAGGAATTGACAGTGATATTGGGGAGTTCCCATCCCTGGATTTGGGATTTCCTCGAGTTTCCAGCCCTGTCCTTGTGGCCAAAAAGT" +
                    "GGCACAGTTTGGTTTTTGAATCCTGTTCGTGACGCCAAAAATATTCTGTGTGGGAAACATTCAAGGGGAGAAGAAAAGACACACACACAATACCTTTAAGGGTAAATAAGCTTTGTGCCACGTAAATGGAAATGCAGATATAATTAGCAAATTATATAATAAG" +
                    "CAAATCAATGTAATAACCAGATTGATATAATAAGCATGTTGATATAATAAGCAAATTGCAATGGGAAGAGGAGAAAGGAAAAGAGATATATATATTTACACTCACCAGACTATGGAGGATTCACCACCAGACTGGGAAGCAACAGCTTGGGCTCCAGAGTCAG" +
                    "CTACTCATTCCTGCACAGATGAGGAGGGTCTAATGAAGCTTCAGCACAATCTGATACCCTAGCTCTTTTGTAACGAGTTGTTTGGCATAAGGCCCAGTCATGAGGGCCATTTGCAACTGGGGTCAAGGAACACAAAACTGTCAACTTGTTTTTGCGATTGTCT" +
                    "ATTGTTTTTCAACAACTAATATATAGAAATAGATTGAAATAGAGATTTCTCTGAAACAGCGCTGAATGAATGCCTCAAGGGCCTCACACAACCTGTTTGAGAACTTGGTGACTACTGAGTTTGTCCACGTTCAATTAAGTTCAAATTTAGTATTTAACTTTTC" +
                    "CTCCACAAATTACCCAGTCTTGATGAATCAGCTCTGTCTAGACATAGGGCAAGGTGAACCCCTTGGGCAGTTACACAACCTCCGCCTTCTGGGTTTAAGCAATTCTCCTGCCTCAGCCTCCGGACTAGCTGGGTCTACAGGTGTGCAGCACCACACCCAGCTA" +
                    "GTTATTTGTACTTTTAGTAGAAATGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCACCTTGGCCTCCCAAAGTGTTGCGATTACAGGCACTTGCCAGTGAACCTGGCCCTAAATGACTTCTTTCTATCTCCTATAA" +
                    "CAATTTGAAATTACTTAAAGGTGGTTTCAAATTGAAAAAATAAAAAGAATTTGGATAAAAATAAAATGTAAACAGTTTTTAAAAATTACAAGAGATTACAAAATATACATGTAAAACCGGAGTGGTCAAAAATGACAAATTTGATTTATTTATAAGGTTTATT" +
                    "AAAATTAGCTTTAGTATTGATAATACACTATTACCAAAGTAAAAGCTGATTTTCTCTTGAAAAAAATTTTATGTATTATTAATATGACAGCAAAATACTTCTGTTCACCTTTTGAATATATTCAAAAAGAGAGAGAGTAAAAAAAAAGATAGAATTTTCCCAT" +
                    "GATCTGGGTTGGGCCTGGCTTAGCTCAGGGAGGAAGCCCTGCCTGAAAAACGCTGCAGCTTAGGCTGTGACTTTTTCTTCACTCAGTCCTGATCACATTTTCTGTTACTCAGGGCCTGAGGGGGCGGGGGCCTTAAGCATTATCCAATCAGAAACGCTGGGCT" +
                    "GACGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATT" +
                    "GAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATC" +
                    "CTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTA" +
                    "CCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAA" +
                    "AAAAGAAAAAAGAAAAAAGAAAAAATACACCCTGGCAGAAAACTTCTCGAATGTAAGTAAAGTGGGAGTTTTCAACTTGCCCACCTCCTTACAAAAACATGCAATAAGAATGTGCTAGATATAAAGCTTACAAATGTGAGAAATACGGGAAGTTTTTCAATTT" +
                    "TAACAAATACTTTCAAAGTCATGTGAAAATTCTCACTAGAAAGAAATCCAATAAAGGCAAGTTGAGGCATAAAAAAAAAAAAAAGAAAAAGCTCAGAATACATATAAGATTGACTAAGTCAGCCAGAAAATAATCCCCTAAAAGAAATTTCTCTCTAAACACC" +
                    "CAATGTGCACAGCTACTCTCAGCATGAGAAACATGAGCTTTATGAAGAAAGGTGGCAGATTTTCAGAAGAATTTTATAAAAGTTTCTTTTCCATCTCTGCTGCTCTCTCATCTCCTAGCCATTGAATGGGGGTTCTATATTGAAATACATCTGACAACTTCCA" +
                    "ACAACACTTTTTGATCAAGAAATAGAATTTGACTATGTTCGTATAGTGGAATATATTAGAACTTGTAACACAGCTAACTGAATAGCTATTATGGTGTTTGGGTGGCCACATCACCTGTCTTTATTTGTCCGGTAATAGCAGCATTCCAATTTAAAGAAATAAA" +
                    "AGATACCAAAATTGTGTTTACTTTTAATTATTCCTATTGAATAAAGTAATAAGCATGTCAGACTGATATCTATCATAACAATAAATTTTGTTTGGATATTATATTAGATATAAATATTTAAGTATGAATAATTTTAATGAACTAGTCATAATGTATGTAGCAT" +
                    "TTTTAAAAATTGTAACTATACTTCAGTTAAAACACTTTATATTTCAAAAGCATAAATAACAATATTAAAATAACAATTTAGGTGATTCATTCAAAGTAAGTATTGGGGCTTTATGTTCATACTATTGTAGAAAATACTGCTTATGGCTCATGCCTGTAATCCC" +
                    "AGCACATTGGGAGGCTGAGGTGGGTAGATCACCTGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGG" +
                    "TTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCA" +
                    "CAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTAT" +
                    "CTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTT" +
                    "TTTTTTTTTTGAGATGGAGTCTCGCTCTGTCACCGAGGCTGGAGTGCAGTGGCACGAACTCGGCTCACTGCAAAAGCTCTGCCTCCCAGCTTCATGCCATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGATTACAGGCAACTGCCACCACGCCCAGCTAATT" +
                    "TTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTTGATCTCCTGACCTCATGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGCAAAACCGCGCCCTGCCCTTATTTAGCTCTTAATATGATTTACATATAT" +
                    "TTCAAAAAAGCAGAGAAAAATATCTACATATAATCTAAATCCCTTAAGAGAGAATAGCAAAAAAATTTTTGGAACTTTTTAAAGAGTTTTGAACTCTTGGACATCTGAATTTTGCACACTATATGCACTTGAAAGAATGTTAATGGGGAAAAAGCAGAAGAGA" +
                    "AAAAGACGTTATAAAAAAATCCATGAGTGCACAACACCTATGAAACAGAATAGAGAGCCCAGAAATAATGCCTTTCACCTGCAACCATCAGATTTCTGACAAAGCTGACAAGAGGAATGTGGGAAGAATTCTCTCTTTCATAAATGGTGCTGGAATAACTATC" +
                    "TACCACTATGTAGAAGACTGAAGTGGACCCCTTCATTACACCATATAAAAAAATCAACTGAAGATAAATTAAGGACTTAAATGTAAAACTTAAAATTATAAGAAACCCTGCAAGATAACCTAGGAAATAGCATTCTAGACACAGAAACAGGTAAAGACTTCAT" +
                    "GATGAAGCTACCAAAAGCAACTGCAACAGAAGTAAATTGACAAATGGGATGTATTTAAACTTAAGAGCTTCTTCACAGCAAAGGAAACTATCAACAGAGTAAACAGACAAACTAGAGAATAAAAGAATATATTTGTAAATTTTGCCTCTGAAAAAGGTCTAAT" +
                    "ATACAGAATTTATTAGGAACTTAAACAAGTTTACAAGAAAAAAACACACTCATTAAAAAGTATGCAAAAAACATGAACAGATGCCTTTCAATAGAAGATACACATAGGGCTAACAAGCATATGAAAAAAAAAATGCTCATTGCTAATCATTAGAGAAATTAGA" +
                    "AGAGAAATCACAGGCTGGGTGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATTACAAGGTCAGGAGATCAAGATCATCCTGGCTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATCAGCCAGGTGTGGCAGTGT" +
                    "GCACCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGTAGGCAGAGGTTGCAGTGAGCTGAGATCACACCACTGCACTCCTGCCTGGGCAACAGAGCAAGACTCCGTCTCAAACACACACACACAGACACACACACACACACAC" +
                    "ACACACACACACACACACACGCAGAGAAACCACAATGAGATACCACCTCATACCAGTCAGAATGGCTATTTTTAAAAAGTCAAAAGATAACAGATGCTGACAAAGTTGCAGAGAAAAGGGAATGCTTATACTCCTCTGGTGGGAGTGTAAATTATTTCAACAA" +
                    "CTGTAAAAATCAGTGTGACAATTCCTCACAGAATGAAAAACAGAATTATCATTCGACTGGGAAACCTCATAATTGAGTATATACCCAAAGAAATATGAAATATTATAAAGACACATCCACATGCATGTTCACTGCAGCACTATTCACAATAGCAAAGACACGG" +
                    "ACAGACTAAATGCCTATCAATGGCAGACTGGATCAAGAAAATATGGTATGGTCAGATGCGGTGGCTCATGCCTGTAATTCCAGCCCTTTGGGAGGCTGAGGCAGGTGGATTGCCTGAGCTTAGAAGTTTGAGACCACTCTGGGCAACATGGCAAAATTTTGTC" +
                    "TCCACAGAAGATACAAAAAAAAAAAAAAAAAA";
    public static final TestDataForInversion forSimpleInversionWithInsertion_LeftBreakpoint_plus;
    public static final TestDataForInversion forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
    public static final TestDataForInversion forSimpleInversionWithHom_leftPlus;
    public static final TestDataForInversion forSimpleInversionWithHom_leftMinus;
    public static final TestDataForInversion forSimpleInversionWithHom_rightPlus;
    public static final TestDataForInversion forSimpleInversionWithHom_rightMinus;

    public static final boolean testDataInitialized;
    static {
        try (final ByteArrayOutputStream outputStream = new ByteArrayOutputStream()) {

            forSimpleInversionWithInsertion_LeftBreakpoint_plus = forSimpleInversionWithInsertion_LeftBreakpoint_plus();
            forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint = forSimpleInversionWithHomology_RightBreakpoint_minus();
            final List<TestDataForInversion> inversion4 = forSimpleInversionWithHomology(outputStream);
            forSimpleInversionWithHom_leftPlus = inversion4.get(0);
            forSimpleInversionWithHom_leftMinus = inversion4.get(1);
            forSimpleInversionWithHom_rightPlus = inversion4.get(2);
            forSimpleInversionWithHom_rightMinus = inversion4.get(3);
            testDataInitialized = true;
        } catch (final Exception ioex) {
            throw new GATKException("Failed to create test data ", ioex);
        }
    }

    static List<TestDataForInversion> getAllTestData() {
        final List<TestDataForInversion> testDataForInversions = Arrays.asList(AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithInsertion_LeftBreakpoint_plus,
                AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint,
                AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_leftPlus,
                AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_leftMinus,
                AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_rightPlus,
                AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_rightMinus);
        return Collections.unmodifiableList(testDataForInversions);
    }

    static List<Tuple2<TestDataForInversion, TestDataForInversion>> getAllTestDataPaired() {
        final List<TestDataForInversion> allTestData = getAllTestData();
        final List<Tuple2<TestDataForInversion, TestDataForInversion>> testDataForSimpleSVs
                = new ArrayList<>(allTestData.size()/2);
        for (int i = 0; i < allTestData.size() - 1; i += 2) {
            testDataForSimpleSVs.add(new Tuple2<>(allTestData.get(i), allTestData.get(i+1)));
        }
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    private static TestDataForInversion
    forSimpleInversionWithInsertion_LeftBreakpoint_plus() {
        // inversion with inserted sequence
        String contigName = "simple_inv_55_with_ins_+";
        final byte[] leftFlank = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(146, (byte)'A');
        final byte[] rightFlankRC = TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(50, (byte)'C');
        final byte[] contigSeq = new byte[leftFlank.length+1+rightFlankRC.length];
        System.arraycopy(leftFlank, 0, contigSeq, 0, leftFlank.length);
        contigSeq[leftFlank.length] = (byte) 'T';
        System.arraycopy(rightFlankRC, 0, contigSeq, leftFlank.length+1, rightFlankRC.length);

        // reference intervals are changed from real values, which were generated on a different chromosome, to a fake but valid value.
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 69149, 69294), 1, 146, TextCigarCodec.decode("146M51S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 69315, 69364), 148, 197, TextCigarCodec.decode("147S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSeq, Arrays.asList(region1, region2));
        final SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final SimpleInterval leftLocation = new SimpleInterval("21", 69294, 69294);
        final SimpleInterval rightLocation = new SimpleInterval("21", 69364, 69364);
        final BreakpointComplications complications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications("", "T");

        final NovelAdjacencyAndAltHaplotype biPathBubble = new NovelAdjacencyAndAltHaplotype(leftLocation, rightLocation,
                StrandSwitch.FORWARD_TO_REVERSE, complications, INTRA_CHR_STRAND_SWITCH_55,  "T".getBytes());

        final List<SvType> breakEnds = Arrays.asList(
                makeBNDType("BND_21_69294_69364_1", Allele.create("AT]21:69364]"), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55),
                makeBNDType("BND_21_69294_69364_2", Allele.create("GA]21:69294]"), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55)
        );

        final List<VariantContext> variants = Arrays.asList(
                makeBND(new SimpleInterval("21:69294-69294"), new SimpleInterval("21:69364-69364"), Allele.create("A"), "T", INV55, true, true, true)
                        .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 50).attribute(MAX_ALIGN_LENGTH, 50).attribute(CONTIG_NAMES, contigName).make(),
                makeBND(new SimpleInterval("21:69294-69294"), new SimpleInterval("21:69364-69364"), Allele.create("G"), "A", INV55, false, true, true)
                        .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 50).attribute(MAX_ALIGN_LENGTH, 50).attribute(CONTIG_NAMES, contigName).make()
        );
        return new TestDataForInversion(region1, region2, contigName, contigSeq, simpleChimera, biPathBubble, breakEnds, variants);
    }

    private static TestDataForInversion
    forSimpleInversionWithHomology_RightBreakpoint_minus() {
        // inversion with strange left breakpoint
        final String contigName = "simple_inv_33_with_hom_-";
        final byte[] contigSequence = LONG_CONTIG1.getBytes(); // 5210 length
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 20262839, 20264443), 1, contigSequence.length - 3603, TextCigarCodec.decode("1970M1I611M1I24M3603H"), false, 60, 8, 2561, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 20248787, 20253040), 1954, contigSequence.length, TextCigarCodec.decode("1953S10M3I9M1I246M2D1572M1I798M5D730M1I347M4I535M"), true, 60, 41, 4068, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSequence, Arrays.asList(region1, region2));

        final SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final SimpleInterval leftLocation = new SimpleInterval("20", 20248787, 20248787);
        final SimpleInterval rightLocation = new SimpleInterval("20", 20263492, 20263492);
        final BreakpointComplications complications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(homologyForLongContig1, "");

        final NovelAdjacencyAndAltHaplotype biPathBubble = new NovelAdjacencyAndAltHaplotype(leftLocation, rightLocation,
                StrandSwitch.REVERSE_TO_FORWARD, complications, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33, EMPTY_BYTE_ARRAY);

        final List<SvType> breakEnds = Arrays.asList(
                makeBNDType("BND_20_20248787_20263492_1", Allele.create("[20:20263492[C"), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33),
                makeBNDType("BND_20_20248787_20263492_2", Allele.create("[20:20248787[A"), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33)
        );

        final List<VariantContext> variants = Arrays.asList(
                makeBND(new SimpleInterval("20:20248787-20248787"), new SimpleInterval("20:20263492-20263492"), Allele.create("C"), "", INV33, true, false, false)
                        .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 2607).attribute(MAX_ALIGN_LENGTH, 2607).attribute(CONTIG_NAMES, contigName).make(),
                makeBND(new SimpleInterval("20:20248787-20248787"), new SimpleInterval("20:20263492-20263492"), Allele.create("A"), "", INV33, false, false, false)
                        .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 2607).attribute(MAX_ALIGN_LENGTH, 2607).attribute(CONTIG_NAMES, contigName).make()
        );
        return new TestDataForInversion(region1, region2, contigName, contigSequence, simpleChimera, biPathBubble, breakEnds, variants);
    }

    /**
     * The following four tests are all going to be for the same inversion, testing if implementations are correct for
     * identifying the breakpoints by looking at different representations of evidence.
     * The inversion we are looking at is
     *
     * '+' strand representation: G....100....G|ACACA|C....100....C               A....100....A|TGTGT|T....100....T
     *
     * 100-bases of 'G' is the left flanking before the homology |ACACA| and the region starting with 100-bases of 'C' and
     * ending with 100-bases of 'A' and maybe (homologyForwardStrandRep uncertainty) the homology |TGTGT| is inverted.
     * 100-bases of 'T' is the right flanking region.
     *
     * Returns a list of four Tuple4's with left flanking evidence '+'/'-' strand representation and right flanking side.
     */
    private static List<TestDataForInversion>
    forSimpleInversionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForInversion> result = new ArrayList<>();

        final byte[] leftLeftPlus = makeDummySequence(100, (byte)'G');
        final byte[] leftLeftMinus = makeDummySequence(100, (byte)'C');
        final byte[] leftRightPlus = makeDummySequence(100, (byte)'C');
        final byte[] leftRightMinus = makeDummySequence(100, (byte)'G');
        final byte[] rightLeftPlus = makeDummySequence(100, (byte)'A');
        final byte[] rightLeftMinus = makeDummySequence(100, (byte)'T');
        final byte[] rightRightPlus = makeDummySequence(100, (byte)'T');
        final byte[] rightRightMinus = makeDummySequence(100, (byte)'A');
        final byte[] leftHomology = "ACACA".getBytes();
        final byte[] rightHomology = "TGTGT".getBytes();
        {// left breakpoint '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftLeftPlus);outputStream.write(leftHomology);outputStream.write(rightLeftMinus);
            byte[] contigSeq = outputStream.toByteArray();

            String contigName = "simple_inv_55_with_hom_+";

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 200101, 200205), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 200501, 200605), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            final SimpleInterval leftLocation = new SimpleInterval("20", 200,200);
            final SimpleInterval rightLocation = new SimpleInterval("20", 605, 605);
            final BreakpointComplications complications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(StringUtil.bytesToString(leftHomology), "");
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftLocation, rightLocation, StrandSwitch.FORWARD_TO_REVERSE, complications, INTRA_CHR_STRAND_SWITCH_55, EMPTY_BYTE_ARRAY);
            final List<SvType> breakEnds = Arrays.asList(
                    makeBNDType("BND_20_200200_200605_1", Allele.create("C]20:200605]"), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55),
                    makeBNDType("BND_20_200200_200605_2", Allele.create("T]20:200200]"), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_55)
            );
            final List<VariantContext> variants = Arrays.asList(
                    makeBND(new SimpleInterval("20:200200-200200"), new SimpleInterval("20:200605-200605"), Allele.create("C"), "", INV55, true, true, true)
                            .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 105).attribute(MAX_ALIGN_LENGTH, 105).attribute(CONTIG_NAMES, contigName).make(),
                    makeBND(new SimpleInterval("20:200200-200200"), new SimpleInterval("20:200605-200605"), Allele.create("T"), "", INV55, false, true, true)
                            .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 105).attribute(MAX_ALIGN_LENGTH, 105).attribute(CONTIG_NAMES, contigName).make()
            );

            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, breakEnds, variants));

            outputStream.reset();
            outputStream.write(rightLeftPlus);outputStream.write(rightHomology);outputStream.write(leftLeftMinus);
            contigSeq = outputStream.toByteArray();
            contigName = "simple_inv_55_with_hom_-";
            region1 = new AlignmentInterval(new SimpleInterval("20", 200501, 200605), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 200101, 200205), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.FORWARD_TO_REVERSE, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, breakEnds, variants));
        }
        {// right breakpoint '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftRightMinus);outputStream.write(rightHomology);outputStream.write(rightRightPlus);
            byte[] contigSeq = outputStream.toByteArray();

            String contigName = "simple_inv_33_with_hom_+";

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 200201, 200305), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 200601, 200705), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            SimpleChimera simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            final SimpleInterval leftLocation = new SimpleInterval("20", 201,201);
            final SimpleInterval rightLocation = new SimpleInterval("20", 606, 606);
            final BreakpointComplications complications = new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(StringUtil.bytesToString(leftHomology), "");
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(leftLocation, rightLocation, StrandSwitch.FORWARD_TO_REVERSE, complications, INTRA_CHR_STRAND_SWITCH_33, EMPTY_BYTE_ARRAY);
            final List<SvType> breakEnds = Arrays.asList(
                    makeBNDType("BND_20_200201_200606_1", Allele.create("[20:200606[A"), true, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33),
                    makeBNDType("BND_20_200201_200606_2", Allele.create("[20:200201[G"), false, BreakEndVariantType.SupportedType.INTRA_CHR_STRAND_SWITCH_33)
            );
            final List<VariantContext> variants = Arrays.asList(
                    makeBND(new SimpleInterval("20:200201-200201"), new SimpleInterval("20:200606-200606"), Allele.create("A"), "", INV33, true, false, false)
                            .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 105).attribute(MAX_ALIGN_LENGTH, 105).attribute(CONTIG_NAMES, contigName).make(),
                    makeBND(new SimpleInterval("20:200201-200201"), new SimpleInterval("20:200606-200606"), Allele.create("G"), "", INV33, false, false, false)
                            .attribute(TOTAL_MAPPINGS, 1).attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(ALIGN_LENGTHS, 105).attribute(MAX_ALIGN_LENGTH, 105).attribute(CONTIG_NAMES, contigName).make()
            );
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, breakEnds, variants));

            outputStream.reset();
            outputStream.write(rightRightMinus);outputStream.write(leftHomology);outputStream.write(leftRightPlus);
            contigSeq = outputStream.toByteArray();

            contigName = "simple_inv_33_with_hom_-";
            region1 = new AlignmentInterval(new SimpleInterval("20", 200601, 200705), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 200201, 200305), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            simpleChimera = new SimpleChimera(contigName, region1, region2, StrandSwitch.REVERSE_TO_FORWARD, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
            result.add(new TestDataForInversion(region1, region2, contigName, contigSeq, simpleChimera, breakpoints, breakEnds, variants));
        }
        return result;
    }
}
