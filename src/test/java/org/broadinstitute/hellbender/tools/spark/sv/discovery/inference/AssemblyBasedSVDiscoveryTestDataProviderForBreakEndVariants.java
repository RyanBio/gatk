package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;

public class AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants {

    // TODO: 5/23/18 use to fill up AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants
    @Test(groups = "sv")
    public void testRefOrderSwitch() {
        AlignmentInterval region1 = new AlignmentInterval(
                // assigned from chr18 to chr21 to use the dict
                new SimpleInterval("chr21", 39477098, 39477363),
                1 ,268,
                TextCigarCodec.decode("236M2I30M108S"), true, 32, 25, 133, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(
                new SimpleInterval("chr21", 39192594, 39192692),
                252 ,350,
                TextCigarCodec.decode("251S99M26S"), true, 32, 1, 94, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "testContig", NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(simpleChimera,
                "TTCCTTAAAATGCAGGTGAATACAAGAATTAGGTTTCAGGTTTTATATATATATTCTGATATATATATATAATATAACCTGAGATATATATATAAATATATATATTAATATATATTAATATATATAAATATATATATATTAATATATATTTATATATAAATATATATATATTAATATATATAAATATATATAAATATATATATATTAATATATATTAATATATAAATATATATATATTAATATATATTAATATATATAAATATATATATTAATATATATAAATATATATATAAATATATATAAATATATAAATATATATATAAATATATATAAATATATATAAATATATATACACACATACATACACATATACATT".getBytes(),
                TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("chr21", 39192594, 39192594));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("chr21", 39477346, 39477346));
        Assert.assertEquals(breakpoints.getComplication().getHomologyForwardStrandRep(), "ATATATAAATATATATA");
        Assert.assertTrue(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().isEmpty());
    }

    public static final class TestDataBreakEndVariants extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera {

        private TestDataBreakEndVariants(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                         final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                         final SimpleChimera manuallyCuratedSimpleChimera,
                                         final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble,
                                         final List<SvType> manuallyCuratedSVTypes,
                                         final List<VariantContext> manuallyCuratedVariants,
                                         final Class<? extends BreakpointsInference> inferencer) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, manuallyCuratedSimpleChimera, manuallyCuratedBiPathBubble, manuallyCuratedSVTypes, manuallyCuratedVariants, inferencer);
        }

        @Override
        public SAMSequenceDictionary getAppropriateDictionary() {
            return TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return inferencer;
        }
    }

    public static final boolean testDataInitialized;

//    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwap_plus;
//    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwap_minus;
//
//    public static final TestDataBreakEndVariants forIntraChromosomeStrandSwitch_plus;
//    public static final TestDataBreakEndVariants forIntraChromosomeStrandSwitch_minus;
//
//    public static final TestDataBreakEndVariants forInterChromosomeSimple_plus;
//    public static final TestDataBreakEndVariants forInterChromosomeSimple_minus;
//
//    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch_plus;
//    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch_minus;

    static {
        testDataInitialized = true;
    }

    public static List<TestDataBreakEndVariants> getAllTestData() {
        return Collections.emptyList();
    }

    public static List<Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants>> getAllTestDataPaired() {
        return Collections.emptyList();
    }

//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwap_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwap_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeStrandSwitch_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeStrandSwitch_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeSimple_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeSimple_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch_minus() {
//
//    }
}
