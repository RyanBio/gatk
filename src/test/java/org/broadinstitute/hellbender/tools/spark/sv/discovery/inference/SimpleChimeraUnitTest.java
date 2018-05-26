package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.TestDataForSimpleSV;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString;

public class SimpleChimeraUnitTest extends AssemblyBasedSVDiscoveryBaseTest {

    @DataProvider
    private Object[][] forSplitPairStrongEnoughEvidenceForCA() {
        final List<Object[]> data = new ArrayList<>(20);

        String alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t60\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        String alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t60\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";

        AlignmentInterval alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        AlignmentInterval alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, true});

        alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t19\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t60\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";
        alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOneSAMString = "asm018495:tig00014\t16\tchr20\t29851171\t60\t72S582M\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGATTGAGCAGTTTTGAATCACTCTTTATGTAGAATCTGCAAGTGGATATTTGGAGCGTTTTGAGACCTACCGTGGAAAAGCAATTATCTTCAGATAAAAACTACACAGAAGCATTCTGAGAAACTGCTTTATGATGTGTGCATTCATCTCACAGAGTTGAACCTTTCTTTTGATTGAGCAGCTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTCGTGTGCTTTGAGTCCTACCATGGAAAAGGAAATATCTTCACATAAAAAATACTCAGAGGAATTCTGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAATTGAACCTATCTTTAGATTGAGCAGTTTAGAATCTCTCTTTTTGCAGTATCTGCAAGTGGATATTTGGAGCCCTTTGCAGCCTGTGGTGGAAAAGAAAATATCTTCACACAAAAACTACTCGGAAGCATTCTTAGAAACTTCTTTTTGATGTGTGCATTCAAATCACAGAGTTGAACCTATATTTTCAT\t*\tSA:Z:chr20,28843254,-,141M513S,60,15;\tMD:Z:5T15C6A553\tRG:Z:GATKSVContigAlignments\tNM:i:3\tAS:i:567\tXS:i:282";
        alignmentTwoSAMString = "asm018495:tig00014\t2064\tchr20\t28843254\t19\t141M513H\t*\t0\t0\tATCTGTAAGTGGATATTTGGAGCCCTTTGTGGCCTATGGTGGAAAAGGAAATATCTTCAAATAAAAAGTATGCAGAAGCATTCTGAGAAACTTTTTTGTGCTGTGTGCATTCATCTCACAGGGTTGAACCTATCTTATGAT\t*\tSA:Z:chr20,29851171,-,72S582M,60,3;\tMD:Z:5C7A0G52C2C0A2A2A11G2A0C6A4G14A13A6\tRG:Z:GATKSVContigAlignments\tNM:i:15\tAS:i:66\tXS:i:0";
        alignmentOne = fromSAMRecordString(alignmentOneSAMString, true);
        alignmentTwo = fromSAMRecordString(alignmentTwoSAMString, true);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOne = new AlignmentInterval(new SimpleInterval("chr1:10000-10028"), 1, 29, TextCigarCodec.decode("29M"), true, 60, 0, 29, ContigAlignmentsModifier.AlnModType.NONE);
        alignmentTwo = new AlignmentInterval(new SimpleInterval("chr1:10201-10501"), 30, 330, TextCigarCodec.decode("301M"), true, 60, 6, 295, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});

        alignmentOne = new AlignmentInterval(new SimpleInterval("chr1:10201-10501"), 30, 330, TextCigarCodec.decode("301M"), false, 60, 6, 295, ContigAlignmentsModifier.AlnModType.NONE);
        alignmentTwo = new AlignmentInterval(new SimpleInterval("chr1:10000-10028"), 1, 29, TextCigarCodec.decode("29M"), false, 60, 0, 29, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{alignmentOne, alignmentTwo, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_MQ, SimpleNovelAdjacencyInterpreter.MORE_RELAXED_ALIGNMENT_MIN_LENGTH, false});


        return data.toArray(new Object[data.size()][]);
    }
    @Test(dataProvider = "forSplitPairStrongEnoughEvidenceForCA", groups = "sv")
    public void testSplitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                       final AlignmentInterval intervalTwo,
                                                       final int mapQThresholdInclusive,
                                                       final int alignmentLengthThresholdInclusive,
                                                       final boolean expected) {
        Assert.assertEquals(SimpleChimera.splitPairStrongEnoughEvidenceForCA(intervalOne, intervalTwo, mapQThresholdInclusive, alignmentLengthThresholdInclusive),
                expected);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testTypeInference_expectException() {
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100100), 1 ,100, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100101, 100200), 101 ,200, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "1",
                AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, b37_seqDict);
        simpleChimera.inferType(null);
    }


    private static final class TestData {
        private final AlignmentInterval one;
        private final AlignmentInterval two;
        private final SAMSequenceDictionary refDict;

        private final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances;
        private final StrandSwitch expectedStrandSwitch;
        private final boolean expectedIsForwardStrandRepresentation;
        private final Tuple2<SimpleInterval, SimpleInterval> expectedCoordinateSortedRefSpans;

        private final boolean isSimpleTranslocation;
        private final boolean isInversionOrInversionBreakpoint;

        private final TypeInferredFromSimpleChimera typeInferred;

        private TestData(final AlignmentInterval one, final AlignmentInterval two, final SAMSequenceDictionary refDict,
                         final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances,
                         final StrandSwitch expectedStrandSwitch, final boolean expectedIsForwardStrandRepresentation,
                         final Tuple2<SimpleInterval, SimpleInterval> expectedCoordinateSortedRefSpans,
                         final boolean isSimpleTranslocation, final boolean isInversionOrInversionBreakpoint,
                         final TypeInferredFromSimpleChimera typeInferred) {
            this.one = one;
            this.two = two;
            this.refDict = refDict;
            this.expectedStrandSwitch = expectedStrandSwitch;
            this.expectedIsForwardStrandRepresentation = expectedIsForwardStrandRepresentation;
            this.expectedCoordinateSortedRefSpans = expectedCoordinateSortedRefSpans;
            this.expectedDistances = expectedDistances;
            this.isSimpleTranslocation = isSimpleTranslocation;
            this.isInversionOrInversionBreakpoint = isInversionOrInversionBreakpoint;
            this.typeInferred = typeInferred;
        }
    }

    private static List<TestData> casesForSimpleSymbolicVariants() {
        final List<TestData> result = new ArrayList<>(20);

        // simple deletion
        TestDataForSimpleSV testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleDeletion_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleDeletion_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // simple insertion
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleInsertion_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleInsertion_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // long range substitution
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forLongRangeSubstitution_fudgedDel_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forLongRangeSubstitution_fudgedDel_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // simple deletion with homology
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forDeletionWithHomology_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forDeletionWithHomology_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // tandem duplication simple contraction
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupContraction_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupContraction_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // tandem duplication simple expansion
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansion_ins_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansion_ins_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // tandem duplication simple expansion with novel insertion
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansionWithNovelIns_dup_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansionWithNovelIns_dup_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));


        // complex small duplication: expansion from 1 unit to 2 units with pseudo-homology
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_1to2_pseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_1to2_pseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));


        // complex small duplication: contraction from 2 units to 1 unit with pseudo-homology
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to1_pseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to1_pseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));


        // complex small duplication: contraction from 3 units to 2 units without pseudo-homology
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_3to2_noPseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_3to2_noPseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));


        // complex small duplication: expansion from 2 units to 3 units without pseudo-homology
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to3_noPseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to3_noPseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        // complex small duplication: duplicated sequence not long enough to be called DUP (with pseudo homology)
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_1to2_short_pseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_1to2_short_pseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));


        // complex small duplication: duplicated sequence not long enough to be called DUP (without pseudo homology)
        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to3_short_noPseudoHom_plus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, true,
                new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to3_short_noPseudoHom_minus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, testData.getAppropriateDictionary(),
                testData.distances, StrandSwitch.NO_SWITCH, false,
                new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, false, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        return result;
    }

    private static List<TestData> casesForBreakEndVariants() {

        final List<TestData> result = new ArrayList<>(20);

//        // same-chr translocation suspect, forward and reverse representation
//        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 61015129, 61015272), 1, 144, TextCigarCodec.decode("144M148H"), true, 60, 1, 139, ContigAlignmentsModifier.AlnModType.NONE);
//        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 60992732, 60992880), 144, 292, TextCigarCodec.decode("143S149M"), true, 60, 0, 149, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 28861368, 28861775), 1, 409, TextCigarCodec.decode("387M1I21M623H"), false, 60, 22, 286, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28896473, 28897229), 276, 1032, TextCigarCodec.decode("275S757M"), false, 60, 1, 752, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // diff-chr translocation suspect without SS
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 24923683, 24923715), 1, 33, TextCigarCodec.decode("33M130H"), true, 60, 0, 33, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 11590055, 11590197), 21, 163, TextCigarCodec.decode("20S143M"), true, 60, 3, 128, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // diff-chr translocation suspect with SS
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 5374092, 5374747), 1, 656, TextCigarCodec.decode("656M322S"), true, 60, 14, 586, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28764673, 28765145), 506, 978, TextCigarCodec.decode("473M505H"), false, 60, 16, 393, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // same-chr reference order switch, but overlaps (hence incomplete picture)
//        intervalOne = new AlignmentInterval(new SimpleInterval("20", 283, 651), 383, 751, TextCigarCodec.decode("382H369M274H"), true, 60, 23, 254, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("20", 1, 413), 613, 1025, TextCigarCodec.decode("612H413M"), true, 60, 0, 413, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, b37_seqDict, testData.distances, true));
//        // same-chr translocation suspect, forward and reverse representation
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, new Tuple2<>(tuple3s.get(i)._1().referenceSpan, tuple3s.get(i)._2().referenceSpan)}); ++i;
//
//        // diff-chr translocation suspect without SS
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        // diff-chr translocation suspect with SS
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        // same-chr reference order switch, but overlaps (hence incomplete picture)
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
        return result;
    }

    private static List<TestData> casesForInversion() {

        final List<TestData> result = new ArrayList<>(20);

        AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.TestDataForInversion testData = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, b37_seqDict, null,
                StrandSwitch.REVERSE_TO_FORWARD, true, new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, true, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_leftPlus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, b37_seqDict, null,
                StrandSwitch.FORWARD_TO_REVERSE, true, new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, true, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_leftMinus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, b37_seqDict, null,
                StrandSwitch.FORWARD_TO_REVERSE, false, new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, true, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_rightPlus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, b37_seqDict, null,
                StrandSwitch.REVERSE_TO_FORWARD, true, new Tuple2<>(testData.firstAlignment.referenceSpan, testData.secondAlignment.referenceSpan),
                false, true, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        testData = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_rightMinus;
        result.add(new TestData(testData.firstAlignment, testData.secondAlignment, b37_seqDict, null,
                StrandSwitch.REVERSE_TO_FORWARD, false, new Tuple2<>(testData.secondAlignment.referenceSpan, testData.firstAlignment.referenceSpan),
                false, true, testData.manuallyCuratedBiPathBubble.getTypeInferredFromSimpleChimera()));

        return result;
    }

    private static List<TestData> casesForInvertedDuplication() {

        final List<TestData> result = new ArrayList<>(20);

        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr21:25625477-25625587"), 1, 111, TextCigarCodec.decode("111M212H"), false, 60, 0, 111, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr21:25625379-25625595"), 107, 323, TextCigarCodec.decode("106S217M"), true, 60, 0, 127, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, null,
                StrandSwitch.FORWARD_TO_REVERSE, true, new Tuple2<>(intervalTwo.referenceSpan, intervalOne.referenceSpan),
                false, true, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33));

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 48513458, 48513545), 1, 88, TextCigarCodec.decode("88M227H"), true, 39, 1, 83, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 48513297, 48513578), 84, 365, TextCigarCodec.decode("83S282M"), false, 60, 0, 282, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new TestData(intervalOne, intervalTwo, TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict, null,
                StrandSwitch.FORWARD_TO_REVERSE, true, new Tuple2<>(intervalTwo.referenceSpan, intervalOne.referenceSpan),
                false, true, TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_55));
        return result;
    }

    @DataProvider(name = "testRepresentationAndSerialization")
    private Object[][] testRepresentationAndSerialization() {
        final List<Object[]> data = new ArrayList<>(50);

        casesForSimpleSymbolicVariants().forEach(obj -> data.add(new Object[]{obj}));
        casesForBreakEndVariants().forEach(obj -> data.add(new Object[]{obj}));
        casesForInversion().forEach(obj -> data.add(new Object[]{obj}));
        casesForInvertedDuplication().forEach(obj -> data.add(new Object[]{obj}));

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "testRepresentationAndSerialization", groups = "sv")
    public void testRepresentationAndSerialization(final TestData testData) {

        final AlignmentInterval alignmentOne = testData.one;
        final AlignmentInterval alignmentTwo = testData.two;
        final SAMSequenceDictionary refDict = testData.refDict;
        final StrandSwitch expectedStrandSwitch = testData.expectedStrandSwitch;
        final boolean expectedIsForwardStrandRepresentation = testData.expectedIsForwardStrandRepresentation;
        final Tuple2<SimpleInterval, SimpleInterval> expectedCoordinateSortedRefSpans = testData.expectedCoordinateSortedRefSpans;
        final SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead expectedDistances = testData.expectedDistances;
        final boolean expectedIsSimpleTranslocation = testData.isSimpleTranslocation;

        Assert.assertEquals(SimpleChimera.determineStrandSwitch(alignmentOne, alignmentTwo), expectedStrandSwitch);
        Assert.assertEquals(SimpleChimera.isForwardStrandRepresentation(alignmentOne, alignmentTwo, expectedStrandSwitch, refDict), expectedIsForwardStrandRepresentation);

        final SimpleChimera simpleChimera = new SimpleChimera(alignmentOne, alignmentTwo, Collections.emptyList(),
                "dummyName", NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, refDict);
        Assert.assertEquals(simpleChimera.getCoordinateSortedRefSpans(refDict), expectedCoordinateSortedRefSpans);
        if (testData.expectedDistances != null) {
            Assert.assertEquals(simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead(),
                                expectedDistances);
        }
        Assert.assertEquals(simpleChimera.isCandidateSimpleTranslocation(), expectedIsSimpleTranslocation);
        Assert.assertEquals(testData.typeInferred, simpleChimera.inferType(testData.refDict));

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, simpleChimera);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final SimpleChimera roundTrip = (SimpleChimera) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, simpleChimera);
        Assert.assertEquals(roundTrip.hashCode(), simpleChimera.hashCode());
    }
}