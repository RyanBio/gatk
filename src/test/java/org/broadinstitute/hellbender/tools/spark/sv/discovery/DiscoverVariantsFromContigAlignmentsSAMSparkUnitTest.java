package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.inferSimpleTypeFromNovelAdjacency;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.SupportedType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.SupportedType.DEL;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.SupportedType.DUP;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.*;

public class DiscoverVariantsFromContigAlignmentsSAMSparkUnitTest extends GATKBaseTest {

    //  for methods used in this CLI are used only in this class

    @Test(groups = "sv")
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region1, region2, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region2, region1, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test(groups = "sv")
    public void testFilterByNextAlignmentMayBeInsertion() {
        final AlignmentInterval overlappingRegion1 = new AlignmentInterval(new SimpleInterval("19", 48699881, 48700034), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval overlappingRegion2 = new AlignmentInterval(new SimpleInterval("19", 48700584, 48700668), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertTrue(DiscoverVariantsFromContigAlignmentsSAMSpark.nextAlignmentMayBeInsertion(overlappingRegion1, overlappingRegion2,  CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, 50,true));
    }

    // TODO: 5/23/18 test
    @Test(groups = "sv")
    public void testParseOneContig() {

    }

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.testDataInitialized) {
            new AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV();
        }
        if(!AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.testDataInitialized) {
            new AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints();
        }
    }


    @DataProvider
    private Object[][] forInferSimpleTypeFromNovelAdjacency() {
        final List<Object[]> data = new ArrayList<>(20);

        {
            // inversion
            data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.manuallyCuratedBiPathBubble, INV.name(), ImmutableSet.of(INV33)});

            data.add(new Object[]{forSimpleInversionWithHom_leftPlus.manuallyCuratedBiPathBubble, INV.name(), ImmutableSet.of(INV55)});

            // simple deletion
            data.add(new Object[]{forSimpleDeletion_plus.manuallyCuratedBiPathBubble, DEL.name(), Collections.emptySet()});

            // simple insertion
            data.add(new Object[]{forSimpleInsertion_minus.manuallyCuratedBiPathBubble, INS.name(), Collections.emptySet()});

            // long range substitution
            data.add(new Object[]{forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble, DEL.name(), Collections.emptySet()});

            // simple deletion with homology
            data.add(new Object[]{forDeletionWithHomology_minus.manuallyCuratedBiPathBubble, DEL.name(), Collections.emptySet()});

            // simple tandem dup contraction from 2 units to 1 unit
            data.add(new Object[]{forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

            // simple tandem dup expansion from 1 unit to 2 units
            data.add(new Object[]{forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

            // simple tandem dup expansion from 1 unit to 2 units and novel insertion
            data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

            // tandem dup expansion from 1 unit to 2 units with pseudo-homology
            data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

            // tandem dup contraction from 2 units to 1 unit with pseudo-homology
            data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

            // tandem dup contraction from 3 units to 2 units
            data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

            // tandem dup expansion from 2 units to 3 units
            data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});
        }

        {
            // inversion
            data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.manuallyCuratedBiPathBubble,
                    inferSimpleTypeFromNovelAdjacency(forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_INV, 14644, INV33});
            data.add(new Object[]{forSimpleInversionWithHom_leftPlus.manuallyCuratedBiPathBubble,
                    inferSimpleTypeFromNovelAdjacency(forSimpleInversionWithHom_leftPlus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_INV, 405, INV55});

            // simple deletion
            data.add(new Object[]{forSimpleDeletion_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleDeletion_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.SupportedType.DEL.name()});

            // simple insertion
            data.add(new Object[]{forSimpleInsertion_minus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleInsertion_minus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_INS, 50, SimpleSVType.SupportedType.INS.name()});

            // long range substitution (i.e. scarred deletion)
            data.add(new Object[]{forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -60, SimpleSVType.SupportedType.DEL.name()});

            // simple deletion with homology
            data.add(new Object[]{forDeletionWithHomology_minus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forDeletionWithHomology_minus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -38, SimpleSVType.SupportedType.DEL.name()});

            // simple tandem dup contraction from 2 units to 1 unit
            data.add(new Object[]{forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -10,
                    DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

            // simple tandem dup expansion from 1 unit to 2 units
            data.add(new Object[]{forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DUP, 10,
                    DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

            // simple tandem dup expansion from 1 unit to 2 units and novel insertion
            data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DUP, 99,
                    DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

            // tandem dup expansion from 1 unit to 2 units with pseudo-homology
            data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DUP, 96,
                    DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

            // tandem dup contraction from 2 units to 1 unit with pseudo-homology
            data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -96,
                    DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});


            // tandem dup contraction from 3 units to 2 units
            data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DEL, -96,
                    DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

            // tandem dup expansion from 2 units to 3 units
            data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble),
                    SYMB_ALT_ALLELE_DUP, 96,
                    DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});
        }
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forInferSimpleTypeFromNovelAdjacency")
    public void testInferSimpleTypeFromNovelAdjacency(final NovelAdjacencyAndAltHaplotype biPathBubble,
                                                      final String expectedTypeString,
                                                      final Set<String> expectedAttributeIDs,
                                                      final String expectedSymbolicAltAlleleStringWithoutBracket,
                                                      final int expectedSvLen,
                                                      final String expectedTypeInfoInIdString) throws IOException {

        final SvType variant = DiscoverVariantsFromContigAlignmentsSAMSpark.inferSimpleTypeFromNovelAdjacency(biPathBubble);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> attributeIDs = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(attributeIDs, expectedAttributeIDs);

        final List<Allele> producedAlleles = AnnotatedVariantProducer.produceAlleles(biPathBubble.getLeftJustifiedLeftRefLoc(), TestUtilsForAssemblyBasedSVDiscovery.b37_reference, variant);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference() && producedAlleles.get(1).isSymbolic());
        Assert.assertEquals(producedAlleles.get(1).toString(), createBracketedSymbAlleleString(expectedSymbolicAltAlleleStringWithoutBracket));

        Assert.assertEquals(variant.getSVLength(), expectedSvLen);

        final String variantId = variant.getInternalVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(INTERVAL_VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedTypeInfoInIdString);
        final String expectedRefContigNameInIdString = biPathBubble.getLeftJustifiedLeftRefLoc().getContig(),
                expectedPOSInfoInIdString = String.valueOf(biPathBubble.getLeftJustifiedLeftRefLoc().getEnd()),
                expectedENDInfoInIdString = String.valueOf(biPathBubble.getLeftJustifiedRightRefLoc().getStart());
        Assert.assertEquals(fields[1], expectedRefContigNameInIdString);
        Assert.assertEquals(fields[2], expectedPOSInfoInIdString);
        Assert.assertEquals(fields[3], expectedENDInfoInIdString);
    }
    @Test(groups = "sv", dataProvider = "forAltAlleleSvLenAndIdProductions")
    public static void testAltAlleleSvLenAndIdProductions(final NovelAdjacencyAndAltHaplotype novelAdjacencyReferenceLocations,
                                                          final SvType simpleType,
                                                          final String expectedSymbolicAltAlleleStringWithoutBracket,
                                                          final int expectedSvLen,
                                                          final String expectedTypeInfoInIdString) throws IOException {

        final List<Allele> producedAlleles = AnnotatedVariantProducer.produceAlleles(novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc(), TestUtilsForAssemblyBasedSVDiscovery.b37_reference, simpleType);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference() && producedAlleles.get(1).isSymbolic());
        Assert.assertEquals(producedAlleles.get(1).toString(), createBracketedSymbAlleleString(expectedSymbolicAltAlleleStringWithoutBracket));

        Assert.assertEquals(simpleType.getSVLength(), expectedSvLen);

        final String variantId = simpleType.getInternalVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(INTERVAL_VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedTypeInfoInIdString);
        final String expectedRefContigNameInIdString = novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc().getContig(),
                expectedPOSInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc().getEnd()),
                expectedENDInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.getLeftJustifiedRightRefLoc().getStart());
        Assert.assertEquals(fields[1], expectedRefContigNameInIdString);
        Assert.assertEquals(fields[2], expectedPOSInfoInIdString);
        Assert.assertEquals(fields[3], expectedENDInfoInIdString);
    }
}
