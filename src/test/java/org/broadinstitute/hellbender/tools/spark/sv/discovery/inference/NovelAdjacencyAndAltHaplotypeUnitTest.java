package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.collect.ImmutableSet;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.SupportedType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.SYMB_ALT_ALLELE_INS;


public class NovelAdjacencyAndAltHaplotypeUnitTest extends AssemblyBasedSVDiscoveryBaseTest {

    @DataProvider(name = "forKryoSerializationAndHashCode")
    private Object[][] forKryoSerializationAndHashCode() {
        final List<Object[]> data = new ArrayList<>();
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testData : baseDataProvider()) {
            data.add(new Object[]{testData.manuallyCuratedBiPathBubble});
        }
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forKryoSerializationAndHashCode")
    public void testKryoSerializerAndHashCode(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) throws IOException {
        try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {

            final Output out = new Output(bos);
            final Kryo kryo = new Kryo();
            kryo.writeClassAndObject(out, novelAdjacencyAndAltHaplotype);
            out.flush();

            try ( final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray()) ) {
                final Input in = new Input(bis);
                @SuppressWarnings("unchecked")
                final NovelAdjacencyAndAltHaplotype roundTrip = (NovelAdjacencyAndAltHaplotype) kryo.readClassAndObject(in);
                Assert.assertEquals(roundTrip, novelAdjacencyAndAltHaplotype);
                Assert.assertEquals(roundTrip.hashCode(), novelAdjacencyAndAltHaplotype.hashCode());
            }
        }
    }

    // TODO: 5/22/18 add test coverage for BND cases and Inversion cases
    @DataProvider(name = "forToSimpleOrBNDTypes")
    private Object[][] forToSimpleOrBNDTypes() {
        final List<Object[]> data = new ArrayList<>(20);

        final Set<String> defaultKeys = Collections.emptySet();

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), defaultKeys) )});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});

        // long range substitution fudged del
        data.add(new Object[]{forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), defaultKeys) )});

        // long range substitution fat insertion
        data.add(new Object[]{forLongRangeSubstitution_fatIns_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});

        // long range substitution del+ins
        data.add(new Object[]{forLongRangeSubstitution_DelAndIns_plus.manuallyCuratedBiPathBubble,
                Arrays.asList( new Tuple2<>(DEL.name(), defaultKeys),
                        new Tuple2<>(INS.name(), defaultKeys))});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), defaultKeys) )});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)) )});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});
        data.add(new Object[]{forSimpleTanDupExpansion_dup_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)) )});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_ins_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)) )});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)) )});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)) )});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)) )});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)) )});

        // short tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_short_pseudoHom_plus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});
        // short tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_short_noPseudoHom_minus.manuallyCuratedBiPathBubble,
                Collections.singletonList( new Tuple2<>(INS.name(), defaultKeys) )});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forToSimpleOrBNDTypes")
    public void testToSimpleOrBNDTypes(final NovelAdjacencyAndAltHaplotype breakpoints,
                                       final List<Tuple2<String, Set<String>>> expectedTypeStringAndAttributeKeys) {

        final List<SvType> variants = breakpoints.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict);
        Assert.assertEquals(variants.size(), expectedTypeStringAndAttributeKeys.size());
        for (int i = 0; i < variants.size(); ++i) {
            final SvType variant = variants.get(i);
            Assert.assertEquals(variant.toString(), expectedTypeStringAndAttributeKeys.get(i)._1);
            final Set<String> attributeIDs = variant.getTypeSpecificAttributes().keySet();
            Assert.assertEquals(attributeIDs, expectedTypeStringAndAttributeKeys.get(i)._2);
        }
    }

    @Test(groups = "sv", dataProvider = "forAltAlleleSvLenAndIdProductions_SDFLACAS")
    public static void testAltAlleleSvLenAndIdProductions_SDFLACAS(final NovelAdjacencyAndAltHaplotype novelAdjacencyReferenceLocations,
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

    @DataProvider(name = "forAltAlleleSvLenAndIdProductions_SDFLACAS")
    private Object[][] forAltAlleleSvLenAndIdProductions_SDFLACAS() {
        final List<Object[]> data = new ArrayList<>(20);

        // no inversion case because new code path doesn't call inversion (instead, BND)

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus.manuallyCuratedBiPathBubble, forSimpleDeletion_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.SupportedType.DEL.name()});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus.manuallyCuratedBiPathBubble, forSimpleInsertion_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 50, SimpleSVType.SupportedType.INS.name()});

        // long range substitution fudged del
        data.add(new Object[]{forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble, forLongRangeSubstitution_fudgedDel_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -60, SimpleSVType.SupportedType.DEL.name()});

        // long range substitution fat ins
        data.add(new Object[]{forLongRangeSubstitution_fatIns_minus.manuallyCuratedBiPathBubble, forLongRangeSubstitution_fatIns_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 60, SimpleSVType.SupportedType.INS.name()});

        // long range substitution linked insertion and deletion
        List<SvType> svTypes = forLongRangeSubstitution_DelAndIns_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict);
        data.add(new Object[]{forLongRangeSubstitution_DelAndIns_plus.manuallyCuratedBiPathBubble, svTypes.get(0),
                SYMB_ALT_ALLELE_DEL, -60, SimpleSVType.SupportedType.DEL.name()});
        data.add(new Object[]{forLongRangeSubstitution_DelAndIns_plus.manuallyCuratedBiPathBubble, svTypes.get(1),
                SYMB_ALT_ALLELE_INS, 55, SimpleSVType.SupportedType.INS.name()});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus.manuallyCuratedBiPathBubble, forDeletionWithHomology_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -38, SimpleSVType.SupportedType.DEL.name()});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble, forSimpleTanDupContraction_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -10,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units that will be called as insertion
        data.add(new Object[]{forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble, forSimpleTanDupExpansion_ins_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 10,
                SimpleSVType.SupportedType.INS.name()});

        // simple tandem dup expansion from 1 unit to 2 units that will be called as duplication
        data.add(new Object[]{forSimpleTanDupExpansion_dup_minus.manuallyCuratedBiPathBubble, forSimpleTanDupExpansion_dup_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DUP, 55,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion that will be called as insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_ins_plus.manuallyCuratedBiPathBubble, forSimpleTanDupExpansionWithNovelIns_ins_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 326,
                SimpleSVType.SupportedType.INS.name()});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion that will be called as duplication
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble, forSimpleTanDupExpansionWithNovelIns_dup_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DUP, 99,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble, forComplexTanDup_1to2_pseudoHom_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble, forComplexTanDup_2to1_pseudoHom_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});


        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble, forComplexTanDup_3to2_noPseudoHom_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble, forComplexTanDup_2to3_noPseudoHom_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // short tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_short_pseudoHom_plus.manuallyCuratedBiPathBubble, forComplexTanDup_1to2_short_pseudoHom_plus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 42,
                SimpleSVType.SupportedType.INS.name()});
        // short tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_short_noPseudoHom_minus.manuallyCuratedBiPathBubble, forComplexTanDup_2to3_short_noPseudoHom_minus.manuallyCuratedBiPathBubble.toSimpleOrBNDTypes(TestUtilsForAssemblyBasedSVDiscovery.b37_reference, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict).get(0),
                SYMB_ALT_ALLELE_INS, 42,
                SimpleSVType.SupportedType.INS.name()});

        return data.toArray(new Object[data.size()][]);
    }

}