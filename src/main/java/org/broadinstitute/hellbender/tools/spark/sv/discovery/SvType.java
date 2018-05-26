package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.EnumUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR;


/**
 * Various types of structural variations
 */
public abstract class SvType {

    protected final String variantId;
    protected final Allele altAllele;
    protected final int svLen;
    protected final Map<String, String> extraAttributes;

    protected static final Map<String, String> noExtraAttributes = Collections.emptyMap();

    protected SvType(final String id, final Allele altAllele, final int len, final Map<String, String> typeSpecificExtraAttributes) {
        variantId = id;
        this.altAllele = altAllele;
        svLen = len;
        extraAttributes = typeSpecificExtraAttributes;
    }

    public final String getInternalVariantId() {
        return variantId;
    }
    public final Allele getAltAllele() {
        return altAllele;
    }
    public final int getSVLength() {
        return svLen;
    }
    public final Map<String, String> getTypeSpecificAttributes() {
        return extraAttributes;
    }

    // TODO: 5/23/18 any better way to do this?
    public static Set<String> getKnownTypes() {
        final SortedSet<String> knownTypes = new TreeSet<>( EnumUtils.getEnumMap(SimpleSVType.SupportedType.class).keySet() );

        knownTypes.add(GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR);

        for (final BreakEndVariantType.SupportedType supportedType : BreakEndVariantType.SupportedType.values()) {
            knownTypes.add(supportedType.name());
        }

        return Collections.unmodifiableSortedSet(knownTypes);
    }

    public static String makeLocationPartOfID(final String chr1, final int pos1, final String chr2, final int pos2) {
        return chr1 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + pos1 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + (chr2.equals(chr1) ? "" : chr2 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR)
                + pos2;
    }

    public static String makeLocationPartOfID(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
        String leftContig = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig();
        String rightContig = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getContig();
        int pos1 = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart();
        int pos2 = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd();
        return makeLocationPartOfID(leftContig, pos1, rightContig, pos2);
    }
}
