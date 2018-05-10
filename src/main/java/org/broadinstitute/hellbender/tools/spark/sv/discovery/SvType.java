package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.EnumUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.*;


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
}
