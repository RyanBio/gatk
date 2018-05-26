package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.BeforeClass;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera;

abstract class AssemblyBasedSVDiscoveryBaseTest extends GATKBaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.testDataInitialized) {
            new AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV();
        }
        if(!AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants.testDataInitialized) {
            new AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants();
        }
        if(!AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.testDataInitialized) {
            new AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants();
        }
    }

    protected List<Tuple2<? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera, ? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> baseDataProviderInPairs() {
        final List<Tuple2<? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera, ? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> data = new ArrayList<>(50);

        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.getAllTestDataPaired() );
        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.getAllTestDataPaired() );
        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants.getAllTestDataPaired() );

        return data;
    }

    protected List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> baseDataProvider() {
        final List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> data = new ArrayList<>(100);

        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV.getAllTestData() );
        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.getAllTestData() );
        data.addAll( AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants.getAllTestData() );

        return data;
    }
}
