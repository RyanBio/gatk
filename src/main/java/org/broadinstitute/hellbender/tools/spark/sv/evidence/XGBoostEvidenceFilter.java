package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import biz.k11i.xgboost.Predictor;
import biz.k11i.xgboost.learner.ObjFunction;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.AbstractFeatureReader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;

/**
 * A class that acts as a filter for BreakpointEvidence.
 * Features are calculated according to evidence type, overlap information, mapping quality, etc.
 * A trained classifier scores the probability the evidence overlaps a breakpoint interval, and passes evidence above
 * the specified threshold.
 */
public final class XGBoostEvidenceFilter implements Iterator<BreakpointEvidence> {
    private static final boolean USE_FAST_MATH_EXP = true;
    private static final List<String> DEFAULT_EVIDENCE_TYPE_ORDER = Arrays.asList(
            "TemplateSizeAnomaly",
            "MateUnmapped", "InterContigPair",
            "SplitRead", "LargeIndel", "WeirdTemplateSize", "SameStrandPair", "OutiesPair"
    );

    private final PartitionCrossingChecker partitionCrossingChecker;

    private final Predictor predictor;
    private final Map<String, Integer> evidenceTypeMap;
    private final double thresholdProbability;
    private final ReadMetadata readMetadata;

    private final EvidenceOverlapChecker evidenceOverlapChecker;
    private final Map<BreakpointEvidence, UnscaledOverlapInfo> rawFeatureCache;

    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;
    private final SVIntervalTree<Double> genomeGaps;
    private final SVIntervalTree<Double> umapS100Mappability;

    public XGBoostEvidenceFilter(
            final Iterator<BreakpointEvidence> evidenceItr,
            final ReadMetadata readMetadata,
            final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params,
            final PartitionCrossingChecker partitionCrossingChecker
    ) {
        if(params.svGenomeGapsFile == null) {
            throw new IllegalArgumentException("XGBoostEvidenceFilter requires non-null sv-genome-gaps-file");
        }
        if(params.svGenomeCentromeresFile == null) {
            throw new IllegalArgumentException("XGBoostEvidenceFilter requires non-null sv-genome-centromeres-file");
        }
        if(params.svGenomeUmapS100File == null) {
            throw new IllegalArgumentException("XGBoostEvidenceFilter requires non-null sv-genome-umap-s100-file");
        }
        predictor = loadPredictor(params.svEvidenceFilterModelFile);
        evidenceTypeMap = loadEvidenceTypeMap(params.svCategoricalVariablesFile);
        genomeGaps = new SVIntervalTree<>();
        loadGenomeIntervals(this.genomeGaps, params.svGenomeGapsFile, readMetadata);
        loadGenomeIntervals(this.genomeGaps, params.svGenomeCentromeresFile, readMetadata);
        umapS100Mappability = new SVIntervalTree<>();
        loadGenomeIntervals(this.umapS100Mappability, params.svGenomeUmapS100File, readMetadata);

        this.partitionCrossingChecker = partitionCrossingChecker;
        thresholdProbability = params.svEvidenceFilterThresholdProbability;
        this.readMetadata = readMetadata;

        evidenceOverlapChecker = new EvidenceOverlapChecker(evidenceItr, readMetadata, params.minEvidenceMapQ);
        rawFeatureCache = new HashMap<>();

        listItr = null;
        treeItr = evidenceOverlapChecker.getTreeItr();
    }


    public static Predictor loadPredictor(final String modelFileLocation) {
        ObjFunction.useFastMathExp(USE_FAST_MATH_EXP);
        try(final InputStream inputStream = getInputStream(modelFileLocation)) {
            return new Predictor(inputStream);
        } catch(Exception e) {
            throw new GATKException(
                    "Unable to load predictor from classifier file " + modelFileLocation + ": " + e.getMessage()
            );
        }
    }

    @VisibleForTesting
    static InputStream getInputStream(final String fileLocation) throws IOException {
        if(fileLocation.startsWith("gatk-resources::")) {
            final InputStream inputStream = XGBoostEvidenceFilter.class.getResourceAsStream(fileLocation.substring(16));
            return AbstractFeatureReader.hasBlockCompressedExtension(fileLocation) ?
                    IOUtils.makeZippedInputStream(new BufferedInputStream(inputStream))
                    : inputStream;
        } else {
            return BucketUtils.openFile(fileLocation);
        }
    }

    @VisibleForTesting
    static Map<String, Integer> loadEvidenceTypeMap(final String categoricalVariablesMapFile) {
        final HashMap<String, Integer> evidenceTypeMap = new HashMap<>();
        if(categoricalVariablesMapFile == null || categoricalVariablesMapFile.isEmpty()) {
            // no categorical variables file, just use default evidence order
            for(int index = 0; index < DEFAULT_EVIDENCE_TYPE_ORDER.size(); ++index) {
                evidenceTypeMap.put(DEFAULT_EVIDENCE_TYPE_ORDER.get(index), index);
            }
            return evidenceTypeMap;
        }
        // override default evidence order with saved order
        try(final InputStream inputStream = getInputStream(categoricalVariablesMapFile)) {
            final JsonNode testDataNode = new ObjectMapper().readTree(inputStream);
            final JsonNode evidenceTypeArrayNode = testDataNode.get("evidence_type");
            if(!evidenceTypeArrayNode.isArray()) {
                throw new IllegalArgumentException("evidenceTypeNode does not encode a valid array");
            }
            int index = 0;
            for(final JsonNode evidenceTypeNode : evidenceTypeArrayNode) {
                evidenceTypeMap.put(evidenceTypeNode.asText(), index);
                ++index;
            }
            return evidenceTypeMap;
        } catch(IOException e) {
            throw new GATKException(
                    "Unable to load evidence-type map from file " + categoricalVariablesMapFile + ": " + e.getMessage()
            );
        }
    }

    /**
     * Load file containing genome_intervals onto an SVIntervalTree. Each interval has weight 1.0. If more than one
     * interval occupies the same [begin, end), add their weights rather than having two intervals at the same spot.
     */
    static void loadGenomeIntervals(final SVIntervalTree<Double> genome_intervals, final String genomeIntervalsFile,
                                    final ReadMetadata readMetadata) {
        // Some BAM files (especially in test suite) have unusual headers that don't include the usual primary contig
        // names. In this case, enforce our own contig ordering.
        boolean goodSamHeader;
        try {
            readMetadata.getContigID("chr1");
            goodSamHeader = true;
        } catch(GATKException e) {
            goodSamHeader = false;
        }
        final boolean mapContigIdWithReadMetadata = goodSamHeader;

        int lineNumber = -1;
        try(final InputStream inputStream = getInputStream(genomeIntervalsFile)) {
            final BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            String line;
            loop_lines:
            while((line = reader.readLine()) != null) {
                ++lineNumber;
                if(line.startsWith("#")) {
                    continue; // skip header
                }
                final String[] splitLine = line.split("\t");
                if(splitLine.length < 3) {
                    throw new GATKException("Error reading " + genomeIntervalsFile + " line " + lineNumber
                            + ": Not enough columns to make genomic interval");
                }

                final int contig;
                if(mapContigIdWithReadMetadata) {
                    contig = readMetadata.getContigID(splitLine[0]);
                } else {
                    switch(splitLine[0]) {
                        case "chrX":
                            contig = 22;
                            break;
                        case "chrY":
                            contig = 23;
                            break;
                        default:
                            try {
                                contig = Integer.parseInt(splitLine[0].substring(3));
                            } catch(NumberFormatException e) {
                                continue loop_lines;
                            }
                    }
                }
                final int begin = Integer.parseInt(splitLine[1]) + 1;
                final int end = Integer.parseInt(splitLine[2]) + 1;
                final SVInterval svInterval = new SVInterval(contig, begin, end);
                final SVIntervalTree.Entry<Double> entry = genome_intervals.find(svInterval);
                // each interval has weight 1.0 for overlap. If multiple intervals exactly coincide, just add weight.
                if(entry == null) {
                    genome_intervals.put(svInterval, 1.0);
                } else {
                    entry.setValue(entry.getValue() + 1.0);
                }
            }
        } catch(IOException e) {
            throw new GATKException("Unable to load genome intervals from file \"" + genomeIntervalsFile + "\" line "
                                    + lineNumber + ": " + e.getMessage());
        }
    }


    @Override
    public boolean hasNext() {
        if ( listItr != null && listItr.hasNext() ) {
            return true;
        }
        listItr = null;
        boolean result = false;
        while ( !result && treeItr.hasNext() ) {
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = treeItr.next();
            final SVInterval curInterval = entry.getInterval();
            final List<BreakpointEvidence> evidenceList = entry.getValue();
            if( isValidated(entry.getValue()) || partitionCrossingChecker.onBoundary(curInterval) ) {
                // already validated (no need to mark validated again) or on partition boundary (punt for now)
                result = true;
            } else if( anyPassesFilter(evidenceList) ) {
                evidenceList.forEach(ev -> ev.setValidated(true));
                result = true;
            }

            if ( result ) {
                listItr = entry.getValue().iterator();
            }
        }
        return result;
    }

    @Override
    public BreakpointEvidence next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("No next element.");
        }
        return listItr.next();
    }

    private boolean isValidated( final List<BreakpointEvidence> evList ) {
        for ( final BreakpointEvidence ev : evList ) {
            if ( ev.isValidated() ) return true;
        }
        return false;
    }

    boolean anyPassesFilter(final List<BreakpointEvidence> evidenceList) {
        for(final BreakpointEvidence evidence : evidenceList) {
            final double evidenceGoodProbability = predictor.predictSingle(getFeatures(evidence));
            if(evidenceGoodProbability > thresholdProbability) {
                return true;
            }
        }
        return false;
    }

    /**
     * Compute features vector for a piece of BreakpointEvidence
     */
    @VisibleForTesting
    EvidenceFeatures getFeatures(final BreakpointEvidence evidence) {
        // create new struct for these two, use CigarOperator to update if instanceof ReadEvidence
        final CigarQualityInfo cigarQualityInfo = new CigarQualityInfo(evidence);
        final double evidenceType = evidenceTypeMap.get(evidence.getClass().getSimpleName());
        final double mappingQuality = (double)getMappingQuality(evidence);
        final double templateSize = getTemplateSize(evidence);

        // calculate these similar to BreakpointDensityFilter, but always calculate full totals, never end early.
        final CoverageScaledOverlapInfo individualOverlapInfo = getIndividualOverlapInfo(evidence);
        final CoverageScaledOverlapInfo clusterOverlapInfo = getClusterOverlapInfo(evidence);

        // calculate properties related to overlap of intervals on the reference genome
        final double referenceGapOverlap = getGenomeIntervalsOverlap(evidence, genomeGaps);
        final double umapS100 = getGenomeIntervalsOverlap(evidence, umapS100Mappability);

        return new EvidenceFeatures(
            new double[]{
                cigarQualityInfo.basesMatched, cigarQualityInfo.referenceLength, evidenceType, mappingQuality, templateSize,
                individualOverlapInfo.numOverlap, individualOverlapInfo.overlapMappingQuality,
                individualOverlapInfo.meanOverlapMappingQuality, individualOverlapInfo.numCoherent,
                individualOverlapInfo.coherentMappingQuality,
                clusterOverlapInfo.numOverlap, clusterOverlapInfo.overlapMappingQuality,
                clusterOverlapInfo.meanOverlapMappingQuality, clusterOverlapInfo.numCoherent,
                clusterOverlapInfo.coherentMappingQuality,
                referenceGapOverlap, umapS100
            }
        );
    }

    private int getMappingQuality(final BreakpointEvidence evidence) {
        // Note: return "max" mapping quality for non-ReadEvidence. Reasoning: some features depend on sum or average of
        // read qualities. Non-ReadEvidence isn't *bad* per se, so give it a good score.
        return evidence instanceof BreakpointEvidence.ReadEvidence ?
            ((BreakpointEvidence.ReadEvidence) evidence).getMappingQuality() : 60;
    }

    private double getTemplateSize(final BreakpointEvidence evidence) {
        if(evidence instanceof BreakpointEvidence.ReadEvidence) {
            // For ReadEvidence, return templateSize as percentile of library's cumulative density function:
            final String readGroup = ((BreakpointEvidence.ReadEvidence) evidence).getReadGroup();
            final String library = readMetadata.getReadGroupToLibraryMap().get(readGroup);
            final LibraryStatistics libraryStatistics = readMetadata.getLibraryStatistics(library);
            final IntHistogram.CDF templateSizeCDF = libraryStatistics.getCDF();
            final int cdfBin = Integer.min(Math.abs(((BreakpointEvidence.ReadEvidence) evidence).getTemplateSize()),
                                           templateSizeCDF.size() - 1);
            return templateSizeCDF.getFraction(cdfBin);
        } else if(evidence instanceof BreakpointEvidence.TemplateSizeAnomaly) {
            // For TemplateSizeAnomaly, return readCount scaled by meanGenomeCoverage
            return (double)(((BreakpointEvidence.TemplateSizeAnomaly) evidence).getReadCount()) / readMetadata.getCoverage();
        } else {
            throw new IllegalArgumentException(
                    "templateSize feature only defined for ReadEvidence and TemplateSizeAnomaly"
            );
        }
    }

    /**
     * For given BreakpointEvidence, calculate number of overlappers, coherent pieces of BreakpointEvidence, and
     * the sums of their mapping qualities. Cache this information so that it can be looked up when computing max or
     * average of these values over the set of overlapping evidence (as in getClusterOverlapInfo).
     */
    private void cacheOverlapInfo(final BreakpointEvidence evidence) {
        int numOverlap = 0;
        int overlapMappingQuality = 0;
        int numCoherent = 0;
        int coherentMappingQuality = 0;
        for(final Iterator<Tuple2<BreakpointEvidence, Boolean>> overlapper_itr = evidenceOverlapChecker.overlappersWithCoherence(evidence);
                overlapper_itr.hasNext();) {
            final Tuple2<BreakpointEvidence, Boolean> itr_results = overlapper_itr.next();
            final BreakpointEvidence overlapper = itr_results._1;
            if(overlapper == evidence) {
                continue; // don't count self-overlap
            }
            numOverlap += 1;
            final int mappingQuality = getMappingQuality(overlapper);
            overlapMappingQuality += mappingQuality;

            final boolean isCoherent = itr_results._2;
            if(isCoherent) {
                numCoherent += 1;
                coherentMappingQuality += mappingQuality;
            }
        }
        rawFeatureCache.put(evidence,
                            new UnscaledOverlapInfo(numOverlap, numCoherent, overlapMappingQuality, coherentMappingQuality));
    }

    private CoverageScaledOverlapInfo getIndividualOverlapInfo(final BreakpointEvidence evidence) {
        if(!rawFeatureCache.containsKey(evidence)) {
            cacheOverlapInfo(evidence);
        }
        final UnscaledOverlapInfo evidenceFeatureCache = rawFeatureCache.get(evidence);
        final double meanOverlapMappingQuality = ((double)evidenceFeatureCache.overlapMappingQuality)
                                               / evidenceFeatureCache.numOverlap;
        return new CoverageScaledOverlapInfo(evidenceFeatureCache.numOverlap, evidenceFeatureCache.numCoherent,
                                             evidenceFeatureCache.overlapMappingQuality, evidenceFeatureCache.coherentMappingQuality,
                                             meanOverlapMappingQuality, readMetadata.getCoverage());
    }

    private CoverageScaledOverlapInfo getClusterOverlapInfo(final BreakpointEvidence evidence) {
        int clusterNumOverlap = 0;
        int clusterNumCoherent = 0;
        int clusterOverlapMappingQuality = 0;
        int clusterCoherentMappingQuality = 0;
        double clusterMeanOverlapMappingQuality = 0.0;
        for (final Iterator<BreakpointEvidence> overlapper_itr = evidenceOverlapChecker.overlappers(evidence); overlapper_itr.hasNext(); ) {
            final BreakpointEvidence overlapper = overlapper_itr.next();
            if (overlapper == evidence) {
                continue; // don't count self-overlap in cluster features
            }
            if(!rawFeatureCache.containsKey(overlapper)) {
                cacheOverlapInfo(overlapper);
            }
            final UnscaledOverlapInfo overlapperFeatureCache = rawFeatureCache.get(overlapper);
            clusterNumOverlap = Math.max(clusterNumOverlap, overlapperFeatureCache.numOverlap);
            clusterNumCoherent = Math.max(clusterNumCoherent, overlapperFeatureCache.numCoherent);
            clusterOverlapMappingQuality = Math.max(clusterOverlapMappingQuality, overlapperFeatureCache.overlapMappingQuality);
            clusterCoherentMappingQuality = Math.max(clusterCoherentMappingQuality, overlapperFeatureCache.coherentMappingQuality);

            final double meanOverlapMappingQuality = ((double) overlapperFeatureCache.overlapMappingQuality)
                                                    / overlapperFeatureCache.numOverlap;
            clusterMeanOverlapMappingQuality = Math.max(clusterMeanOverlapMappingQuality, meanOverlapMappingQuality);
        }

        return new CoverageScaledOverlapInfo(clusterNumOverlap, clusterNumCoherent, clusterOverlapMappingQuality,
                clusterCoherentMappingQuality, clusterMeanOverlapMappingQuality, readMetadata.getCoverage());
    }

    private static <T> void addToTree( final SVIntervalTree<List<T>> tree,
                                       final SVInterval interval,
                                       final T value ) {
        final SVIntervalTree.Entry<List<T>> entry = tree.find(interval);
        if ( entry != null ) {
            entry.getValue().add(value);
        } else {
            final List<T> valueList = new ArrayList<>(1);
            valueList.add(value);
            tree.put(interval, valueList);
        }
    }

    private static class EvidenceOverlapChecker {
        private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
        private final ReadMetadata readMetadata;
        private final int minEvidenceMapQ;

        EvidenceOverlapChecker(final Iterator<BreakpointEvidence> evidenceItr, final ReadMetadata readMetadata,
                               final int minEvidenceMapQ) {
            this.readMetadata = readMetadata;
            this.minEvidenceMapQ = minEvidenceMapQ;
            evidenceTree = new SVIntervalTree<>();
            while ( evidenceItr.hasNext() ) {
                final BreakpointEvidence evidence = evidenceItr.next();
                addToTree(evidenceTree, evidence.getLocation(), evidence);
            }
        }

        Iterator<BreakpointEvidence> overlappers(final BreakpointEvidence evidence) {
            return new OverlapperIterator(evidence, evidenceTree);
        }

        // returns iterator to overlapping evidence, paired with a Boolean that is true if the evidence is also coherent
        Iterator<Tuple2<BreakpointEvidence, Boolean>> overlappersWithCoherence(final BreakpointEvidence evidence) {
            return new OverlapAndCoherenceIterator(evidence, evidenceTree, readMetadata, minEvidenceMapQ);
        }

        Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> getTreeItr() {
            return evidenceTree.iterator();
        }

        private static class OverlapperIterator  implements Iterator<BreakpointEvidence> {
            private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
            private Iterator<BreakpointEvidence> listItr;

            OverlapperIterator(final BreakpointEvidence evidence,
                               final SVIntervalTree<List<BreakpointEvidence>> evidenceTree) {
                treeItr = evidenceTree.overlappers(evidence.getLocation());
                listItr = null;
            }

            @Override
            public boolean hasNext() {
                if ( listItr != null && listItr.hasNext() ) {
                    return true;
                }
                while(treeItr.hasNext()) {
                    listItr = treeItr.next().getValue().iterator();
                    if(listItr.hasNext()) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public BreakpointEvidence next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("No next element.");
                }
                return listItr.next();
            }
        }

        /**
         * Implements iterator of BreakpointEvidence that overlaps a specified piece of evidence. Each call to next()
         * returns overlapping BreakpointEvidence paired with a Boolean that is true if the overlapper is also coherent
         * with the evidence.
         */
        private static class OverlapAndCoherenceIterator implements Iterator<Tuple2<BreakpointEvidence, Boolean>> {
            private final ReadMetadata readMetadata;
            private final int minEvidenceMapQ;
            private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
            private Iterator<BreakpointEvidence> listItr;
            private final boolean checkCoherence;
            private final Boolean isForwardStrand;
            private final List<StrandedInterval> distalTargets;

            OverlapAndCoherenceIterator(final BreakpointEvidence evidence,
                                        final SVIntervalTree<List<BreakpointEvidence>> evidenceTree,
                                        final ReadMetadata readMetadata,
                                        final int minEvidenceMapQ) {
                this.readMetadata = readMetadata;
                this.minEvidenceMapQ = minEvidenceMapQ;
                treeItr = evidenceTree.overlappers(evidence.getLocation());
                isForwardStrand = evidence.isEvidenceUpstreamOfBreakpoint();
                checkCoherence = (isForwardStrand != null) && evidence.hasDistalTargets(readMetadata, minEvidenceMapQ);
                if(checkCoherence) {
                    distalTargets = new ArrayList<>(evidence.getDistalTargets(readMetadata, minEvidenceMapQ));
                } else {
                    distalTargets = null;
                }
            }

            @Override
            public boolean hasNext() {
                if ( listItr != null && listItr.hasNext() ) {
                    return true;
                }
                while(treeItr.hasNext()) {
                    listItr = treeItr.next().getValue().iterator();
                    if(listItr.hasNext()) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public Tuple2<BreakpointEvidence, Boolean> next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("No next element.");
                }
                final BreakpointEvidence overlapper = listItr.next();
                Boolean isCoherent = false;
                if(checkCoherence) {
                    Boolean overlapperStrand = overlapper.isEvidenceUpstreamOfBreakpoint();
                    if(overlapperStrand == isForwardStrand && overlapper.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                        for(final StrandedInterval distalTarget : distalTargets) {
                            for(final StrandedInterval overlapperDistalTarget : new ArrayList<>(overlapper.getDistalTargets(readMetadata, minEvidenceMapQ))) {
                                if(distalTarget.getStrand() == overlapperDistalTarget.getStrand()
                                        && distalTarget.getInterval().overlaps(overlapperDistalTarget.getInterval())) {
                                    isCoherent = true; // could break, but everything has 1 distal target, so false optimization
                                }
                            }
                        }
                    }
                }
                return new Tuple2<>(overlapper, isCoherent);
            }
        }
    }

    private static class UnscaledOverlapInfo {
        final int numOverlap;
        final int numCoherent;
        final int overlapMappingQuality;
        final int coherentMappingQuality;
        UnscaledOverlapInfo(final int numOverlap, final int numCoherent, final int overlapMappingQuality,
                            final int coherentMappingQuality) {
            this.numOverlap = numOverlap;
            this.numCoherent = numCoherent;
            this.overlapMappingQuality = overlapMappingQuality;
            this.coherentMappingQuality = coherentMappingQuality;
        }
    }

    /**
     * Class that takes raw overlap info and automatically scales it to meanGenomeCoverage, storing the result for later retrieval
     */
    private static class CoverageScaledOverlapInfo {
        final double numOverlap;
        final double overlapMappingQuality;
        final double meanOverlapMappingQuality;
        final double numCoherent;
        final double coherentMappingQuality;

        CoverageScaledOverlapInfo(final int numOverlap, final int numCoherent, final int overlapMappingQuality,
                                  final int coherentMappingQuality, final double meanOverlapMappingQuality,
                                  final double coverage) {
            this.numOverlap = ((double)numOverlap) / coverage;
            this.overlapMappingQuality = ((double)overlapMappingQuality) / coverage;
            this.numCoherent = ((double)numCoherent) / coverage;
            this.coherentMappingQuality = ((double)coherentMappingQuality) / coverage;
            this.meanOverlapMappingQuality = meanOverlapMappingQuality;
        }
    }

    private static class CigarQualityInfo {
        final int basesMatched;
        final int referenceLength;

        CigarQualityInfo(final BreakpointEvidence evidence) {
            int numMatched = 0;
            int refLength = 0;
            if(evidence instanceof BreakpointEvidence.ReadEvidence) {
                final String cigarString = ((BreakpointEvidence.ReadEvidence) evidence).getCigarString();
                for (final CigarElement element : TextCigarCodec.decode(cigarString).getCigarElements()) {
                    final CigarOperator op = element.getOperator();
                    if (op.consumesReferenceBases()) {
                        refLength += element.getLength();
                        if (op.consumesReadBases()) {
                            numMatched += element.getLength();
                        }
                    }
                }
            }
            basesMatched = numMatched;
            referenceLength = refLength;
        }
    }

    final double getGenomeIntervalsOverlap(final BreakpointEvidence evidence,
                                           final SVIntervalTree<Double> genomeIntervals) {
        final SVInterval svInterval = evidence.getLocation();
        double overlap = 0.0;
        for(final Iterator<SVIntervalTree.Entry<Double>> overlapperItr = genomeIntervals.overlappers(svInterval);
            overlapperItr.hasNext();) {
            final SVIntervalTree.Entry<Double> overlapEntry = overlapperItr.next();
            final SVInterval overlapper = overlapEntry.getInterval();
            final double overlapLength = Double.min(svInterval.getEnd(), overlapper.getEnd())
                    - Double.max(svInterval.getStart(), overlapper.getStart());
            overlap += overlapEntry.getValue() * overlapLength / svInterval.getLength();
        }
        return overlap;
    }
}
