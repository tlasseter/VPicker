package com.Sts.Types.PreStack;

import javax.media.opengl.GL;

import com.Sts.DBTypes.*;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsGLDraw;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsVolumeConstructorTimer;
import com.Sts.Utilities.StsSeismicFilter;
import com.Sts.Utilities.StsToolkit;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.DataCube.StsPreStackLineFileBlocks;
import com.Sts.Types.*;
import com.Sts.IO.*;
import com.Sts.SeismicAttributes.*;

import java.util.*;

/**
 * A superGather consists of a single gather or a group of gathers in a rectangle or cross configuration.
 * Row and col are the indices of this superGather.  Computations are done on a gather-by-gather basis,
 * so gatherRow and gatherCol are the indices of the current gather being computed.  gatherData is the
 * trace data for the current gather at gatherRow & gatherCol.
 */
public class StsGather
{
    /** superGather which owns this gather */
    public StsSuperGather superGather;
    /** volume the superGather and this gather belong to */
    public StsPreStackLineSet lineSet;
    /** row for the current gather being processed; ranges about the superGatherRow */
    public int gatherRow = 0;
    /** col for the current gather being processed ; ranges about the superGatherCol */
    public int gatherCol = 0;
    /** true if gather has been successfully initialized */
    public boolean initialized = false;
    /** line this gather is on */
    public StsPreStackLine line;
    /** indicates this line has been NMOed */
    //    public boolean isNMOed = false;
    /** index of current gather on the line */
    public int nLineGather = -1;
    /** index of first trace in gather isVisible */
    public int nFirstGatherLineTrace;
    /** index of last trace in gather isVisible */
    public int nLastGatherLineTrace;
    /** number of traces in gather */
    public int nGatherTraces;
    /** pointTypes for gather */
    byte[][] gatherPointTypes;
    /** traces for this gather */
    public StsGatherTraceSet gatherTraceSet = null;
    /** fractional max error allowed in curve fitting wavelets */
    double maxError = 0.01f;
    /** sorted indices for absOffset axis display */
    int[] absOffsetSortedLineIndices = null;//	float maxAmplitude = 0.0f;
    /** estimated maxSemblanceCurve based on stat analysis */
    StsPoint[] maxSemblanceCurve;//	double[] traceOffsets = null;
    /** index location of trace nearest zero offset in this gather */
    int nGatherMinOffsetTrace = 0;
    /** range for each trace which is not muted */
    private double[][] muteRange = null;
    /** velocity trace from velocity volume */
    //    public float[] traceVelocities;
    /** interval and rms velocities computed for this gather over sample points inside velocity profile for this gather. */
    public Velocities velocities = new Velocities();
    /** input interval and rms velocities computed for this gather over sample points inside velocity profile for this gather. */
    public Velocities inputVelocities;
    /** lock object while computing gather data */
    protected boolean computingGather = false;
    /** semblance data is being computed */
    public boolean computingSemblance = false;
    /** trace offsets for all traces on line containing this gather */
    public double[] traceOffsets = null;
    /**
     * offset values where NMO curves, mutes, etc are drawn. Points are drawn at every offset at midpoint between every offset,
     * so the number of values is 2*nTraces-1.  For the pair of offsets straddling zero, the midpoint is not drawn but is replaced
     * by the value at zero.
     */
    public double[] drawOffsets;
    /** indicates files needed for this gather are available */
    public boolean fileOk = false;
    /** current velocityProfile */
    public StsVelocityProfile velocityProfile = null;
    /** input velocityProfile (if available) */
    StsVelocityProfile inputVelocityProfile = null;
    /** used in normalizing a trace */
    double maxAmplitude;
    /** indicates display is flattened */
    boolean flatten = false;
    /** scratch array for filteredData [nTraces][nSamples] */
    // double[][] data;
    /** scratch array for vrms-squared */
    double[] rVsqs;
    /** scratch array for t0-squared */
    double[] t0sqs;
    /** stretch mute curve with an offset for each t0 sample in a gather */
    double[] stretchMuteOffsets;
    /** stretch mute factor used in computing current stretchMuteTimes */
    double stretchFactor = -1;
    private float stretchMuteMinOffset;

    protected int nInterpIntervals;

    static final boolean debug = false;

    static final boolean debugVelTime = false;
    static final float debugVelocity = 6.5262f;
    static final float debugTime = 940.f;
    static int debugVelocityIndex;
    static int debugTimeIndex;
    static int debugSuperGatherRow;
    static int debugSuperGatherCol;

    public final static byte nullByte = StsParameters.nullByte;
    final static double roundOff = StsParameters.roundOff;

    public final static int TOP_MUTE_INDEX = -1;
    public final static int BOT_MUTE_INDEX = -2;
    public final static int STRETCH_MUTE_INDEX = -3;
    public final static int NONE_INDEX = -99;

    public final static byte ORDER_2ND = StsSemblanceComputeProperties.ORDER_2ND;
    public final static byte ORDER_4TH = StsSemblanceComputeProperties.ORDER_4TH;
    public final static byte ORDER_6TH = StsSemblanceComputeProperties.ORDER_6TH;
    public final static byte ORDER_6TH_OPT = StsSemblanceComputeProperties.ORDER_6TH_OPT;
    public final static byte ORDER_FAST = StsSemblanceComputeProperties.ORDER_FAST;

    static StsColor outsideMuteColor = StsColor.GRAY;
    //    transient StsColor muteColor = StsColor.GREEN;

    public StsGather()
    {
    }

    public StsGather(StsSuperGather superGather)
    {
        this.superGather = superGather;
        setSuperGather(superGather);
    }

    public void setSuperGather(StsSuperGather superGather)
    {
        this.superGather = superGather;
        this.velocityProfile = superGather.velocityProfile;
        this.inputVelocityProfile = superGather.inputVelocityProfile;
        lineSet = superGather.lineSet;
        setStretchMuteFactor();
        setStretchMuteMinOffset();
    }

    private void setStretchMuteMinOffset()
    {
        stretchMuteMinOffset = lineSet.getWiggleDisplayProperties().getStretchMuteMinOffset();
    }

    private void setStretchMuteFactor()
    {
        double percentStretch = lineSet.getWiggleDisplayProperties().getStretchMute();
        stretchFactor = percentStretch / 100 + 1;
    }

    public StsGather(StsPreStackLineSet lineSet)
    {
        this.lineSet = lineSet;
    }

    public boolean isInsideRowColRange()
    {
        return superGather.isInsideRowColRange(gatherRow, gatherCol);
    }

    /** clear data associated with this gather and return true if gather is ok for this new row and col */
    public boolean initializeRowCol(int row, int col)
    {
        boolean changed = gatherRow != row || gatherCol != col;
        if (!changed) return true;
        this.gatherRow = row;
        this.gatherCol = col;
        gatherDataChanged();
        return initializeGather();
    }

    public void initializeTimeRange()
    {
        velocities.initializeTimeRange();
        if (inputVelocities != null) inputVelocities.initializeTimeRange();
    }

    /** return true if gather hasn't changed or initialized ok */
    public boolean initializeLineGather(int nLine, int nGather)
    {
        boolean changed = false;
        if (lineSet.isInline)
        {
            changed = gatherRow != nLine || gatherCol != nGather;
            this.gatherRow = nLine;
            this.gatherCol = nGather;
        }
        else
        {
            changed = gatherRow != nGather || gatherCol != nLine;
            this.gatherRow = nGather;
            this.gatherCol = nLine;
        }
        if (!changed) return true;
        gatherDataChanged();
        return initializeGather();
    }

    public boolean checkInitializeGather()
    {
        if (initialized) return true;
        return initializeGather();
    }

    synchronized protected boolean initializeGather()
    {
        if (debug) debugThread("initializeGather()");
        flatten = lineSet.lineSetClass.getFlatten();
        nGatherTraces = 0;
        absOffsetSortedLineIndices = null;
        if (!lineSet.hasGatherTraces(gatherRow, gatherCol)) return false;
        line = lineSet.getDataLine(gatherRow, gatherCol);
        if (line == null) return false;
        if (lineSet.isInline)
            nLineGather = gatherCol - line.minGatherIndex;
        else
            nLineGather = gatherRow - line.minGatherIndex;
        if (nLineGather < 0 || nLineGather >= line.nCols) return false;

        if (nLineGather == 0)
            nFirstGatherLineTrace = 0;
        else
            nFirstGatherLineTrace = line.nLastTraceInGathers[nLineGather - 1] + 1;
        nLastGatherLineTrace = line.nLastTraceInGathers[nLineGather];
        nGatherTraces = nLastGatherLineTrace - nFirstGatherLineTrace + 1;
        nGatherMinOffsetTrace = findGatherMinOffsetTrace();

        nInterpIntervals = StsTraceUtilities.computeInterpolationInterval(line.getZInc(), 5);

		traceOffsets = line.getTraceOffsets();
        return true;
    }

    private void debugThread(String message)
    {
        System.out.println(Thread.currentThread().getName() + " class " + StsToolkit.getSimpleClassname(this) + " " + message +
            " computingSemblance " + computingSemblance + " gatherRow " + gatherRow + " gatherCol " + gatherCol);
    }

    public StsPoint[] getVelocityProfilePoints()
    {
        StsPreStackVelocityModel velocityModel = lineSet.getVelocityModel();
        if (velocityModel == null) return null;
        return velocityModel.getVelocityProfilePoints(gatherRow, gatherCol);
    }

    /** requires that velocity has been previously set and optionally the f factor. */
    public double computeOffsetTime(double x)
    {
        return velocities.nmoCoefficients.computeOffsetTime(x);
    }

    public boolean isComputingGather()
    {
        return computingGather;
    }

    public float[] computeStackedTrace(int nDisplaySlices)
    {
        int index;
        double f;
        if (line == null || line.traceOffsets == null)
        {
            //			System.out.println("System error! StsSuperGather.computeStackedTrace.traceOffsets are null.");
            return null;
        }
        try
        {
            float[] stackedData = new float[nDisplaySlices];
            if (StsPreStackLineSet3d.debugTimer) StsVolumeConstructorTimer.calcVelTimer.start();
            StsPoint[] velocityPoints = getVelocityProfilePoints();
            if (StsPreStackLineSet3d.debugTimer)
                StsVolumeConstructorTimer.calcVelTimer.stopAccumulateIncrementCountPrintInterval("getVelocityProfilePoints");
            if (velocityPoints == null) return null;

            checkComputeMuteRange();

            if (gatherTraceSet == null)
            {
                gatherTraceSet = new StsGatherTraceSet(this);
            }
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.start();
            int nSlices = line.nSlices;
            float zMin = line.zMin;
            float zInc = line.zInc;

            if (rVsqs == null || rVsqs.length < nDisplaySlices)
            {
                rVsqs = new double[nDisplaySlices];
                t0sqs = new double[nDisplaySlices];
            }

            float t0 = zMin;
            double[] vrms = velocities.computeVelocities(velocityProfile, zMin, zInc, nDisplaySlices);
            if (vrms == null) return null;

            for (int n = 0; n < nDisplaySlices; n++, t0 += zInc)
            {
                rVsqs[n] = 1.0 / (vrms[n] * vrms[n]);
                t0sqs[n] = t0 * t0;
            }

            double[] xsqs = new double[nGatherTraces];
            int nTrace = nFirstGatherLineTrace;
            for (int n = 0; n < nGatherTraces; n++, nTrace++)
            {
                double x = line.traceOffsets[nTrace];
                xsqs[n] = x * x;
            }

            double[] sumValues = new double[nDisplaySlices];
            int[] nValues = new int[nDisplaySlices];

            nTrace = nFirstGatherLineTrace;
            StsGatherTrace[] gatherTraces = gatherTraceSet.getGatherTraces();
            for (StsGatherTrace trace : gatherTraces)
            {
                double xsq = trace.xsq;
                t0 = zMin;
                double value;
                double v0, v1;

                for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
                {
                    double tOffset = Math.sqrt(t0sqs[s] + xsq * rVsqs[s]);
                    if (!isInsideMuteRange(trace, tOffset, t0)) continue;

                    if (!lineSet.isNMOed)
                        f = (tOffset - zMin) / zInc;
                    else
                        f = (t0 - zMin) / zInc;

                    if (f < 0.0 || f > nSlices - 1) continue;

                    index = StsMath.floor(f);
                    if (index < 0)
                    {
                        System.out.println("System error! StsSuperGather.computeStackedTrace.index < 0");
                    }
                    if (index == nSlices - 1)
                    {
                        value = trace.data[index];
                    }
                    else
                    {
                        f = f - index;
                        v0 = trace.data[index];
                        v1 = trace.data[index + 1];
                        value = v0 + f * (v1 - v0);
                    }
                    sumValues[s] += value;
                    nValues[s]++;
                }
            }
            for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
            {
                if (nValues[s] == 0) continue;
                stackedData[s] = (float) (sumValues[s] / nValues[s]);
            }
            //			float debugValue = stackedData[100];
            //			System.out.println("row " + row + " col " + col + " value " + debugValue);
            if (StsVolumeConstructor.runTimer)
            {
                StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("computeStackTrace");
                StsVolumeConstructorTimer.stackTimer.printElapsedTime();
                StsVolumeConstructorTimer.calcVelTimer.printElapsedTime();
            }
            return stackedData;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeStackedTrace", e);
            return null;
        }
    }

    public StsGatherTrace[] getConstructGatherTraces()
    {
        if (gatherTraceSet == null)
            gatherTraceSet = new StsGatherTraceSet(this);
        return gatherTraceSet.getGatherTraces();
    }

    public StsGatherTrace[] constructGatherTraces()
    {
        if(line == null)
            return null;
        int nSlices = line.nSlices;
        StsGatherTrace[] traces = new StsGatherTrace[nGatherTraces];
        float[] dataF = getGatherDataf();
        double[] traceOffsets = getGatherOffsets();
        if (dataF != null)
        {
            for(int t = 0, i = 0; t < nGatherTraces; t++, i += nSlices)
                traces[t] = new StsGatherTrace(dataF, i, nSlices, traceOffsets[t]);
        }
        return traces;
    }

    public float[] computeSemblanceTrace(int nDisplaySlices)
    {
        int index;
        double f;
        if (line == null || line.traceOffsets == null)
        {
            //			System.out.println("System error! StsSuperGather.computeStackedTrace.traceOffsets are null.");
            return null;
        }
        try
        {
            float[] stackedData = new float[nDisplaySlices];
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.calcVelTimer.start();
            StsPoint[] velocityPoints = getVelocityProfilePoints();
            if (StsVolumeConstructor.runTimer)
                StsVolumeConstructorTimer.calcVelTimer.stopAccumulateIncrementCountPrintInterval("getVelocityProfilePoints");
            if (velocityPoints == null) return null;

            checkComputeMuteRange();

            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.start();
            int nSlices = line.nSlices;
            float zMin = line.zMin;
            float zInc = line.zInc;

            if (rVsqs == null || rVsqs.length < nDisplaySlices)
            {
                rVsqs = new double[nDisplaySlices];
                t0sqs = new double[nDisplaySlices];
            }

            float t0 = zMin;
            double[] vrms = velocities.computeVelocities(velocityProfile, zMin, zInc, nDisplaySlices);
            if (vrms == null) return null;

            for (int n = 0; n < nDisplaySlices; n++, t0 += zInc)
            {
                rVsqs[n] = 1.0 / (vrms[n] * vrms[n]);
                t0sqs[n] = t0 * t0;
            }

            double[] xsqs = new double[nGatherTraces];
            int nTrace = nFirstGatherLineTrace;
            for (int n = 0; n < nGatherTraces; n++, nTrace++)
            {
                double x = line.traceOffsets[nTrace];
                xsqs[n] = x * x;
            }

            double[] sumValues = new double[nDisplaySlices];
            double[] sumSqValues = new double[nDisplaySlices];
            int[] nValues = new int[nDisplaySlices];

            int index0 = 0;
            StsGatherTrace[] gatherTraces = gatherTraceSet.getGatherTraces();
            for (StsGatherTrace trace : gatherTraces)
            {
                double xsq = trace.xsq;
                t0 = zMin;
                double value;
                double v0, v1;
                float[] data = trace.data;
                for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
                {
                    double tOffset = Math.sqrt(t0sqs[s] + xsq * rVsqs[s]);
                    if (!isInsideMuteRange(trace, tOffset, t0)) continue;

                    if (!lineSet.isNMOed)
                        f = (tOffset - zMin) / zInc;
                    else
                        f = (t0 - zMin) / zInc;

                    if (f < 0.0 || f > nSlices - 1) continue;

                    index = StsMath.floor(f);
                    if (index < 0)
                    {
                        System.out.println("System error! StsSuperGather.computeStackedTrace.index < 0");
                    }
                    if (index == nSlices - 1)
                    {
                        value = data[index];
                    }
                    else
                    {
                        f = f - index;
                        v0 = data[index];
                        v1 = data[index + 1];
                        value = v0 + f * (v1 - v0);
                    }
                    sumValues[s] += value;
                    sumSqValues[s] += value * value;
                    nValues[s]++;
                }
            }
            for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
            {
                if (nValues[s] >= 4 && sumSqValues[s] != 0.0)
                    stackedData[s] = (float) (sumValues[s] * sumValues[s] / sumSqValues[s] / nValues[s]);
            }
            //			float debugValue = stackedData[100];
            //			System.out.println("row " + row + " col " + col + " value " + debugValue);
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data");
            return stackedData;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeStackedTrace", e);
            if (StsVolumeConstructor.runTimer)
                StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data exception");
            return null;
        }
    }

    public void computeSemblanceTrace(StsSuperGather.Semblance semblance)
    {
        int index;
        double f;
        if (line == null || line.traceOffsets == null)
        {
            return;
        }
        try
        {
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.calcVelTimer.start();
            StsPoint[] velocityPoints = getVelocityProfilePoints();
            if (StsVolumeConstructor.runTimer)
                StsVolumeConstructorTimer.calcVelTimer.stopAccumulateIncrementCountPrintInterval("getVelocityProfilePoints");
            if (velocityPoints == null) return;

            float[] gatherData = getGatherDataf();
            if (gatherData == null) return;
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.start();
            int nSlices = line.nSlices;
            float zMin = line.zMin;
            float zInc = line.zInc;

            float t0 = zMin;
            int nDisplaySlices = semblance.nValues;
            for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
            {
                float v = StsMath.interpolateValue(velocityPoints, t0, 1, 0);
                double rvsq = 1.0 / (v * v);

                int index0 = 0;
                StsGatherTrace[] gatherTraces = gatherTraceSet.getGatherTraces();
                for (StsGatherTrace trace : gatherTraces)
                {
                    double x = trace.x;

                    double value;
                    double v0, v1;
                    double tOffset = Math.sqrt(t0 * t0 + x * x * rvsq);

                    if (!isInsideMuteRange(trace, tOffset, t0)) continue;

                    f = (tOffset - zMin) / zInc;
                    if (f < 0.0 || f > nSlices - 1) continue;

                    index = StsMath.floor(f);
                    if (index < 0)
                    {
                        System.out.println("System error! StsSuperGather.computeStackedTrace.index < 0");
                    }
                    if (index == nSlices - 1)
                    {
                        value = gatherData[index0 + index];
                    }
                    else
                    {
                        f = f - index;
                        v0 = gatherData[index0 + index];
                        v1 = gatherData[index0 + index + 1];
                        value = v0 + f * (v1 - v0);
                    }
                    semblance.addSemblanceValue(s, value);
                }
            }
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data");
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeSemblanceTrace", e);
            if (StsVolumeConstructor.runTimer)
                StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data exception");
        }
    }

    public class Velocities
    {
        /** interval velocities; one for each tInc from lineSet.zMin to lineSet.zMax. total number is nStandardValues */
        public double[] vInt = null;
        /** rms velocities; one for each tInc from lineSet.zMin to lineSet.zMax. total number is nStandardValues */
        public double[] vrms = null;
        ///** time datum for velocity profile; interval velocity is constant from here to tMin; generally the top of seismic data */
        //		public double tDatum;
        /** first time point on semblanceRange intervals after top point of profile */
        public double tMin;
        /** semblanceRangeProperties increment which is multiple of data increment */
        public double tInc;
        /** last time point on semblanceRange intervals before bottom point of profile */
        public double tMax;
        /** number of velocity points from tMin to tMax (inclusive) with increment tInc */
        public int nValues;

        /** order of NMO polynomial */
        public byte order = ORDER_2ND;
        /** subclasses for calculating NMO coefficients for various orders */
        protected NMOCoefficients nmoCoefficients = null;
        /** 2nd momemt of velocity profile */
        public double[] cnmo2 = null;
        /** 4th moment of velocity profile */
        public double[] cnmo4 = null;
        /** 6th moment of velocity profile */
        public double[] cnmo6 = null;
        /** index of current t0 value in arrays */
        public int index;
        /** increment interpolation factor for current t0 in arrays */
        public double indexInc;
        /** velocity at current index and indexInc */
        public double velocity;
        /** t zero at current index and indexInc */
        public double tzero;
        /** offset times for each tzero for each trace; dimensioned [nTraces][nSamples] */
        public OffsetTimes offsetTimes = new OffsetTimes();

        public Velocities()
        {
        }

        /*
             public void checkInitialize(int nValues, double tMin, double tInc)
             {
                 this.nValues = nValues;
                 this.tMin = tMin;
                 this.tInc = tInc;
                 this.tMax = tMin + (nValues - 1) * tInc;
             }
        */
        public void initialize(int index)
        {
            this.index = index;
            this.indexInc = 0;
            velocity = vrms[index];
            tzero = tMin + index * tInc;
            initialize();
        }

        private void initialize()
        {
            initialize(1.0);
        }

        public void initialize(double factor)
        {
            nmoCoefficients = getNMOCoefficients(order);
            nmoCoefficients.setFactor(factor);
            nmoCoefficients.initialize();
        }

        void setTZero(double tzero)
        {
            this.tzero = tzero;
        }

        void setVelocity(double v)
        {
            double indexF = StsMath.arrayIndexF(v, vrms);
            index = StsMath.floor(indexF);
            indexInc = indexF - index;
            velocity = v;
            tzero = tMin + indexF * tInc;
            initialize();
        }

        double setVelocityFromTZero(double t0)
        {
            return setVelocityFromTZero(t0, 1.0);
        }

        double setVelocityFromTZero(double t0, double factor)
        {
            tzero = t0;
            double indexF = (t0 - tMin) / tInc;
            if (indexF < 0.0)
            {
                index = 0;
                indexInc = 0;
                velocity = vrms[0];
            }
            else if (indexF >= nValues - 1)
            {
                index = nValues - 1;
                indexInc = 0;
                velocity = vrms[index];
            }
            else
            {
                index = StsMath.floor(indexF);
                indexInc = indexF - index;
                velocity = vrms[index] + indexInc * (vrms[index + 1] - vrms[index]);
            }
            initialize(factor);
            return velocity;
        }

        double computeOffsetTime(double offset)
        {
            return nmoCoefficients.computeOffsetTime(offset);
        }

        public boolean checkSetOrder()
        {
            byte order = lineSet.semblanceComputeProperties.order;
            return checkSetOrder(order);
        }

        public boolean checkSetOrder(byte order)
        {
            if (order == this.order && nmoCoefficients != null) return false;
            this.order = order;
            nmoCoefficients = getNMOCoefficients(order);
            return true;
        }

        protected boolean checkVelocityProfileChanged(StsVelocityProfile currentVelocityProfile)
        {
            // no useful references to velocityProfile, so commenting out for now
            // since this is called very late in the display loop and results in arrays being nulled
            // and calculation running all over again.  TJL 1/21/07

            // Added back in because without it BB semblance textrues do not update on add point or edit point.
            // Ever.... SAJ 1/22/07

            if (currentVelocityProfile == velocityProfile) return false;
            velocityProfile = currentVelocityProfile;
            clearVelocities();
            return true;
        }

        public void clearVelocities()
        {
            vInt = null;
            vrms = null;
            cnmo2 = null;
            cnmo4 = null;
            cnmo6 = null;
            clearOffsetTimes();
            //            traceVelocities = null;
            //            semblanceBackbone = null;
            if (debug) System.out.println("velocity profile changed: arrays cleared");
        }

        private void clearOffsetTimes()
        {
            if (offsetTimes == null) return;
            offsetTimes.times = null;
        }

        public double[] getIntVelocities()
        {
            return vInt;
        }

        public double[] getRmsVelocities()
        {
            return vrms;
        }

        // TODO: evaluate whether we need to rebuild complete velocity profile to get a tzero when we have scaled entire profile by velMult
        public double computeTzero(StsGatherTrace trace, double velMult, double targetTime)
        {
            return computeTzero(trace, targetTime / velMult);
        }

        public double[] getOffsetTimes(int nGatherTrace)
        {
            if (!checkComputeOffsetTimes()) return null;
            return offsetTimes.times[nGatherTrace];
        }

        /** given a particular trace in the gather, get the tzero time for this offset time in the array of toffset versus tzero) */
        public double computeTzero(StsGatherTrace trace, double toffset)
        {
            int traceIndex = trace.index;
            float indexF = StsMath.arrayIndexF(toffset, getOffsetTimes(traceIndex));
            return tMin + indexF * tInc;
        }

        public double computeTOffset(int nGatherTrace, double t0)
        {
            return offsetTimes.computeTOffset(nGatherTrace, t0);
            /*
                if(t0 < tMin)
                    return velocities.tMin;
                else if(t0 >= tMax)
                    return velocities.tMax;
                else
                {
                    double indexF = (t0 - tMin) / tInc;
                    int index = StsMath.floor(indexF);
                    double f = indexF - index;
                    double tt0 = offsetTimes.times[nGatherTrace][index];
                    double tt1 = offsetTimes.times[nGatherTrace][index + 1];
                    return tt0 + f * (tt1 - tt0);
                }
            */
        }

        public double computeTOffset(double offset, double t0)
        {
            return offsetTimes.computeTOffset(offset, t0);
        }

        public double computeTZero(double offset, double tOffset, double t0Estimated)
        {
            if (tOffset < tMin)
                return velocities.tMin;
            else if (tOffset >= tMax)
                return velocities.tMax;

            double[] gatherOffsets = getGatherOffsets();
            double offsetIndexF = StsMath.arrayIndexF(offset, gatherOffsets);
            int offsetIndex = (int) offsetIndexF;
            offsetIndex = StsMath.minMax(offsetIndex, 0, gatherOffsets.length - 2);
            double offsetF = offsetIndexF - offsetIndex;
            double tt0 = offsetTimes.computeTZeroFromT(offsetIndex, tOffset, offset, t0Estimated);
            double tt1 = offsetTimes.computeTZeroFromT(offsetIndex + 1, tOffset, offset, t0Estimated);
            return tt0 + offsetF * (tt1 - tt0);
        }

        public double getVelocityFromTZero(double t0, double offset)
        {
            if (vrms == null)
            {
                if (velocityProfile == null) return StsParameters.largeFloat;
                return velocityProfile.getVelocityFromTzero((float) t0);
            }
            double indexF = (t0 - tMin) / tInc;
            if (indexF <= 0.0)
                return vrms[0];
            int index = StsMath.floor(indexF);
            index = Math.min(index, nValues - 2);

            double f = indexF - index;
            return vrms[index] + f * (vrms[index + 1] - vrms[index]);
        }

        public void initializeTimeRange()
        {
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            tInc = rangeProperties.zInc;
            tMin = rangeProperties.zMin;
            tMax = rangeProperties.zMax;
            nValues = (int) Math.round((tMax - tMin) / tInc) + 1;
        }

        public boolean checkVelocityRangeChanged()
        {
            boolean changed = false;
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            double newTInc = rangeProperties.zInc;
            if (newTInc != superGather.semblanceZInc)
                System.out.println("Duh");
            if (newTInc != tInc)
            {
                tInc = newTInc;
                changed = true;
            }

            double newTMin = rangeProperties.zMin;
            if (newTMin != tMin)
            {
                tMin = newTMin;
                changed = true;
            }

            double newTMax = rangeProperties.zMax;
            if (newTMax != tMax)
            {
                tMax = newTMax;
                changed = true;
            }

            int newNValues = (int) Math.round((tMax - tMin) / tInc) + 1;
            if (nValues != newNValues)
            {
                nValues = newNValues;
                changed = true;
            }
            if (changed) clearVelocities();
            return changed;
        }

        public void velocityRangeChanged(Velocities velocities)
        {
            tInc = velocities.tInc;
            tMin = velocities.tMin;
            tMax = velocities.tMax;
            nValues = velocities.nValues;
        }

        public boolean recompute()
        {
            return vInt == null;
        }

        public synchronized boolean computeVelocities(StsVelocityProfile velocityProfile)
        {
            StsPoint[] points = getProfilePoints(velocityProfile);
            if (points == null) return false;
            checkSetOrder();
            clearOffsetTimes();

            //            if(vInt != null)return true;
            int nPoints = points.length;
            if (nPoints == 0) return false;
            if (debug)
                StsException.systemDebug(this, "computeVelocities", " last velocity is " + points[nPoints - 1].v[0]);
            try
            {
                //				checkClearSemblanceBackboneParameters(points);
                if (vInt == null || vInt.length != nValues)
                    vInt = new double[nValues];
                if (vrms == null || vrms.length != nValues)
                    vrms = new double[nValues];
                double vp1 = points[0].v[0];
                if (nPoints == 1)
                {
                    Arrays.fill(vInt, vp1);
                    Arrays.fill(vrms, vp1);
                    //TODO Overkill here: constant offset times for constant velocity are easy to compute
                    //TODO but this is only a single point so no big deal
                    if (!computeNMOTerms()) return false;
                    return true;
                }
                double tp1 = points[0].v[1];
                int nNextPoint = 1;
                double t1 = superGather.seismicVolumeZMin;
                double t0 = 0;
                int n = 0;
                for (n = 0; n < nValues; n++)
                {
                    t0 = t1;
                    t1 += tInc;
                    if (t1 > tp1) break;
                    vrms[n] = vp1;
                    vInt[n] = vp1;
                }
                t1 = t0;
                double tIncR = 1 / tInc;
                double vp0 = 0, tp0 = 0;
                double v1 = vp1;
                double v1sqt1 = v1 * v1 * t1 * tIncR;
                double dv = 0;
                for (; n < nValues; n++)
                {
                    t0 = t1;
                    t1 += tInc;
                    while (t1 > tp1 && nNextPoint < nPoints)
                    {
                        vp0 = vp1;
                        tp0 = tp1;
                        vp1 = points[nNextPoint].v[0];
                        tp1 = points[nNextPoint].v[1];
                        double f = (t1 - tp0) / (tp1 - tp0);
                        v1 = vp0 + f * (vp1 - vp0);
                        double dvdt = (vp1 - vp0) / (tp1 - tp0);
                        dv = dvdt * tInc;
                        nNextPoint++;
                    }
                    double v0sqt0 = v1sqt1;
                    v1sqt1 = v1 * v1 * t1 * tIncR;
                    vrms[n] = v1;
                    //vInt[n] = (float) Math.sqrt(v1sqt1 - v0sqt0);  //this equation makes vInt increase with depth
                    v1 += dv;
                    // a single vInt could be computed from this equation:
                    vInt[n] = Math.sqrt((vp1 * vp1 * tp1 - vp0 * vp0 * tp0) / (tp1 - tp0));    //users have strongly urged for single vInt between points - SWC 10/14/2009
                    //                    System.out.println("vInt: " + vInt[n] + " new vInt: " + vIntx);
                }
                if (!computeNMOTerms()) return false;
                return true;
            }
            catch (Exception e)
            {
                StsException.outputException("StsViewSemblance2d.checkComputeVelocities() failed.", e, StsException.WARNING);
                return false;
            }
        }

        public boolean computeOffsetTimes()
        {
            return offsetTimes.computeOffsetTimes();
        }

        public double[] computeVelocities(StsVelocityProfile velocityProfile, float zMin, float zInc, int nValues)
        {
            StsPoint[] points = getProfilePoints(velocityProfile);
            if (points == null) return null;
            int nPoints = points.length;
            try
            {
                double[] vrms = new double[nValues];
                double vp1 = points[0].v[0];
                double tp1 = points[0].v[1];
                int nNextPoint = 1;
                double t1 = zMin;
                double t0 = 0;
                int n;
                for (n = 0; n < nValues; n++)
                {
                    t0 = t1;
                    t1 += zInc;
                    if (t1 > tp1) break;
                    vrms[n] = vp1;
                }
                t1 = t0;
                double vp0, tp0;
                double v1 = vp1;

                double dv = 0;
                for (; n < nValues; n++)
                {
                    t0 = t1;
                    t1 += zInc;
                    while (t1 > tp1 && nNextPoint < nPoints)
                    {
                        vp0 = vp1;
                        tp0 = tp1;
                        vp1 = points[nNextPoint].v[0];
                        tp1 = points[nNextPoint].v[1];
                        double f = (t1 - tp0) / (tp1 - tp0);
                        v1 = vp0 + f * (vp1 - vp0);
                        double dvdt = (vp1 - vp0) / (tp1 - tp0);
                        dv = dvdt * zInc;
                        nNextPoint++;
                    }
                    vrms[n] = v1;
                    v1 += dv;
                }
                return vrms;
            }
            catch (Exception e)
            {
                StsException.outputException("StsViewSemblance2d.checkComputeVelocities() failed.", e, StsException.WARNING);
                return null;
            }
        }


        public synchronized boolean checkComputeConstantIntervalVelocities(StsVelocityProfile velocityProfile)
        {
            if (velocityProfile == null) return false;
            if (!checkComputeVelocities()) return false;

            if (vInt != null) return true;
            StsPoint[] points = getProfilePoints(velocityProfile);
            checkSetOrder();
            int nPoints = points.length;
            if (debug)
                System.out.println("Recomputing interval velocities. last velocity is " + points[nPoints - 1].v[0]);
            try
            {
                if (vInt == null || vInt.length != nValues)
                    vInt = new double[nValues];
                if (vrms == null || vrms.length != nValues)
                    vrms = new double[nValues];

                if (nPoints == 1)
                {
                    double v = points[0].v[0];
                    Arrays.fill(vInt, v);
                    Arrays.fill(vrms, v);
                    return true;
                }
                double t0 = superGather.seismicVolumeZMin;
                double t1 = points[0].v[1];
                double v0 = points[0].v[0];
                double v1 = v0;
                double t = t0;
                int n = 0;
                while (t < t1)
                {
                    vrms[n] = v0;
                    vInt[n] = v0;
                    t += tInc;
                    n++;
                }
                double vi = 0;
                for (int i = 1; i < nPoints; i++)
                {
                    t0 = t1;
                    t1 = points[i].v[1];
                    v0 = v1;
                    v1 = points[i].v[0];
                    vi = Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / (t1 - t0));
                    while (t < t1)
                    {
                        vInt[n] = vi;
                        double w = t0 / t;
                        vrms[n] = Math.sqrt(w * v0 * v0 + (1 - w) * vi * vi);
                        n++;
                        t += tInc;
                    }
                }
                t0 = t1;
                v0 = v1;
                for (; n < nValues; n++, t += tInc)
                {
                    vInt[n] = vi;
                    double w = t0 / t;
                    vrms[n] = Math.sqrt(w * v0 * v0 + (1 - w) * vi * vi);
                }
                if (!computeNMOTerms()) return false;
                offsetTimes.computeOffsetTimes();
                return true;
            }
            catch (Exception e)
            {
                StsException.outputException("StsViewSemblance2d.checkComputeVelocities() failed.", e, StsException.WARNING);
                return false;
            }
        }

        StsPoint[] getProfilePoints(StsVelocityProfile velocityProfile)
        {
            int nPoints = 0;
            StsPoint[] points = null;
            if (velocityProfile != null)
            {
                points = velocityProfile.getProfilePoints();
                if (points != null) return points;
            }
            points = getVelocityProfilePoints();
            return points;
        }

        public double[] checkComputeLimitIntervalVelocities(StsVelocityProfile velocityProfile)
        {
            StsPoint[] points = velocityProfile.getProfilePoints();
            if (points == null) return null;
            int nPoints = points.length;
            if (nPoints < 2) return null;

            try
            {
                int nSlices = line.nSlices;
                double[] vInt = new double[nSlices];
                double vp1 = points[0].v[0];
                double tp1 = points[0].v[1];
                int nextPoint = 1;
                float tInc = line.zInc;
                float t = line.zMin + tInc / 2;
                int n = 0;
                for (; n < nSlices; n++, t += tInc)
                {
                    if (t > tp1) break;
                    vInt[n] = vp1;
                }
                double vp0 = 0, tp0 = 0;
                double dvdt = 0;
                for (; n < nSlices; n++, t += tInc)
                {
                    while (t > tp1 && nextPoint < nPoints)
                    {
                        vp0 = vp1;
                        tp0 = tp1;
                        vp1 = points[nextPoint].v[0];
                        tp1 = points[nextPoint].v[1];
                        dvdt = (vp1 - vp0) / (tp1 - tp0);
                        nextPoint++;
                    }
                    double f = (t - tp0) / (tp1 - tp0);
                    double v = vp0 + f * (vp1 - vp0);
                    //					vInt[n] = Math.sqrt(v*v + 2*dvdt*v*t);
                    vInt[n] = Math.sqrt(v * v + 2 * dvdt * v * t + dvdt * dvdt * tInc * tInc);
                }
                return vInt;
            }
            catch (Exception e)
            {
                StsException.outputException("StsViewSemblance2d.checkComputeLimitIntervalVelocities() failed.", e, StsException.WARNING);
                return null;
            }
        }

        boolean computeNMOTerms()
        {
            //            if(cnmo2 != null) return true;
            //            if(order == ORDER_2ND) return true;
            if (debug) System.out.println("Recomputing NMO arrays with order " + order);
            try
            {
                if (cnmo2 == null || cnmo2.length != nValues)
                {
                    cnmo2 = new double[nValues];
                    cnmo4 = new double[nValues];
                    cnmo6 = new double[nValues];
                }
                double tsum = 0.0;
                double vsum2 = 0.0;
                double vsum4 = 0.0;
                double vsum6 = 0.0;
                double vrms2 = 0.0;
                double vrms4;
                double vrms6;
                for (int n = 0; n < nValues - 1; n++)
                {
                    double vi = vInt[n];
                    vsum2 += Math.pow(vi, 2) * tInc;
                    vsum4 += Math.pow(vi, 4) * tInc;
                    vsum6 += Math.pow(vi, 6) * tInc;
                    tsum += tInc;
                    vrms2 = vsum2 / tsum;
                    vrms4 = vsum4 / tsum;
                    vrms6 = vsum6 / tsum;
                    cnmo2[n] = 1.0 / vrms2;
                    cnmo4[n] = (vrms2 * vrms2 - vrms4) / (4 * tsum * tsum * Math.pow(vrms2, 4));
                    cnmo6[n] = (2 * vrms4 * vrms4 - vrms2 * vrms6 - vrms2 * vrms2 * vrms4) / (Math.pow(tsum, 4) * Math.pow(vrms2, 7));
                }
                if (debug) System.out.println("last vrms2 " + vrms2);
                return true;
            }
            catch (Exception e)
            {
                StsException.outputException("StsSuperGather.computeNMOTerms() failed.", e, StsException.WARNING);
                return false;
            }
        }

        StsGather.Velocities.NMOCoefficients getNMOCoefficients(byte order)
        {
            if (nmoCoefficients != null && nmoCoefficients.order == order)
                return nmoCoefficients;
            switch (order)
            {
                case ORDER_2ND:
                    return new StsGather.Velocities.NMOCoefficients2ndOrder();
                case ORDER_4TH:
                    return new StsGather.Velocities.NMOCoefficients4thOrder();
                case ORDER_6TH:
                    return new StsGather.Velocities.NMOCoefficients6thOrder();
                case ORDER_6TH_OPT:
                    return new StsGather.Velocities.NMOCoefficients6thOrderOpt();
                case ORDER_FAST:
                    return new StsGather.Velocities.NMOCoefficients2ndOrder();
                default:
                    return null;
            }
        }

        abstract class NMOCoefficients
        {
            double f2, f4, f6;
            double c2, c4, c6;
            byte order = -1;

            abstract void setFactor(double f);

            abstract void initialize();

            abstract double computeOffsetTime(double offset);
        }

        class NMOCoefficients2ndOrder extends StsGather.Velocities.NMOCoefficients
        {
            NMOCoefficients2ndOrder()
            {
                order = ORDER_2ND;
            }

            void setFactor(double f)
            {
                f2 = 1 / (f * f);
            }

            void initialize()
            {
                c2 = f2 / (velocity * velocity);
            }

            double computeOffsetTime(double x)
            {
                return Math.sqrt(tzero * tzero + x * x * c2);
            }
        }

        class NMOCoefficients4thOrder extends StsGather.Velocities.NMOCoefficients
        {
            NMOCoefficients4thOrder()
            {
                order = ORDER_4TH;
            }

            void setFactor(double f)
            {
                f2 = 1 / (f * f);
                f4 = f2 * f2;
            }

            void initialize()
            {
                if (indexInc == 0.0)
                {
                    c2 = f2 * cnmo2[index];
                    c4 = f4 * cnmo4[index];
                }
                else
                {
                    c2 = f2 * (cnmo2[index] + indexInc * (cnmo2[index + 1] - cnmo2[index]));
                    c4 = f4 * (cnmo4[index] + indexInc * (cnmo4[index + 1] - cnmo4[index]));
                }
            }

            double computeOffsetTime(double x)
            {
                double tSq = tzero * tzero + x * x * c2;
                double dtSq = Math.pow(x, 4) * c4;
                if (dtSq > -tSq) tSq += dtSq;
                return Math.sqrt(tSq);
            }
        }

        class NMOCoefficients6thOrder extends StsGather.Velocities.NMOCoefficients
        {
            NMOCoefficients6thOrder()
            {
                order = ORDER_6TH;
            }

            void setFactor(double f)
            {
                f2 = 1 / (f * f);
                f4 = f2 * f2;
                f6 = f2 * f4;
            }

            void initialize()
            {
                if (indexInc == 0.0)
                {
                    c2 = f2 * cnmo2[index];
                    c4 = f4 * cnmo4[index];
                    c6 = f6 * cnmo6[index];
                }
                else
                {
                    c2 = f2 * (cnmo2[index] + indexInc * (cnmo2[index + 1] - cnmo2[index]));
                    c4 = f4 * (cnmo4[index] + indexInc * (cnmo4[index + 1] - cnmo4[index]));
                    c6 = f6 * (cnmo6[index] + indexInc * (cnmo6[index + 1] - cnmo6[index]));
                }
            }

            double computeOffsetTime(double x)
            {
                double tSq, dtSq;
                try
                {
                    tSq = tzero * tzero + x * x * c2;
                    dtSq = Math.pow(x, 4) * c4;
                    if (dtSq < -tSq) return Math.sqrt(tSq);
                    tSq += dtSq;
                    dtSq = Math.pow(x, 6) * c6;
                    if (dtSq < -tSq) return Math.sqrt(tSq);
                    tSq += dtSq;
                    if (tSq < 0.0) return 0.0;
                    return Math.sqrt(tSq);
                }
                catch (Exception e)
                {
                    StsException.systemError("StsSuperGather.computeTOffset() failed.");
                    return 0.0;
                }
            }
        }

        class NMOCoefficients6thOrderOpt extends StsGather.Velocities.NMOCoefficients
        {
            NMOCoefficients6thOrderOpt()
            {
                order = ORDER_6TH_OPT;
            }

            void setFactor(double f)
            {
                f2 = 1 / (f * f);
                f4 = f2 * f2;
                f6 = f2 * f4;
            }

            void initialize()
            {
                if (indexInc == 0.0)
                {
                    c2 = f2 * cnmo2[index];
                    c4 = f4 * cnmo4[index];
                    c6 = f6 * cnmo6[index];
                }
                else
                {
                    c2 = f2 * (cnmo2[index] + indexInc * (cnmo2[index + 1] - cnmo2[index]));
                    c4 = f4 * (cnmo4[index] + indexInc * (cnmo4[index + 1] - cnmo4[index]));
                    c6 = f6 * (cnmo6[index] + indexInc * (cnmo6[index + 1] - cnmo6[index]));
                }
            }

            /** Sun paper version of optimized 6th order */
            double computeOffsetTime(double x)
            {
                double tSq = tzero * tzero + x * x * c2;
                double dtSq = Math.pow(x, 4) * c4;
                // factor of 4 below is from Bill's code
                // double dtSq = 4 * Math.pow(x, 4) * c4;
                if (dtSq < -tSq) return Math.sqrt(tSq);
                tSq += dtSq;
                double t = Math.sqrt(tSq);
                if (t < tzero) t = tzero;
                // Sun equation has an 8 in the denominator; Bill's parameter is 5.332; 32 is the
                t = t + 5.332 * Math.pow(x, 6) * c6 / (8.0 * t);
                return t;
            }

            // Bill's version of optimized 6th order: bad
            double computeOffsetTimeBillBad(double x)
            {
                double tSq = tzero * tzero + x * x * c2;
                // factor of 4 below is from Bill's code
                double dtSq = 4 * Math.pow(x, 4) * c4;
                if (dtSq < -tSq) return Math.sqrt(tSq);
                tSq += dtSq;
                double t = Math.sqrt(tSq);
                if (t < tzero) t = tzero;
                // Sun equation has an 8 in the denominator; Bill's parameter is 5.332; 32 is also required for Bill's equation
                t = t + 32 * 5.332 * Math.pow(x, 6) * c6 / (8.0 * t);
                return t;
            }
        }

        private void computeMaxSemblanceCurve(byte[] semblanceBytes, int nSlices, int nVel, float vMin, float vInc, float zMin, float zInc)
        {
            maxSemblanceCurve = new StsPoint[nSlices];
            float z = zMin;
            for (int s = 0; s < nSlices; s++, z += zInc)
            {
                int i = s;
                double avg = 0;
                double wtSum = 0;
                double vel = vMin;
                for (int v = 0; v < nVel; v++, i += nSlices, vel += vInc)
                {
                    int semblance = StsMath.signedByteToUnsignedInt(semblanceBytes[i]);
                    avg += vel * semblance;
                    wtSum += semblance;
                }
                avg /= wtSum;
                i = s;
                double std = 0;
                vel = vMin;
                for (int v = 0; v < nVel; v++, i += nSlices)
                {
                    int semblance = StsMath.signedByteToUnsignedInt(semblanceBytes[i]);
                    double dif = (avg - vel) * semblance;
                    std += dif * dif;
                }
                std = Math.sqrt(std / (wtSum * wtSum));
                float velRange = vInc * (nVel - 1);
                //			if(std < velRange/20.0f)
                maxSemblanceCurve[s] = new StsPoint((float) avg, z);
                //			else
                //				maxSemblanceCurve[s] = new StsPoint(0.0f, z);
            }
        }

        private void computeMaxSemblanceCurve(byte[] semblanceBytes, int nSemblancePoints, int nVel, float vMin, float vInc, double[] zeroOffsetTimes)
        {
            maxSemblanceCurve = new StsPoint[nSemblancePoints];
            for (int s = 0; s < nSemblancePoints; s++)
            {
                int i = s;
                double avg = 0;
                double wtSum = 0;
                double vel = vMin;
                for (int v = 0; v < nVel; v++, i += nSemblancePoints, vel += vInc)
                {
                    int semblance = StsMath.signedByteToUnsignedInt(semblanceBytes[i]);
                    avg += vel * semblance;
                    wtSum += semblance;
                }
                avg /= wtSum;
                i = s;
                double std = 0;
                vel = vMin;
                for (int v = 0; v < nVel; v++, i += nSemblancePoints)
                {
                    int semblance = StsMath.signedByteToUnsignedInt(semblanceBytes[i]);
                    double dif = (avg - vel) * semblance;
                    std += dif * dif;
                }
                std = Math.sqrt(std / (wtSum * wtSum));
                float velRange = vInc * (nVel - 1);
                //			if(std < velRange/20.0f)
                maxSemblanceCurve[s] = new StsPoint((float) avg, (float) zeroOffsetTimes[s]);
                //			else
                //				maxSemblanceCurve[s] = new StsPoint(0.0f, (float)semblanceTimes[s]);
            }
        }

        public void displayMaxSemblanceCurve(GL gl)
        {
            StsGLDraw.drawLineStrip2d(gl, StsColor.WHITE, maxSemblanceCurve);
            //		StsGLDraw.drawDottedLineStrip2d(gl, StsColor.GREEN, maxSemblanceCurve, 2);
        }

        public double getScaledValue(byte byteValue)
        {
            double f = (double) StsMath.signedByteToUnsignedInt(byteValue);
            return line.dataMin + (f / 254) * (line.dataMax - line.dataMin);
        }

        public class OffsetTimes
        {
            double[][] times;
            int nTimes;

            OffsetTimes()
            {
            }

            /**
             * working from the bottom of the trace up, for that offset assign the offsetTime for each tZero as long as
             * they are monotonically decreasing upward.  If not assign previous offsetTime.
             */
            public synchronized boolean computeOffsetTimes()
            {
                try
                {
                    computingGather = true;
                    velocities.checkSetOrder(order);
                    int nTimes = nValues;
                    if (times == null || times.length < nGatherTraces || times[0].length < nTimes)
                        times = new double[nGatherTraces][nTimes];
                    for (int i = nTimes - 1; i >= 0; i--)
                    {
                        initialize(i);
                        for (int g = 0, n = nFirstGatherLineTrace; g < nGatherTraces; g++, n++)
                        {
                            double offset = traceOffsets[n];
                            double time = nmoCoefficients.computeOffsetTime(offset);
                            if (i < nTimes - 1 && time > times[g][i + 1])
                                times[g][i] = times[g][i + 1];
                            else
                                times[g][i] = time;
                        }
                    }
                    return true;
                }
                catch (Exception e)
                {
                    StsException.outputException("OffsertTimes.computeOffsetTimes() failed.", e, StsException.WARNING);
                    return false;
                }
                finally
                {
                    computingGather = false;
                }
            }

            public double computeTOffset(int nTrace, double t0)
            {
                double t;
                if (t0 < tMin) return t0;

                if (t0 > tMax)
                    return times[nTrace][nValues - 1] + t0 - tMax;

                double indexF = (t0 - tMin) / tInc;
                int index = StsMath.floor(indexF);
                if (index >= nValues - 1) index--;
                double f = indexF - index;
                double tt0 = times[nTrace][index];
                double tt1 = times[nTrace][index + 1];
                if (tt0 == tt1) return tt0;
                return tt0 + f * (tt1 - tt0);
            }

            public double computeTOffset(double offset, double t0)
            {
                if (t0 < tMin)
                    return velocities.tMin;
                else if (t0 >= tMax)
                    return velocities.tMax;

                double[] gatherOffsets = getGatherOffsets();
                double offsetIndexF = StsMath.arrayIndexF(offset, gatherOffsets);
                int offsetIndex = (int) offsetIndexF;
                offsetIndex = StsMath.minMax(offsetIndex, 0, gatherOffsets.length - 2);
                double offsetF = offsetIndexF - offsetIndex;

                double indexF = (t0 - tMin) / tInc;
                int index = StsMath.floor(indexF);
                double f = indexF - index;
                double tt00 = offsetTimes.times[offsetIndex][index];
                double tt01 = offsetTimes.times[offsetIndex][index + 1];
                double tt0 = tt00 + f * (tt01 - tt00);
                double tt10 = offsetTimes.times[offsetIndex + 1][index];
                double tt11 = offsetTimes.times[offsetIndex + 1][index + 1];
                double tt1 = tt10 + f * (tt11 - tt10);
                return tt0 + offsetF * (tt1 - tt0);
            }

            public double computeOffsetTimeAtTraceIndexFandTzero(float traceIndexF, double t0)
            {
                double t;
                if (t0 < tMin) return t0;

                int traceIndex = (int) traceIndexF;
                double traceF = traceIndexF - traceIndex;
                if (t0 > tMax)
                {
                    double tt0 = times[traceIndex][nValues - 1] + t0 - tMax;
                    double tt1 = times[traceIndex + 1][nValues - 1] + t0 - tMax;
                    return tt0 + traceF * (tt1 - tt0);
                }
                else
                {
                    double t0IndexF = (t0 - tMin) / tInc;
                    int t0Index = (int) t0IndexF;
                    if (t0Index >= nValues - 1) t0Index--;
                    double t0f = t0IndexF - t0Index;
                    double tt0 = times[traceIndex][t0Index];
                    double tt1 = times[traceIndex][t0Index + 1];
                    double tt00 = tt0 + t0f * (tt1 - tt0);
                    tt0 = times[traceIndex + 1][t0Index];
                    tt1 = times[traceIndex + 1][t0Index + 1];
                    double tt11 = tt0 + t0f * (tt1 - tt0);
                    return tt00 + traceF * (tt11 - tt00);
                }
            }

            /**
             * Given a point defining offset and offsetTime, find an interpolated value in the offsetTimes table
             * for trace number nTrace along a common 2nd degree NMO curve (same tZero and velocity).
             * Start with an estimated tZero value (tz) and for the bounding interval in the table trace,
             * compute the tOffset at the given offset for both points of the bounding interval.  If these bound
             * the input value of offsetTime, interpolate tZero and return.  If not move up or down in the trace
             * until the bounding interval is found or the limits are met.
             */
            public double computeTZeroFromT(int nTrace, double t, double x, double tz)
            {
                double[] offsetTimes = times[nTrace]; // vector of offset times for this trace
                double tSq = t * t;
                double xSq = x * x;
                double indexF = (tz - tMin) / tInc;
                int index = StsMath.minMax(StsMath.floor(indexF), 0, nValues - 2);
                if (index >= nValues - 1) index--;
                double tz0 = tMin + index * tInc;
                double tzSq0 = tz0 * tz0;
                double t0 = offsetTimes[index];
                double tSq0 = t0 * t0;
                index++;
                double tz1 = tMin + index * tInc;
                double tzSq1 = tz1 * tz1;
                double t1 = offsetTimes[index];
                double tSq1 = t1 * t1;
                while (!StsMath.betweenInclusive(tSq, tSq0, tSq1))
                {
                    if (tSq > tSq1)
                    {
                        index++;
                        if (index >= nValues) return tz1;
                        tz0 = tz1;
                        tz1 += tInc;
                        tzSq1 = tz1 * tz1;
                        tSq0 = tSq1;
                        t1 = offsetTimes[index];
                        tSq1 = t1 * t1;
                    }
                    else // tSq < tSq0
                    {
                        index--;
                        if (index < 0) return tz0;
                        tz1 = tz0;
                        tz0 -= tInc;
                        tzSq1 = tzSq0;
                        tzSq0 = tz0 * tz0;
                        tSq1 = tSq0;
                        t0 = offsetTimes[index];
                        tSq0 = t0 * t0;
                    }
                }
                double f = (tSq - tSq0) / (tSq1 - tSq0);
                double tzSq = tzSq0 + f * (tzSq1 - tzSq0);
                // debug check: see if squared interpolation better than unsquared (it is slightly better)
                /*
                double vSq = vSq0 + f*(vSq1 - vSq0);
                double tNew = Math.sqrt(tzSq + xSq/vSq);
                double diff = tNew - t;
                StsException.systemDebug(this, "computeTZeroFromT", "Sq diff error: " + diff);
                double t0 = Math.sqrt(tSq0);
                double t1 = Math.sqrt(tSq1);
                f = (t - t0)/(t1 - t0);
                double v0 = Math.sqrt(vSq0);
                double v1 = Math.sqrt(vSq1);
                double v = v0 + f*(v1 - v0);
                tNew = Math.sqrt(tzSq + xSq/(v*v));
                diff = tNew - t;
                StsException.systemDebug(this, "computeTZeroFromT", "Sq diff error: " + diff);
                */
                return Math.sqrt(tzSq);
            }
        }

        //        /**
        //         * From Oz p.439, stretch mute factor s = dT/dt = T/t where T is offset time and t is tzero time.
        //         * The stretch mute is a percent stretch with 0% as no stretch (which is at x=0).
        //         * So stretch mute = (s - 1)*100.
        //         *
        //         * stretchMutePoints are at each offset and midpoint between offsets, so number of points is 2*nTraces-1.
        //         * For a set of traces, the mute values are at 0, 2, 4...., so they are at: 2*(nTrace-1)
        //         */
        //        public boolean checkComputeStretchMuteTimes()
        //        {
        //            if(!checkComputeVelocitiesAndOffsetTimes()) return false;
        //
        //            if (drawOffsets == null && !computeDrawOffsets())
        //            {
        //                StsException.systemError(this, "checkComputeStretchMuteTimes", "drawOffsets are null!");
        //                return false;
        //            }
        //            int nDrawOffsets = drawOffsets.length;
        //            if (nDrawOffsets != 2 * nGatherTraces - 1)
        //            {
        //                StsException.systemError(this, "getComputeStretchMuteTimes", "nDrawOffsets and nGatherTraces not consistent.");
        //                return false;
        //            }
        //            double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
        //            double smf = stretchMute / 100 + 1;  // this is the stretchMuteFactor
        //            if(debug) StsException.systemDebug(this, "checkComputeStretchMuteTimes", "stretchMuteTimes: " + stretchMuteTimes +
        //                      " smf: " + smf + " currentStretchMuteFactor: " + currentStretchMuteFactor);
        //            if(stretchMuteTimes != null && smf == currentStretchMuteFactor) return true;
        //
        //            boolean hasMuteTimes = (stretchMuteTimes != null);
        //            if (!hasMuteTimes)
        //            {
        //                stretchMuteTimes = new double[drawOffsets.length];
        //                Arrays.fill(stretchMuteTimes, StsParameters.doubleNullValue);
        //                hasMuteTimes = false;
        //            }
        //            // if stretchMuteTimes haven't previously been computed (hasMuteTimes==false), we use an estimate of estimatedStretchMuteTime to
        //            // compute the correct stretchMute time at a givem smf particular offset where we have the offsetTimes.
        //            // Once we have the first estimatedStretchMuteTime, we use it as the estimate for the next.
        //            // If we have previously computed stretchMuteTimes (hasMuteTimes==true), we use the current estimatedStretchMuteTime
        //            // to compute the correct one at the given smf.
        //            double estimatedStretchMuteTime = stretchMuteTimes[0];
        //            for (int t = 0, n = 0; t < nGatherTraces; t++, n += 2)
        //            {
        //                if (hasMuteTimes) estimatedStretchMuteTime = stretchMuteTimes[n];
        //                double[] offsetTimes = getOffsetTimes(t);
        //                // TODO possibly replace table lookup with search of NMO curve equation (tough for 4th & 6th order). TJL 3/13/09
        //                estimatedStretchMuteTime = computeStretchMuteTime(smf, estimatedStretchMuteTime, offsetTimes, t);
        //                stretchMuteTimes[n] = estimatedStretchMuteTime;
        //            }
        //            // interpolate the estimatedStretchMuteTime at an offset midway between each pair of traces
        //            for (int t = 0, n = 1; t < nGatherTraces-1; t++, n += 2)
        //            {
        //                double muteTime0 = stretchMuteTimes[n-1];
        //                double muteTime1 = stretchMuteTimes[n+1];
        //                stretchMuteTimes[n] = 0.5*(muteTime1 + muteTime0);
        //            }
        //            currentStretchMuteFactor = smf;
        //            return true;
        //        }
        /*
               public void adjustStretchMute(StsPoint pressedPick, StsPoint pick, StsVelocityProfile velocityProfile, StsPreStackLineSet lineSet)
               {
                   String offsetAxisType = lineSet.getWiggleDisplayProperties().getOffsetAxisType();
                   double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
                   double s = stretchMute/100 + 1;
                   float x = pick.v[0];
                   float offsetIndexF = getTraceIndexFromOffset(x);
                   int offsetIndex = (int)offsetIndexF;
                   float offsetF = offsetIndexF - offsetIndex;
                   double[] offsetTimes0 = getOffsetTimes(offsetIndex);
                   double[] offsetTimes1 = getOffsetTimes(offsetIndex + 1);
                   double stretchMuteTime0 = getStretchMuteTime(offsetIndex);
                   double stretchMuteTime1 = getStretchMuteTime(offsetIndex + 1);
                   double slope0 = computeStretchMuteSlope(offsetTimes0, stretchMuteTime0);
                   double slope1 = computeStretchMuteSlope(offsetTimes1, stretchMuteTime1);
                   double slope = slope0 + offsetF*(slope1 - slope0);
                   double dt = pick.v[1] - pressedPick.v[1];
                   double sNew = slope*dt + s;
                   double tPressed = pressedPick.v[1];
                   boolean flatten = lineSet.lineSetClass.getFlatten();
                   double c;
                   if (flatten)
                   {
                       double tOffset = velocities.computeTOffset(x, t);
                       c = tOffset / t;
                   }
                   else
                   {
                       c = stretchMute / 100 + 1;
                       double t0Estimated = t / c;
                       double t0 = velocities.computeTZero(x, t, t0Estimated);
                       t0 = Math.max(t0, 0.00001f);
                       c = t / t0;
                   }
                   double newStretchMute = Math.min((c - 1) * 100, StsWiggleDisplayProperties.stretchMuteMax);
                   if (debug) System.out.println("adjustStretchMute from " + stretchMute + " to: " + newStretchMute);
                   lineSet.getWiggleDisplayProperties().setStretchMute((float) newStretchMute);
               }
        */
        //        private double computeStretchMuteSlope(double[] offsetTimes, double t)
        //        {
        //            float indexF = StsMath.binarySearch(offsetTimes, t);
        //            int index = (int)indexF;
        //            double t0 = offsetTimes[index];
        //            double tz0 = tMin + tInc*index;
        //            double smf0 = t0/tz0;
        //            double t1 = offsetTimes[index+1];
        //            double tz1 = tz0 + tInc;
        //            double smf1 = t1/tz1;
        //            return (smf1 - smf0)/(t1 - t0);
        //        }

        //        public double getStretchMuteTime(int nGatherTrace)
        //        {
        //            return stretchMuteTimes[nGatherTrace*2 + 1];
        //        }

        //        public boolean checkComputeStretchMuteTimesX(StsPoint[] profilePoints)
        //        {
        //            if (drawOffsets == null)
        //            {
        //                System.out.println("drawOffsets are null!");
        //                return false;
        //            }
        //            int nDrawOffsets = drawOffsets.length;
        //            if (nDrawOffsets != 2 * nGatherTraces - 1)
        //            {
        //                StsException.systemError(this, "getComputeStretchMuteTimes", "nDrawOffsets and nGatherTraces not consistent.");
        //                return false;
        //            }
        //            double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
        //            double smf = stretchMute / 100 + 1;
        //            // double smf = Math.sqrt(sm * sm - 1);
        //
        //            if(stretchMuteTimes != null && smf == currentStretchMuteFactor) return true;
        //
        //            boolean hasMuteTimes = (stretchMuteTimes != null);
        //            if (!hasMuteTimes)
        //            {
        //                stretchMuteTimes = new double[drawOffsets.length];
        //                Arrays.fill(stretchMuteTimes, StsParameters.doubleNullValue);
        //                hasMuteTimes = false;
        //            }
        //
        //            double[] offsetTimes1 = getOffsetTimes(0);
        //            double muteTimeEst = stretchMuteTimes[0];
        //            double[] offsetTimes = new double[nValues];
        //            for (int t = 0, n = 0; t < nGatherTraces; t++)
        //            {
        //                if (hasMuteTimes) muteTimeEst = stretchMuteTimes[n];
        //                double[] offsetTimes0 = offsetTimes1;
        //                double stretchMuteTime = computeStretchMuteTime(smf, muteTimeEst, offsetTimes0, t);
        //                if (!hasMuteTimes) muteTimeEst = stretchMuteTime;
        //                stretchMuteTimes[n] = stretchMuteTime;
        //                n++;
        //                if (t == nGatherTraces - 1) break;
        //                offsetTimes1 = getOffsetTimes(t);
        //                if (hasMuteTimes) muteTimeEst = stretchMuteTimes[n];
        //                StsMath.interpolate(offsetTimes0, offsetTimes1, 0.5, nValues, offsetTimes);
        //                stretchMuteTime = computeStretchMuteTime(smf, muteTimeEst, offsetTimes, t);
        //                if (!hasMuteTimes) muteTimeEst = stretchMuteTime;
        //                stretchMuteTimes[n] = stretchMuteTime;
        //                n++;
        //            }
        //            currentStretchMuteFactor = smf;
        //            return true;
        //        }

        //        /** We are given a stretch mute factor s, an estimated offset time tOffsetEst, and a vector of offsetTimes
        //         *  spaced at equal intervals of tZero at a specified offset.
        //         *  We wish to find the offsetTime which corresponds to this stretch mute factor, s.
        //         *  From Oz, p 439.  stretch factor s = dT/dt = T/t where T is offset time and t is tzero time.
        //         *  s monotonically decreases with offsetTime: ds/dT = (1/t)*(1 - s*s); since s is always >= 1,
        //         *  ds/dT is always < 0.
        //         *  Given the estimated offsetTime, find the interval which bounds it.
        //         *  Compute the stretchMute factor s0 and s1 for the the two times defining the interval, t0 and t1.
        //         *  If s0 and s1 bound s, we've found the correct interval.  If not move up or down to the next
        //         *  interval using the observation above that ds/dt < 0.
        //         *  Once the correct interval is found, we want to solve for the offset time which gives us the
        //         *  desired value of s.  Assuming t and tz vary linearly in the interval, we have:
        //         *  t = t0 + f*(t1 - t0),  tz = tz0 + f*(tz1 - tz0), and s = t/tz
        //         *  Solve these for f and return the interpolated value of t.
        //         *  (These equations might be more linear in tz and t squared).
        //         *
        //         * @param s the stretch mute factor
        //         * @param tOffsetEst
        //         * @param offsetTimes spaced in constant increments of tZero from tMin to tMax
        //         * @param nOffset
        //         * @return
        //         */
        //        private double computeStretchMuteTime(double s, double tOffsetEst, double[] offsetTimes, int nOffset)
        //        {
        //            double tz;
        //            int index;
        //            if (tOffsetEst != StsParameters.doubleNullValue)
        //            {
        //                tz = StsMath.minMax(tOffsetEst / s, tMin, tMax);
        //                index = (int) Math.round((tz - tMin) / tInc);
        //            }
        //            else
        //            {
        //                index = nValues / 2;
        //            }
        //            index = StsMath.minMax(index, 0, nValues-1);
        //            tz = tMin + tInc * index;
        //            double t = offsetTimes[index];
        //            double s1 = t / tz;
        //            double t1 = t;
        //            double t0 = t1;
        //            double dtz, tz0, tz1;
        //            // if this s (s1) is smaller than the desired value (s) decrease offset time to find it
        //            // if larger, increase offset time to find it
        //            if (s1 < s)
        //            {
        //                // if we are at top of trace, use top interval; prev is at index 1 and next is at index 0
        //                // so start off with next set to index 1
        //                if(index == 0)
        //                {
        //                    index = 1;
        //                    tz = tMin + tInc;
        //                    t1 = offsetTimes[index];
        //                    s1 = t1/ tz;
        //                }
        //                while(index > 0)
        //                {
        //                    t0 = t1;
        //                    index--;
        //                    t1 = offsetTimes[index];
        //                    tz -= tInc;
        //                    s1 = t1/ tz;
        //                    if(s1 >= s) break;
        //                }
        //                dtz = -tInc;
        //                tz1 = tz;
        //                tz0 = tz + tInc;
        //            }
        //            else // s1 >= s
        //            {
        //                // if we are at bottom of trace, use bottom interval; prev is index nValues-2 and next is nValues-1
        //                // so set next to index nValues-2
        //                if(index == nValues - 1)
        //                {
        //                    index--;
        //                    tz = tMax - tInc;
        //                    t1 = offsetTimes[index];
        //                    s1 = t1/ tz;
        //                }
        //                while (index < nValues - 1)
        //                {
        //                    t0 = t1;
        //                    index++;
        //                    t1 = offsetTimes[index];
        //                    tz += tInc;
        //                    s1 = t1/ tz;
        //                    if(s1 < s) break;
        //                }
        //                dtz = tInc;
        //                tz1 = tz;
        //                tz0 = tz - tInc;
        //            }
        //            double f = (s*tz0 - t0)/(t1 - t0 - s*(tz1 - tz0));
        //            t = t0 + f * (t1 - t0);
        //            // tz = tz0 + f* dtz;
        //            /*
        //            // double check stretchMute
        //            double v0 = vrms[index];
        //            double v1 = vrms[index+1];
        //            double v = v0 + f*(v1-v0);
        //            double x = drawOffsets[2*nOffset];
        //            double tCheck = Math.sqrt(tz*tz + x*x/(v*v));
        //
        //            if(!StsMath.sameAs(t, tCheck, 0.001))
        //            {
        //                StsException.systemError(this, "computeStretchMuteTime", "tOffset error at offset: " + x + " tZero: " + tz + " tOffset: " + t + " tOffsetCheck: " + tCheck);
        //            }
        //            double sCheck = t/tz;
        //            if(!StsMath.sameAs(s, sCheck, 0.001))
        //            {
        //                StsException.systemError(this, "computeStretchMuteTime", "stretchMute calc error at offset: " + x + " tZero: " + tz + " tOffset: " + t + " stretchFactor: " + s + " stretchFactorCheck: " + sCheck);
        //            }
        //            */
        //        /*
        //            double sNew = t / tz;
        //            if(!StsMath.sameAs(1.0, sNew / s, 0.001))
        //            {
        //                System.out.println("s not computed properly at mute time: " + t + " s: " + s + " sNew: " + sNew);
        //            }
        //        */
        //            return t;
        //        }

        //        private double computeStretchMuteTimeX(double smf, double muteTimeEst, double[] offsetTimes)
        //        {
        //            double t0;
        //            int index;
        //            if (muteTimeEst != StsParameters.doubleNullValue)
        //            {
        //                t0 = StsMath.minMax(muteTimeEst / smf, tMin, tMax);
        //                index = (int) Math.round((t0 - tMin) / tInc);
        //            }
        //            else
        //            {
        //                index = nValues / 2;
        //            }
        //            index = StsMath.minMax(index, 0, nValues-1);
        //            t0 = tMin + tInc * index;
        //            double muteTime = offsetTimes[index];
        //            double smfNext = muteTime/t0;
        //            double muteTimeNext = muteTime;
        //            double smfPrev = smfNext;
        //            double muteTimePrev = muteTimeNext;
        //
        //            // if this smf (smfNext) is smaller than the desired value (smf) decrease offset time to to find it
        //            // if larger, increase offset time to find it
        //            if (smfNext < smf)
        //            {
        //                // if we are at top of trace, use top interval; prev is at index 0 and next is at index 1
        //                if(index == 0)
        //                {
        //                    index++;
        //                    t0 += tInc;
        //                    muteTimeNext = offsetTimes[index];
        //                    smfNext = muteTimeNext/t0;
        //                }
        //                else
        //                {
        //                    while (index > 0 && smfPrev < smf)
        //                    {
        //                        index--;
        //                        muteTimePrev = muteTimeNext;
        //                        muteTimeNext = offsetTimes[index];
        //                        smfPrev = smfNext;
        //                        t0 -= tInc;
        //                        smfNext = muteTimeNext/t0;
        //                    }
        //                }
        //            }
        //            else if(smfNext >= smf)
        //            {
        //                // if we are at bottom of trace, use bottom interval; prev is index nValues-1 and next is nValues-2
        //                if(index == nValues - 1)
        //                {
        //                    index--;
        //                    t0 -= tInc;
        //                    muteTimeNext = offsetTimes[index];
        //                    smfNext = muteTimeNext/t0;
        //                }
        //                else
        //                {
        //                    while (index < nValues - 1 && smfNext >= smf)
        //                    {
        //                        index++;
        //                        muteTimePrev = muteTimeNext;
        //                        muteTimeNext = offsetTimes[index];
        //                        smfPrev = smfNext;
        //                        t0 += tInc;
        //                        smfNext = muteTimeNext/t0;
        //                    }
        //                }
        //            }
        //            double f = (smf - smfPrev) / (smfNext - smfPrev);
        //            return muteTimePrev + f * (muteTimeNext - muteTimePrev);
        //        }
    }

    public synchronized boolean checkComputeGatherData(boolean computeWavelets)
    {
        try
        {
            if (gatherTraceSet != null) return true;
            if (debug) StsException.systemDebug(this, "computeGatherData", "Computing gatherData for gather view");
            computingGather = true;
            gatherTraceSet = new StsGatherTraceSet(this);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsPreStackLine.computeGatherData() failed.", e, StsException.WARNING);
            return false;
        }
        finally
        {
            computingGather = false;
        }
    }

    /**
     * create a vector points[nPoints][2]: time and amplitude for each point.
     * Insert the first point; if it is zero, insert another zero point just before the first non-zero point
     * but only if it is greater than the time value of the third point.  So if the first point is zero, and the second point
     * is not, we insert all points.  If the first n points are zero, we insert a zero at index 0 and at index n-1, and then
     * all points after that.
     *
     * @param planeData
     * @param index
     * @return
     */
    public double[][] getWigglePoints(Object planeData, int index)
    {
        int s = 0, n = 0;
        try
        {
            int nSlices = line.nSlices;
            //        System.out.println("index =" + index + " nCroppedSlices=" + nCroppedSlices);

            int minSlice = 0;
            double[][] points = new double[nSlices][3];
            float tMin = line.getZCoor(minSlice);
            float tInc = line.getZInc();
            int nPoints = nSlices;
            if (planeData instanceof byte[])
            {
                byte[] bytes = (byte[]) planeData;
                float horizScale = 1.0f / 254.0f;
                float t = tMin;
                for (s = 0; s < nSlices; s++, t += tInc)
                {
                    points[s][0] = t;
                    points[s][1] = getAmplitude(bytes[index++], horizScale);
                }
            }
            else // is float[]
            {
                float[] floats = (float[]) planeData;
                boolean started = false;
                double v1 = floats[index++];
                if (Math.abs(v1) < roundOff) v1 = 0;

                int nLastNonZero = 0;
                float t = tMin;
                // insert first point regardless
                points[0][0] = t;
                points[0][1] = v1;
                double max = v1;
                n = 1;
                s = 1;
                t += tInc;
                // find first non-zero point; if s > 1, insert a zero point just before it and then insert this point
                for (; s < nSlices; s++, t += tInc)
                {
                    double v = floats[index++];
                    if (Math.abs(v) < roundOff) v = 0;
                    if (v == 0) continue;
                    if (s > 1)
                    {
                        points[n][0] = t - tInc;
                        //    points[n][1] = 0;
                        n++;
                    }
                    points[n][0] = t;
                    points[n][1] = v;
                    max = Math.max(max, Math.abs(v));
                    nLastNonZero = n;
                    n++;
                    break;
                }
                s++;
                t += tInc;
                // insert all points, but keep track of last nonZero index
                for (; s < nSlices; s++, t += tInc)
                {
                    points[n][0] = t;
                    double v = floats[index++];
                    double vAbs = Math.abs(v);
                    if (vAbs > roundOff)
                    {
                        points[n][1] = v;
                        max = Math.max(max, vAbs);
                        nLastNonZero = n;
                    }
                    n++;
                }
                nPoints = n;
                int nTrimmedPoints = nLastNonZero + 1;
                if (nPoints == 1)
                {
                    points = new double[2][3];
                    points[0][0] = tMin;
                    points[1][0] = tMin + tInc * (nSlices - 1);
                    nPoints = 2;
                }
                else if (nTrimmedPoints < nSlices)
                {
                    points = (double[][]) StsMath.trimArray(points, nTrimmedPoints); // include a zero point at the start & end
                    nPoints = nTrimmedPoints;
                }
                if (max > 0)
                {
                    double rMax = 1.0 / max;
                    for (int i = 0; i < nPoints; i++)
                        points[i][1] *= rMax;
                }
            }

            // compute slopes (0 curvature as end conditions)
            for (n = 1; n < nPoints - 1; n++)
                points[n][2] = (points[n + 1][1] - points[n - 1][1]) / 2;
            points[0][2] = 1.5 * (points[1][1] - points[0][1]) - points[1][2] / 2;
            points[nPoints - 1][2] = 1.5 * (points[nPoints - 1][1] - points[nPoints - 2][1]) - points[nPoints - 2][2] / 2;
            return points;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getWigglePoints", " index " + index + " slice " + s + " n " + n, e);
            return null;
        }
    }

    public float[] getGatherDataf()
    {
        float[] fdata;
        try
        {
            if (line == null) return null;

            StsPreStackLineFileBlocks fileBlocks = line.fileBlocks;
            fdata = fileBlocks.getGatherDataf(StsObject.getCurrentModel().getProject().getBlocksMemoryManager(), 4, nLineGather, nFirstGatherLineTrace, nLastGatherLineTrace, line.nSlices);
            if (fdata == null) return null;

            int nCols = nLastGatherLineTrace - nFirstGatherLineTrace + 1;
            byte dataType = StsFilterProperties.PRESTACK;
            fdata = StsSeismicFilter.filter(dataType, fdata, nCols, line.nSlices, line.getZInc(), lineSet.filterProperties, lineSet.agcProperties, line.dataMin, line.dataMax);
            if (superGather != null && superGather.datumShift != 0.0f)
                fdata = applyDatumShift(fdata, nCols, line.nSlices, line.getZInc(), superGather.datumShift);
            return fdata;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getGatherDataf", e);
            return null;
        }
    }

    /**
     * Stack gather at this constant velocity.
     * For the stretchMute to be applied:
     * From Oz p.439, stretch mute factor smf = dT/dt = T/t where T is offset time and t is tzero time.
     * The stretch mute is a percent stretch with 0% as no stretch (which is at x=0).
     * So stretch mute = (s - 1)*100.
     * For each offset, compute the minimum tzero muteTime tzm. Times less than this are muted.
     * Since tsq = tzsq + xsq/vsq
     * muteTime = F*x where F = 1/(Math.sqrt(s*s - 1)*velocity)
     */
    public float[] computeConstantVelocityStackTrace(int nDisplaySlices, double velocity)
    {
        int index;
        double f;
        if (line == null || line.traceOffsets == null)
        {
            //			System.out.println("System error! StsSuperGather.computeStackedTrace.traceOffsets are null.");
            return null;
        }
        try
        {
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.start();

            maxAmplitude = 0.0;
            float[] gatherData = getGatherDataf();
            if (gatherData == null) return null;
            int nSlices = line.nSlices;
            float zMin = superGather.seismicVolumeZMin;
            //            float zMax = rangeProperties.zMax;
            float zInc = superGather.seismicVolumeZInc;
            float t0 = zMin;

            double rvsq = 1.0 / (velocity * velocity);

            double[] sumValues = new double[nDisplaySlices];
            int[] nValues = new int[nDisplaySlices];
            float[] stackedData = new float[nDisplaySlices];

            int nTrace = nFirstGatherLineTrace;
            int index0 = 0;
            boolean isNMOed = lineSet.getIsNMOed();

            StsGatherTrace[] gatherTraces = gatherTraceSet.getGatherTraces();
            for (StsGatherTrace trace : gatherTraces)
            {
                double x = Math.abs(line.traceOffsets[nTrace++]);
                double xovsq = x * x * rvsq;
                t0 = zMin;
                double value;
                double v0, v1;
                for (int s = 0; s < nDisplaySlices; s++, t0 += zInc)
                {
                    double tOffset = Math.sqrt(t0 * t0 + xovsq);
                    if (!isInsideMuteRange(trace, tOffset, t0)) continue;
                    if (!isNMOed)
                        f = (tOffset - zMin) / zInc;
                    else
                    {
                        float t0i = (float) inputVelocities.computeTZero(x, tOffset, t0);
                        f = (t0i - zMin) / zInc;
                    }

                    if (f < 0.0 || f > nSlices - 1) continue;

                    index = StsMath.floor(f);
                    if (index < 0)
                    {
                        System.out.println("System error! StsSuperGather.computeStackedTrace.index < 0");
                    }
                    if (index == nSlices - 1)
                    {
                        value = gatherData[index0 + index];
                    }
                    else
                    {
                        f = f - index;
                        v0 = gatherData[index0 + index];
                        v1 = gatherData[index0 + index + 1];
                        value = v0 + f * (v1 - v0);
                    }
                    sumValues[s] += value;
                    nValues[s]++;
                }
            }
            for (int s = 0; s < nDisplaySlices; s++)
            {
                if (nValues[s] > 0) stackedData[s] = ((float) sumValues[s] / nValues[s]);
                // maxAmplitude = Math.max(maxAmplitude, Math.abs(stackedData[s]));
            }

            //StsMath.normalizeAmplitude(stackedData, (float) maxAmplitude);
            if (StsVolumeConstructor.runTimer) StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data");
            return stackedData;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSuperGather.computeConstantVelocityStackTrace() failed.", e, StsException.WARNING);
            if (StsVolumeConstructor.runTimer)
                StsVolumeConstructorTimer.stackTimer.stopAccumulateIncrementCountPrintInterval("stack data exception");
            return null;
        }
    }

    private int findGatherMinOffsetTrace()
    {
        if (nGatherTraces <= 0) return -1;
        if (nGatherTraces == 1) return 0;

        //		double[] traceOffsets = line.traceOffsets;
        if (line.traceOffsets == null) return -1;
        int nTrace = nGatherTraces / 2;
        if (nTrace + nFirstGatherLineTrace > nLastGatherLineTrace || nTrace + nFirstGatherLineTrace >= line.traceOffsets.length)
        {
            System.out.println("Bad juju: traceOffsets.length " + line.traceOffsets.length);
            System.out.println("   nFirstGatherTrace " + nFirstGatherLineTrace + " nLastGatherTrace " + nLastGatherLineTrace);
            System.out.println(" nTrace " + nTrace + " nGathers " + nGatherTraces);
            return -1;
        }
        double offset = line.traceOffsets[nTrace + nFirstGatherLineTrace];
        if (offset < 0.0f)
        {
            while (nTrace < nGatherTraces - 1)
            {
                nTrace++;
                double nextOffset = line.traceOffsets[nTrace + nFirstGatherLineTrace];
                if (nextOffset > 0.0f)
                {
                    if (nextOffset <= -offset)
                    {
                        return nTrace;
                    }
                    else
                    {
                        return nTrace - 1;
                    }
                }
                if (nTrace == (nGatherTraces - 1))
                    return nTrace;
                offset = nextOffset;
            }
        }
        else if (offset > 0.0f)
        {
            while (nTrace > 0)
            {
                nTrace--;
                double nextOffset = line.traceOffsets[nTrace + nFirstGatherLineTrace];
                if (nextOffset < 0.0f)
                {
                    if (-nextOffset <= offset)
                    {
                        return nTrace;
                    }
                    else
                    {
                        return nTrace + 1;
                    }
                }
                if (nTrace == 0)
                    return nTrace;
                offset = nextOffset;
            }
        }
        return nTrace;
    }

    /** the gather row or col have changed: null out all arrays */
    public void gatherDataChanged()
    {
        gatherTraceSet = null;
        //        traceVelocities = null;
    }

    /** trace data at this gather has been changed: null out gather data */
    public void gatherTracesChanged()
    {
        gatherTraceSet = null;
    }

    // TODO check that datumShift is not being applied if no datumShift is available
    public float[] applyDatumShift(float[] data, int nTraces, int nSlices, float sampleInterval, double datumShift)
    {
        // Shift all the traces and clip the bottoms
        int padSamples = Math.abs((int) Math.round(datumShift / sampleInterval));
        if (padSamples == 0) return data;
        float[] padTrace = new float[padSamples];
        for (int i = 0; i < padSamples; i++)
            padTrace[i] = 0.0f;

        int bufferOffset = 0;
        int traceLeft = nSlices - padSamples;
        if (datumShift > 0)
        {
            for (int i = 0; i < nTraces; i++)
            {
                // Shift valid samples down by datum shift amount
                System.arraycopy(data, bufferOffset, data, bufferOffset + padSamples, traceLeft);
                // Add 0.0 samples to top of trace
                System.arraycopy(padTrace, 0, data, bufferOffset, padSamples);
                bufferOffset = bufferOffset + nSlices;
            }
        }
        else
        {
            for (int i = 0; i < nTraces; i++)
            {
                // Shift valid samples up by datum shift amount
                System.arraycopy(data, bufferOffset + padSamples, data, bufferOffset, traceLeft);
                // Add 0.0 samples to bottom of trace
                System.arraycopy(padTrace, 0, data, bufferOffset + traceLeft, padSamples);
                bufferOffset = bufferOffset + nSlices;
            }
        }
        return data;
    }

    public double getStretchMuteTZeroFromOffset(double offset)
    {
        if (stretchMuteOffsets == null || lineSet == null) return 0;
        double t0Inc = lineSet.getZInc();
        double t0Min = lineSet.getZMin();

        for (int i = 0; i < stretchMuteOffsets.length; i++)
        { //offsets start small and get bigger
            if (offset < stretchMuteOffsets[i])
            {
                return t0Min + i * t0Inc;
            }
        }
        return lineSet.getZMax(); //offset is off the charts, so mute whole trace!!
    }

    /** computes velocity arrays for currentVelocities and inputVelocities (if exists) */
    public boolean checkComputeVelocities()
    {
        boolean profileChanged = velocities.checkVelocityProfileChanged(velocityProfile);
        boolean rangeChanged = velocities.checkVelocityRangeChanged();
        if (profileChanged || rangeChanged || velocities.recompute())
            if (!velocities.computeVelocities(velocityProfile)) return false;
        if (lineSet.inputVelocityModel == null) return true;
        if (inputVelocities == null)
        {
            inputVelocities = new Velocities();
            inputVelocities.initializeTimeRange();
        }

        if (!rangeChanged && !inputVelocities.recompute()) return true;
        if (rangeChanged) inputVelocities.velocityRangeChanged(velocities);
        return inputVelocities.computeVelocities(inputVelocityProfile);
    }

    private boolean checkComputeOffsetTimes()
    {
        if (velocities.offsetTimes.times == null)
        {
            if (!velocities.computeOffsetTimes()) return false;
        }
        if (lineSet.inputVelocityModel == null) return true;
        if (inputVelocities == null) return true;
        if (inputVelocities.offsetTimes.times != null) return true;
        return inputVelocities.computeOffsetTimes();
    }

    public boolean checkComputeVelocitiesAndOffsetTimes()
    {
        if (!checkComputeVelocities()) return false;
        return checkComputeOffsetTimes();
    }

    static int nMaxIterations = 100;

    public void checkAdjustFlattening(double[][] drawPoints, StsGatherTrace trace)
    {
        if (velocityProfile == null) return;

        int nDrawPoints = drawPoints.length;
        if (flatten) // is not NMOed: apply the current velocity profile NMO (velocities.offsetTimes.times table)
        {
            if (lineSet.isNMOed) // is NMOed: deNMO using initialVelocities.offsetTimes table and NMO using velocities.offsetTimes
            {
                double[] initialTimes = inputVelocities.getOffsetTimes(trace.index);
                for (int n = 1; n < nDrawPoints; n++)
                {
                    float t0 = (float) drawPoints[n][0];
                    if (t0 < velocities.tMin)
                        drawPoints[n][0] = velocities.tMin;
                    else if (t0 >= velocities.tMax)
                        drawPoints[n][0] = velocities.tMax;
                    else
                    {
                        double indexF = (t0 - velocities.tMin) / velocities.tInc;
                        int index = StsMath.floor(indexF);
                        double f = indexF - index;
                        double tOffset = initialTimes[index] + f * (initialTimes[index + 1] - initialTimes[index]);
                        drawPoints[n][0] = velocities.computeTzero(trace, tOffset);
                        //                        if(nGatherTrace == 10)
                        //                            System.out.println("flattened and NMOed point: t0 " + t0 + " tOffset " + tOffset + " adjusted t0 " + drawPoints[n][0]);
                    }
                }
            }
            else // not NMOed: NMO using velocities.offsetTimes table
            {
                double[] times = velocities.getOffsetTimes(trace.index);
                int nTimes = times.length;
                double t0 = times[0];
                double t1 = times[1];
                int i = 2;
                for (int n = 1; n < nDrawPoints; n++)
                {
                    // Added 0,0,0 point so trace draws to top
                    float tx = (float) drawPoints[n][0];
                    while ((t0 == t1 || tx > t1) && i < nTimes - 1)
                    {
                        t0 = t1;
                        t1 = times[i++];
                    }
                    if (t0 == t1) continue;
                    if (tx < t0)
                        drawPoints[n][0] = 0.0;
                    else
                    {
                        double indexF = i - 2 + (tx - t0) / (t1 - t0);
                        drawPoints[n][0] = velocities.tMin + indexF * velocities.tInc;
                    }
                }
            }
        }
        else
        if (lineSet.isNMOed) // unflatten: if not NMOed, do nothing; otherwise if NMOed, deNMO with inputVelocities.offsetTimes table
        {
            double[] initialTimes = inputVelocities.getOffsetTimes(trace.index);
            for (int n = 1; n < nDrawPoints; n++)
            {
                float t0 = (float) drawPoints[n][0];
                if (t0 < velocities.tMin)
                    drawPoints[n][0] = velocities.tMin;
                else if (t0 >= velocities.tMax)
                    drawPoints[n][0] = velocities.tMax;
                else
                {
                    double indexF = (t0 - velocities.tMin) / velocities.tInc;
                    int index = StsMath.floor(indexF);
                    double f = indexF - index;
                    drawPoints[n][0] = initialTimes[index] + f * (initialTimes[index + 1] - initialTimes[index]);
                }
            }
        }
    }

    protected double[] getOffsetTimes(int nGatherTrace)
    {
        return velocities.getOffsetTimes(nGatherTrace);
    }

    private void computeSortedTraceIndicesByAbsOffset()
    {
        if (nGatherTraces == 0) return;
        absOffsetSortedLineIndices = new int[nGatherTraces];
        int nLineMinOffsetTrace = nGatherMinOffsetTrace + nFirstGatherLineTrace;
        absOffsetSortedLineIndices[0] = nLineMinOffsetTrace;
        int nMinus = nLineMinOffsetTrace - 1;
        int nPlus = nLineMinOffsetTrace + 1;
        int n = 1;

        double minusOffset;
        double plusOffset;

        try
        {
            if (nMinus >= nFirstGatherLineTrace)
                minusOffset = Math.abs(line.traceOffsets[nMinus]);
//			minusOffset = Math.abs(line.traceOffsets[nMinus]);
            else
                minusOffset = StsParameters.largeFloat;

            if (nPlus <= nLastGatherLineTrace)
                plusOffset = Math.abs(line.traceOffsets[nPlus]);
//			plusOffset = Math.abs(line.traceOffsets[nPlus]);
            else
                plusOffset = StsParameters.largeFloat;

            while (n < nGatherTraces)
            {
                if (minusOffset < plusOffset)
                {
                    absOffsetSortedLineIndices[n++] = nMinus--;
                    if (nMinus >= nFirstGatherLineTrace)
                        minusOffset = Math.abs(line.traceOffsets[nMinus]);
//					minusOffset = Math.abs(line.traceOffsets[nMinus]);
                    else
                        minusOffset = StsParameters.largeFloat;
                }
                else if (plusOffset != StsParameters.largeFloat)
                {
                    absOffsetSortedLineIndices[n++] = nPlus++;
                    if (nPlus <= nLastGatherLineTrace)
                        plusOffset = Math.abs(line.traceOffsets[nPlus]);
//					plusOffset = Math.abs(line.traceOffsets[nPlus]);
                    else
                        plusOffset = StsParameters.largeFloat;
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSuperGather.computeSortedTraceIndicesByAbsOffset() failed. gatherRow " + gatherRow + " gatherCol " + gatherCol +
                " nFirstGatherTrace: " + nFirstGatherLineTrace + " nLastGatherTrace " + nLastGatherLineTrace +
                " nMinus " + nMinus + " nPlus " + nPlus, e, StsException.WARNING);
        }
    }

    private int[] getTraceDisplayRange(float[] offsetDisplayRange, int offsetAxisType)
    {
        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
        {
            if (traceOffsets == null) return null;

            int nFirstDisplayTrace, nLastDisplayTrace;
            for (nFirstDisplayTrace = nFirstGatherLineTrace; nFirstDisplayTrace < nLastGatherLineTrace; nFirstDisplayTrace++)
                if (traceOffsets[nFirstDisplayTrace] >= offsetDisplayRange[0]) break;
            for (nLastDisplayTrace = nLastGatherLineTrace; nLastDisplayTrace > nFirstGatherLineTrace; nLastDisplayTrace--)
                if (traceOffsets[nLastDisplayTrace] <= offsetDisplayRange[1]) break;
            int nFirstGatherIndex = nFirstDisplayTrace - nFirstGatherLineTrace;
            int nLastGatherIndex = nLastDisplayTrace - nFirstGatherLineTrace;
            return new int[]
                {nFirstGatherIndex, nLastGatherIndex};
        }
        else // offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_ABS_INDEX
        {
            int nFirstTrace, nLastTrace;
            for (nFirstTrace = 0; nFirstTrace < nGatherTraces; nFirstTrace++)
                if (nFirstTrace >= offsetDisplayRange[0]) break;
            for (nLastTrace = nGatherTraces - 1; nLastTrace >= 0; nLastTrace--)
                if (nLastTrace <= offsetDisplayRange[1]) break;
            return new int[]
                {nFirstTrace, nLastTrace};
        }
    }

    private int[] getTraceOffsetDisplayRange(float[] offsetDisplayRange)
    {
        if (traceOffsets == null) return null;

        int nFirstDisplayTrace, nLastDisplayTrace;
        for (nFirstDisplayTrace = nFirstGatherLineTrace; nFirstDisplayTrace < nLastGatherLineTrace; nFirstDisplayTrace++)
            if (traceOffsets[nFirstDisplayTrace] >= offsetDisplayRange[0]) break;
        for (nLastDisplayTrace = nLastGatherLineTrace; nLastDisplayTrace > nFirstGatherLineTrace; nLastDisplayTrace--)
            if (traceOffsets[nLastDisplayTrace] <= offsetDisplayRange[1]) break;
        return new int[]
            {nFirstDisplayTrace, nLastDisplayTrace};
    }

    private int[] getTraceIndexDisplayRange(float[] offsetDisplayRange)
    {
        int nFirstDisplayTrace, nLastDisplayTrace, nTrace;
        for (nTrace = 0, nFirstDisplayTrace = nFirstGatherLineTrace; nFirstDisplayTrace < nLastGatherLineTrace; nFirstDisplayTrace++, nTrace++)
            if (nTrace >= offsetDisplayRange[0]) break;
        for (nTrace = nGatherTraces - 1, nLastDisplayTrace = nLastGatherLineTrace; nLastDisplayTrace > nFirstGatherLineTrace; nLastDisplayTrace--, nTrace--)
            if (nTrace <= offsetDisplayRange[1]) break;
        return new int[]
            {nFirstDisplayTrace, nLastDisplayTrace};
    }

    private float getAmplitude(byte byteValue, float horizScale)
    {
        if (byteValue == nullByte) return 0.0f;
        int signedInt = StsMath.unsignedByteToSignedInt(byteValue);
        return horizScale * (signedInt - line.signedIntDataZero);
    }

    /**
     * From Oz p.439, stretch factor s = dT/dt = T/t where T is offset time and t is tzero time.
     * So s = SQRT(1 + x*x/(t*t*v*v).
     * The stretch mute is a percent stretch with 0% as no stretch (which is at x=0).
     * So stretch mute = (s - 1)*100.
     *
     * @param pick
     * @param velocityProfile
     * @param lineSet
     */
    public void adjustStretchMute(StsPoint pressedPick, StsPoint pick, StsVelocityProfile velocityProfile, StsPreStackLineSet lineSet)
    {
        double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
        float x = pick.v[0];
        x = getOffsetForTypeValue(x);
        double t = pick.v[1];
        if (debug) StsException.systemDebug(this, "adjustStretchMute", "offset: " + x + " offsetTime: " + t);
        boolean flatten = lineSet.lineSetClass.getFlatten();
        double c;
        if (flatten)
        {
            double tOffset = velocities.computeTOffset(x, t);
            c = tOffset / t;
        }
        else
        {
            c = stretchMute / 100 + 1;
            double t0Estimated = t / c;
            double t0 = velocities.computeTZero(x, t, t0Estimated);
            if (debug)
                StsException.systemDebug(this, "adjustStretchMute", "t0Estimated: " + t0Estimated + " t0: " + t0);
            t0 = Math.max(t0, 0.00001f);
            double cNew = t / t0;
            if (debug)
                StsException.systemDebug(this, "adjustStretchMute", "stretchMuteFactor adjusted from " + c + " to: " + cNew);
            c = t / t0;
        }
        double newStretchMute = Math.min((c - 1) * 100, StsWiggleDisplayProperties.stretchMuteMax);
        if (debug)
            StsException.systemDebug(this, "adjustStretchMute", "stretchMute adjusted from " + stretchMute + " to: " + newStretchMute);
        lineSet.getWiggleDisplayProperties().setStretchMute((float) newStretchMute);
        stretchMuteChanged();
    }

    /*
        protected void checkMuteChanged(StsVelocityProfile velocityProfile)
        {
            if(velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
                muteRange = null;
        }

        public void checkComputeMuteRange()
        {
            if(muteRange == null) checkComputeMuteRange();
        }
    */

    public void muteRangeChanged()
    {
        muteRange = null;
    }

    public boolean checkComputeMuteRange()
    {
        StsPoint topMute, bottomMute;
        muteRange = null;  //temp. hack to fix mute bug - somehow top/bottom mute getting set to bad values and muting whole record
        return true;
        /*
        if (muteRange != null) return true;
        if (lineSet.velocityModel == null)
        {
            muteRange = null;
            return false;
        }
        if (velocityProfile == null || velocityProfile.getNProfilePoints() == 0)
        {
            muteRange = null;
            return false;
        }

        if (nGatherTraces <= 0) return false;

        topMute = velocityProfile.getTopMute();
        bottomMute = velocityProfile.getBottomMute();

        // hacque:  ignore if mute vertical range is zero
        if (bottomMute.v[0] == StsParameters.nullValue)
        {
            muteRange = null;
            return false;
        }
        if (drawOffsets == null) computeDrawOffsets();

        muteRange = new double[nGatherTraces][2];
        double[] traceOffsets = line.getTraceOffsets();
        int nTrace = nFirstGatherLineTrace;
        for (int n = 0; n < nGatherTraces; n++)
        {
            double offset = traceOffsets[nTrace++];

            double topMuteTime, botMuteTime;

            if (topMute != null)
                topMuteTime = computeNmoMute(topMute, offset);
            else
                topMuteTime = 0.0;

            if (bottomMute != null)
                botMuteTime = computeNmoMute(bottomMute, offset);
            else
                botMuteTime = lineSet.zMax;

            muteRange[n][0] = topMuteTime;
            muteRange[n][1] = botMuteTime;
        }
        return true;
        */
    }

//    public boolean isInsideMuteRange(int nGatherTrace, double time, double t0)
//    {
//    	double thisOffset = Math.abs(traceOffsets[nFirstGatherLineTrace + nGatherTrace]);
//    	double tMin = lineSet.getZMin();
//    	double tInc = lineSet.getZInc();
//    	int t0Index = (int)Math.floor(t0/tInc + tMin);
//    	if (t0Index > stretchMuteOffsets.length-1) return false;
//    	double stretchOffset = stretchMuteOffsets[t0Index];
//    	return isInsideMuteRange(nGatherTrace, time, thisOffset, stretchOffset);
//    }

    /**
     * returns false if data should be muted
     *
     * @param trace = the gather trace
     * @param t     = time of sample in question
     * @param t0    = zero-offset time
     * @return
     */
    public boolean isInsideMuteRange(StsGatherTrace trace, double t, double t0)
    {
        double offset = trace.x;
        if (muteRange == null) return isInsideStretchMuteRange(t, t0, offset);
        return t >= muteRange[trace.index][0] && t <= muteRange[trace.index][1] && isInsideStretchMuteRange(t, t0, offset);
    }

    /**
     * for display purposes, sometimes we need to check for mute solely based on trace index and time
     * <p/>
     * returns false if data should be muted
     *
     * @param nGatherTrace
     * @param t
     * @return
     */
    private boolean isInsideMuteRange(int nGatherTrace, double t)
    {
        if (muteRange == null) return true;
        return t >= muteRange[nGatherTrace][0];
    }

//    /**
//     * returns false if data should be muted
//     * 
//     * @param nTrace = Trace number within Gather
//     * @param stretchMuteOffsets = stretch mute offsets / one per time sample
//     * @return
//     */
//    public boolean isInsideStretchMuteRange(int nGatherTrace, double t0, double[] stretchMuteOffsets)
//    {
//    	if (stretchMuteOffsets == null || traceOffsets == null) return true;
//    	double thisOffset = Math.abs(traceOffsets[nFirstGatherLineTrace + nGatherTrace]);
//    	double tMin = lineSet.getZMin();
//    	double tInc = lineSet.getZInc();
//    	int t0Index = (int)Math.floor(t0/tInc + tMin);
//    	if (t0Index > stretchMuteOffsets.length-1) return false;
//    	double stretchOffset = stretchMuteOffsets[t0Index];
//    	return thisOffset <= stretchOffset;
//    }


    private float computeNmoMute(StsPoint mute, double offset)
    {
        float v = mute.v[0];
        float t0 = mute.v[1];
        return (float) Math.sqrt(t0 * t0 + offset * offset / (v * v));
    }

    public void drawTextureTileTimeSurface(StsTextureTile tile, GL gl)
    {
//		int nGatherTraces = nLastGatherTrace - nFirstGatherTrace + 1;
        if (nGatherTraces < tile.rowMin) return;

        if (traceOffsets == null) return;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;

        int rowMin = tile.rowMin;
        int rowMax = Math.min(nGatherTraces - 1, tile.rowMax);
        double rowTexCoor = tile.minRowTexCoor + dRowTexCoor * rowMin;

        double tMin = line.getZMin();
        double tInc = line.getZInc();
        double x1;
        int nTrace = rowMin + nFirstGatherLineTrace;
        x1 = traceOffsets[nTrace++];
//		x1 = line.traceOffsets[nTrace++];
        for (int row = rowMin + 1; row <= rowMax; row++, rowTexCoor += dRowTexCoor, nTrace++)
        {
            double x0 = x1;
            x1 = traceOffsets[nTrace];

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            double t = tMin + tile.colMin * tInc;

            for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
            {
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex2d(x0, t);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex2d(x1, t);
            }
            gl.glEnd();
        }
    }

    public boolean displayAttributeOnGather(StsGLPanel3d glPanel3d, String attName, float[][] axisRange, StsColor color, double shift)
    {
        double x = 0.0;

        boolean isNMOed = lineSet.getIsNMOed();
        boolean flatten = lineSet.lineSetClass.getFlatten();
        GL gl = glPanel3d.getGL();

        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        try
        {
            gl.glDisable(GL.GL_LIGHTING);

            if (line == null)
                return false;
            double[] attribute = line.getAttributeArray(attName, gatherRow, gatherCol);
            if (attribute == null)
            {
                StsMessageFiles.logMessage("Requested attribute: " + attName + " cannot be accessed for gather at (" + gatherRow + "," + gatherCol + ") from " + line.getName());
                return true;
            }

            // Check attribute range against axis range
            double min = Double.MAX_VALUE;
            double max = -Double.MAX_VALUE;
            for (int i = 0; i < attribute.length; i++)
            {
                if (min > attribute[i]) min = attribute[i];
                if (max < attribute[i]) max = attribute[i];
            }
            // Selected attribute is outside time axis range
            if (((min > axisRange[1][0]) && (max > axisRange[1][0])) ||
                ((min < axisRange[1][1]) && (max < axisRange[1][1])))
            {
                new StsMessage(glPanel3d.window, StsMessage.WARNING, "Entire attribute range outside axis range.\n Verify that " + attName + " is in time/depth units?");
                return false;
            }

            if (!isNMOed)
            {
                if (!flatten)
                {
                    int[] gatherIndexRange = this.getTraceDisplayRange(axisRange[0], offsetAxisType);
                    if (gatherIndexRange == null) return false;
                    int nFirstGatherIndex = gatherIndexRange[0];
                    int nLastGatherIndex = gatherIndexRange[1];
                    int nFirstLineIndex = nFirstGatherIndex + nFirstGatherLineTrace;

                    if (offsetAxisType != StsWiggleDisplayProperties.OFFSET_AXIS_ABS_INDEX)
                    {
                        if (traceOffsets == null) return false;

                        int attIdx = 0;
                        for (int nGatherIndex = nFirstGatherIndex, nLineIndex = nFirstLineIndex; nGatherIndex <= nLastGatherIndex; nGatherIndex++, nLineIndex++)
                        {
                            if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
                                x = traceOffsets[nLineIndex];
                            else
                                x = nGatherIndex;

                            if (isInsideMuteRange(nGatherIndex, attribute[attIdx]))
                                color.setGLColor(gl);
                            else
                                outsideMuteColor.setGLColor(gl);

                            // Selected attribute is outside time axis range
                            if ((attribute[attIdx] > axisRange[1][0]) || (attribute[attIdx] < axisRange[1][1]))
                                continue;

                            StsGLDraw.drawPoint2d((float) x, (float) attribute[attIdx] + (float) shift, StsColor.BLUE, gl, 4);
                            attIdx++;
                        }
                    }
                    else
                    {
                        if (absOffsetSortedLineIndices == null)
                            computeSortedTraceIndicesByAbsOffset();

                        if (traceOffsets == null) return false;

                        int idx = 0, attIdx = 0;

                        for (int i = 0; i < nGatherTraces; i++)
                        {
                            if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_ABS_INDEX)
                                idx = absOffsetSortedLineIndices[i];
                            else
                                idx = i + nFirstGatherLineTrace;

                            attIdx = idx - nFirstGatherLineTrace;

                            // Selected attribute is outside time axis range
                            if ((attribute[attIdx] > axisRange[1][0]) || (attribute[attIdx] < axisRange[1][1]))
                                continue;

                            x = traceOffsets[idx];

                            if (isInsideMuteRange(i, attribute[attIdx]))
                                color.setGLColor(gl);
                            else
                                outsideMuteColor.setGLColor(gl);

                            // Draw point
                            StsGLDraw.drawPoint2d((float) i, (float) attribute[attIdx] + (float) shift, StsColor.BLUE, gl, 4);
                        }
                    }
                }
            }
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch (Exception e)
        {
            StsException.outputException("Error displaying attribute on gather", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    protected void checkSetOrder()
    {
        byte order = lineSet.semblanceComputeProperties.order;
        checkSetOrder(order);
    }

    protected void checkSetOrder(byte order)
    {
        velocities.checkSetOrder(order);
        if (inputVelocities != null) inputVelocities.checkSetOrder(order);
    }


    public float getOffset(String offsetAxisType, float offsetValue)
    {
        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE_STRING)
            return offsetValue;
        int nTrace = (int) StsMath.minMax(offsetValue, 0, (float) nGatherTraces);
        if (traceOffsets == null) return 0.0f;
        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_ABS_INDEX_STRING)
        {
            if (absOffsetSortedLineIndices == null) computeSortedTraceIndicesByAbsOffset();
            return (float) Math.abs(traceOffsets[absOffsetSortedLineIndices[nTrace]]);
        }
        else
            return (float) traceOffsets[nTrace + nFirstGatherLineTrace];
    }

    protected boolean computeDrawOffsets()
    {
        if (nGatherTraces <= 0) return false;
        double[] traceOffsets = line.getTraceOffsets();
        if (traceOffsets == null) return false;
        int nValues = nGatherTraces * 2 - 1;
        drawOffsets = new double[nValues];
        int nTrace = nFirstGatherLineTrace;
        double x1 = traceOffsets[nTrace++];
        int n = 0;
        for (int i = 1; i < nGatherTraces; i++, nTrace++)
        {
            double x0 = x1;
            x1 = traceOffsets[nTrace];
            drawOffsets[n++] = x0;
            if (x1 > 0.0 && x0 < 0.0)
                drawOffsets[n++] = 0.0;
            else
                drawOffsets[n++] = (x0 + x1) / 2;
        }
        drawOffsets[n] = x1;
        return true;
    }

    public double[] getGatherOffsets()
    {
        if (nGatherTraces <= 0) return null;
        double[] gatherOffsets = new double[nGatherTraces];
        int nTrace = nFirstGatherLineTrace;
        double[] traceOffsets = line.getTraceOffsets();
        for (int i = 0; i < nGatherTraces; i++, nTrace++)
            gatherOffsets[i] = traceOffsets[nTrace];
        return gatherOffsets;
    }

//    /**
//     * compute stretch mute at every offset and midpoint between every offset. For two offsets straddling zero,
//     * compute the mute at zero instead of the midpoint. So the number of mute values is 2*nTraces-1.
//     * <p/>
//     * From Oz p.439, stretch factor s = dT/dt = T/t where T is offset time and t is tzero time.
//     * So s = SQRT(1 + x*x/(t*t*v*v).
//     * The stretch mute is a percent stretch with 0% as no stretch (which is at x=0).
//     * So stretch mute = (s - 1)*100.
//     * <p/>
//     * Given velocity is linear on an interval from t0 to t1: v = v0 + a(t-t0) where a = (t1 - t0)/(v1 - v0)
//     * Stretch factor equation can be rewritten:
//     * x = SQRT(s*s-1)*t*v = SQRT(s*s-1)*t*(v0 + a(t-t0))
//     * which is a quadratic equation in t given x, s, and a:
//     * a*t*t + (v0-a*t0)*t - x/SQRT(s*s-1) = 0
//     * <p/>
//     * We solve this equation for each offset x to find the t value for the stretch mute
//     *
//     * stretchMutePoints are at each offset and midpoint between offsets, so number of points is 2*nTraces-1.
//     * For a set of traces, the mute values are at 0, 2, 4...., so they are at: 2*(nTrace-1)
//     */
//    public double[] compute2ndOrderStretchMuteTimes(StsVelocityProfile velocityProfile, double vFactor)
//    {
////        if (stretchMuteTimes != null) return stretchMuteTimes;
////        if(!checkComputeDrawOffsets()) return null;
//        if (velocityProfile.getNProfilePoints() == 0)
//            return null;
//
//        if (drawOffsets == null && !computeDrawOffsets()) return null;
//
//        double[] stretchMuteTimes = new double[drawOffsets.length];
//        double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
//        double c = stretchMute / 100 + 1;
//        double coef = Math.sqrt(c * c - 1);
//        StsPoint[] profilePoints = velocityProfile.getProfilePoints();
//        int nPoints = profilePoints.length;
//        if (nPoints < 2) return null;
//
//        double v0 = profilePoints[0].v[0]*vFactor;
//        double tz0 = 0.0;
//        double z0 = 0.0;
//        double v1 = profilePoints[0].v[0]*vFactor;
//        double tz1 = profilePoints[0].v[1];
//        double z1 = v1 * tz1;
//        double a = 0.0;
//        double b = v0;
//        int p = 0;
//        boolean isFlattened = lineSet.lineSetClass.getFlatten();
//        for (int n = 0; n < drawOffsets.length; n++)
//        {
//        	double x =  Math.abs(drawOffsets[n]);
//            double z = x / coef;
//            while (z > z1 && p < nPoints - 1)
//            {
//                v0 = v1;
//                tz0 = tz1;
//           //   z0 = z1;
//                p++;
//                v1 = profilePoints[p].v[0]*vFactor;
//                tz1 = profilePoints[p].v[1];
//                z1 = v1 * tz1;
//                if (tz1 > tz0)
//                {
//                    a = (v1 - v0) / (tz1 - tz0);
//                    b = v0 - a * tz0;
//                }
//            }
//            double tz;
//            if (a != 0.0)
//            {
//                if(b < 0)
//                    StsException.systemDebug(this, "compute2ndOrderStretchMuteTimes", "b is negative.");
//                tz = (-b + Math.sqrt(b * b + 4 * a * z)) / (2 * a);
//            }
//            else
//                tz = z / b;
//
//            if (isFlattened)
//                stretchMuteTimes[n] = tz;
//            else
//                stretchMuteTimes[n] = tz*c; //* c;  multiplying by ( stretchMute / 100 + 1 ) doesn't make sense - output is supposed to be tzero anyway
//        }
//        return stretchMuteTimes;
//    }

    /**
     * avoid double solution to quadratic equation by solving for mute offsets instead of mute times
     * <p/>
     * From Oz p.439, stretch = (t - t0) / t0 where t0 equals zero-offset time
     * So stretch = SQRT(1 + x*x/(t0*t0*v*v) - 1.
     * solve for x = t0*v*sqrt((stretch + 1)^2 - 1))
     * stretch factor = stretch + 1 = percent_stretch/100 + 1
     * x = t0 * v * sqrt(SFactor*SFactor - 1)
     * <p/>
     * returns x values for each t0 time sample in the data
     */
    private double[] compute2ndOrderStretchMuteOffsets(StsVelocityProfile velocityProfile, double vFactor)
    {
        if (velocityProfile == null || velocityProfile.getNProfilePoints() == 0)
            return null;

        if (vFactor == 0) vFactor = 1;

        int nTZeroes = lineSet.nSlices;
        float t0 = lineSet.getZMin();
        float tInc = lineSet.getZInc();
        double[] stretchMuteOffsets = new double[nTZeroes];
        setStretchMuteFactor();
        setStretchMuteMinOffset();
        double coef = Math.sqrt(stretchFactor * stretchFactor - 1);

        for (int n = 0; n < nTZeroes; n++, t0 += tInc)
        {
            double v = velocityProfile.getVelocityFromTzero(t0);
            double x = t0 * v * coef;
            if (x < stretchMuteMinOffset) x = stretchMuteMinOffset;
            stretchMuteOffsets[n] = x;
        }
        return stretchMuteOffsets;
    }

//    public double[] computeFocusStretchMuteTimes(StsVelocityProfile velocityProfile, double vFactor)
//    {
//        if (velocityProfile.getNProfilePoints() == 0)
//            return null;
//
//        if (drawOffsets == null && !computeDrawOffsets()) return null;
//
//        double[] stretchMuteTimes = new double[drawOffsets.length];
//        double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
//        double sFactor = Math.sqrt((stretchMute / 100 + 1)*(stretchMute / 100 + 1)-1);
//        StsPoint[] profilePoints = velocityProfile.getProfilePoints();
//        int nPoints = profilePoints.length;
//        if (nPoints < 2) return null;
//        double[] muteOffsets = new double[nPoints];
//        double[] muteTimes =  new double[nPoints];
//        
//        //generate mute offset/time pairs
//        for (int n = 0; n < nPoints; n++)
//        {
//        	double t0 = profilePoints[n].v[1]; //time is same
//        	double vel = profilePoints[n].v[0]*vFactor; //time is same
//        	double muteOffset = (t0*vel*sFactor);
//            muteTimes[n] = t0;
//            muteOffsets[n] = muteOffset;
//        }
//        
//        //interpolate generated offset/time pairs to find times for actual offsets
//        //assumes offsets are increasing!!
//        double minTableMuteOffset = muteOffsets[0];
//        double maxTableMuteOffset = muteOffsets[muteOffsets.length-1];
//        int tableIndex=0;
//        for (int i=0; i < stretchMuteTimes.length; i++) {
//        	double absTraceOffset = Math.abs(drawOffsets[i]);
//        	if (absTraceOffset == 0) {
//        		stretchMuteTimes[i] = 0;
//        		continue;
//        	}
//        	
//        	//find index in table of offset just before absTraceOffset
//        	if (absTraceOffset < minTableMuteOffset) { //test for before first table offset
//        		tableIndex = 0;
//        	} else if (absTraceOffset >= minTableMuteOffset && absTraceOffset <= maxTableMuteOffset) { //test for inside table offsets
//        		for (int n= 0; n < muteOffsets.length-1; n++) {
//            		if ( absTraceOffset >= muteOffsets[n] && absTraceOffset <= muteOffsets[n+1]) {
//            			tableIndex = n;
//                		continue;
//            		}
//        		}
//        	} else { //test for after table offsets
//        		tableIndex = muteOffsets.length-2;
//        	}
//        	double tableMuteOffset1 = muteOffsets[tableIndex];
//        	double tableMuteOffset2 = muteOffsets[tableIndex+1];
//        	double muteT1 = muteTimes[tableIndex];
//        	double muteT2 = muteTimes[tableIndex+1];
//        	double offsetDiff = absTraceOffset - tableMuteOffset1;
//        	//calculate mute!!
//        	double muteTime = muteT1 + (offsetDiff)*(muteT2 - muteT1)/(tableMuteOffset2-tableMuteOffset1);
//        	stretchMuteTimes[i] = muteTime;
//        	continue;
//        }
//        return stretchMuteTimes;
//    }

//    public double[] computeConstantVelocityStretchMuteTimes(double v)
//    {
//
//        double stretchMute = lineSet.getWiggleDisplayProperties().stretchMute;
//        int nTrace = nFirstGatherLineTrace;
//        int nTraces = nLastGatherLineTrace - nFirstGatherLineTrace + 1;
//        double[] stretchMuteTimes = new double[nTraces];
//        double c = stretchMute / 100 + 1;
//        double f;
//        if(lineSet.isNMOed)
//            f = c/(Math.sqrt(c*c-1)*v);
//        else
//            f = 1/(Math.sqrt(c*c-1)*v);
//        for (int n = 0; n < nTraces; n++, nTrace++)
//            stretchMuteTimes[n] = f*Math.abs(traceOffsets[nTrace]);
//        return stretchMuteTimes;
//    }

    public boolean checkComputeStretchMuteOffsets()
    {
        if (stretchMuteOffsets == null)
        {
            compute2ndOrderStretchMuteOffsets();
        }
        if (stretchMuteOffsets == null) return false;
        return true;
    }

    private void compute2ndOrderStretchMuteOffsets()
    {
        stretchMuteOffsets = compute2ndOrderStretchMuteOffsets(velocityProfile, 1);
    }

//	public double[] computeStretchMutesWithVelocityFactor(StsVelocityProfile velocityProfile, double vFactor)
//    {
//        /*
//        if(vFactor == 1.0)
//        {
//            if(!checkComputeStretchMuteTimes()) return null;
//            return this.stretchMuteTimes;
//        }
//        else // velocity profile is multiplied by vFactor for residual semblance or VVS
//        */
//        {
//            double[] stretchMuteOffsets = compute2ndOrderStretchMuteOffsets(velocityProfile, vFactor);
//        	//double[] stretchMuteTimes = computeFocusStretchMuteTimes(velocityProfile, vFactor);
//            if(stretchMuteOffsets == null) muteRange = null;
//            return stretchMuteOffsets;
//        }
//    }

    public float getOffsetForTypeValue(float offsetValue)
    {
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
            return offsetValue;
        // return Math.abs(offsetValue);

        offsetValue = StsMath.minMax(offsetValue, 0.0f, superGather.nSuperGatherTraces - 1);
        int index = (int) offsetValue;
        return (float) superGather.getSortedTraceOffset(index);
    }

    /**
     * Given the offsetValue, the units of which depends on type, return the interpolating index between traces
     * which bounds this value.  If type is VALUE, the offsetValue is offset; if it is INDEX, it's already the
     * interpolating index so just return; otherwise type is ABS_INDEX, traces must be resorted into absolute offsets,
     * so return the trace index whose offset is just below this one plus a fraction which is fractional distance
     * to next trace in the absolute value order.
     *
     * @param offsetValue
     * @return
     */
    public float getTraceIndexFloatFromTypeValue(float offsetValue)
    {
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
            return StsMath.arrayIndexF((double) offsetValue, traceOffsets, nFirstGatherLineTrace, nLastGatherLineTrace);

        offsetValue = StsMath.minMax(offsetValue, 0.0f, nGatherTraces - 1);

        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_INDEX)
            return offsetValue;

        int index = (int) offsetValue;
        double f = offsetValue - index;

        if (traceOffsets == null) return 0.0f;

        if (absOffsetSortedLineIndices == null) computeSortedTraceIndicesByAbsOffset();
        return (float) (absOffsetSortedLineIndices[index] - nFirstGatherLineTrace + f);
    }

    public float getTraceIndexFromOffset(float offsetValue)
    {
        double indexF;
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
            return StsMath.arrayIndexF((double) offsetValue, traceOffsets, nFirstGatherLineTrace, nLastGatherLineTrace);
        else
            return StsMath.minMax(offsetValue, 0.0f, nGatherTraces - 1);
    }

    public double[] getTraceOrderedAttributes(StsMappedDoubleBuffer lineAttributesArrayBuffer)
    {
        int nAttributes = lineSet.nAttributes;
        int nLineTraces = line.nLineTraces;
        double[] attributes = new double[nGatherTraces * nAttributes];
        int position = nFirstGatherLineTrace;
        double[] gatherAttributeValues = new double[nGatherTraces];
        for (int a = 0; a < nAttributes; a++, position += nLineTraces)
        {
            lineAttributesArrayBuffer.position(position);
            lineAttributesArrayBuffer.get(gatherAttributeValues);
            int i = a;
            for (int t = 0; t < nGatherTraces; t++, i += nAttributes)
                attributes[i] = gatherAttributeValues[t];

        }
        return attributes;
    }

    public void stretchMuteChanged()
    {
        stretchMuteOffsets = null;
        setStretchMuteFactor();
        setStretchMuteMinOffset();
        this.checkComputeStretchMuteOffsets();
    }

//	public boolean isInsideMuteRange(float t0, double[] stretchMuteTimes2, int n) {
//		return t0 >= stretchMuteTimes2[n];
//	}

    /**
     * yilmaz 440
     * <p/>
     * stretch factor = tOffset/tZero
     * <p/>
     * based on most simple definition of stretch
     * <p/>
     * doesn't depend on what method is used to calculate NMO
     * <p/>
     * sets t0 to 1st sample if t0 = 0
     *
     * @param t      = time at offset
     * @param t0
     * @param offset
     * @return false if data needs to be muted
     */
    public boolean isInsideStretchMuteRange(double t, double t0, double offset)
    {
        offset = Math.abs(offset);
        if (offset < stretchMuteMinOffset)
        {
            return true;
        }
        if (t0 == 0) t0 = lineSet.getZInc();
        double stretch = t / t0;
        return stretch < stretchFactor;
    }

//	public double[] computeConstantVelocityStretchMuteOffsets(double v) {
//		
//        int nTZeroes = lineSet.nCroppedSlices;
//        float t0 = lineSet.getZMin();
//        float tInc = lineSet.getZInc();
//        double[] stretchMuteOffsets = new double[nTZeroes];
//        double percentStretch = lineSet.getWiggleDisplayProperties().stretchMute;
//        double sFactor = percentStretch / 100 + 1;
//        double coef = Math.sqrt(sFactor * sFactor - 1);
//        
//        for (int n=0; n < nTZeroes; n++, t0 += tInc)
//        {
//        	double x = t0 * v * coef;
//        	stretchMuteOffsets[n] = x;
//        }
//        return stretchMuteOffsets;
//	}
}
