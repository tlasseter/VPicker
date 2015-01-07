package com.Sts.Types.PreStack;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import javax.media.opengl.GL;

import com.Sts.Actions.Wizards.VelocityAnalysis.StsVelocityAnalysisEdit;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.Seismic.StsTraceUtilities;
import com.Sts.Types.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

/**
 * A superGather consists of a single gather or a group of gathers in a rectangle or cross configuration.
 * Row and col are the indices of this superGather.  Computations are done on a gather-by-gather basis,
 * so gatherRow and gatherCol are the indices of the current gather being computed.  gatherData is the
 * trace data for the current gather at gatherRow & gatherCol.
 */
public class StsSuperGather
{
    /** model singleton */
    StsModel model;
    /** gathers in superGather.  See GatherSet innerClass for details. */
    public StsGather[] gathers;
    /** center gather of superGather */
    public StsGather centerGather;
    /** type of superGather: rectangle, cross, or single */
    byte gatherType = StsPreStackLineSetClass.SUPER_SINGLE;
    /** number of gather cols in superGather rectangle, cross, or single pattern */
    public int nSuperGatherCols = -1;
    /** number of gather rows in superGather rectangle, cross, or single pattern */
    public int nSuperGatherRows = -1;
    /** number of gathers in superGather (nRows*nCols for rectangle, nRows+nCols-1 for cross, 1 for single */
    public int nGathersInSuperGather = -1;
    /** set which owns this gather */
    public StsPreStackLineSet lineSet;
    /** StsClass managing lineSets; used for getting various properties */
    public StsPreStackLineSetClass lineSetClass;
    /** line this superGather is on; row index on grid for 3d, line number for 2d */
    public int superGatherRow = -1;
    /** sequence number for this superGather; col index on grid for 3d, sequence number on line for 2d */
    public int superGatherCol = -1;
    /** wiggle properties used by this superGather */
    public StsWiggleDisplayProperties wiggleProperties;
    /** indicates this line has been NMOed */
    //    public boolean isNMOed = false;
    /** semblance data converted to bytes for display */
    public byte[] semblanceBytes;
    /** semblance backbone data converted to bytes for display */
    public byte[] residualSemblance;
    /** constant velocity stack traces data to be displayed */
    public float[][][] cvsTraces;
    /** variable velocity stack traces data to be displayed */
    public float[][][] vvsTraces;
    /** current profile points used for semblanceBytes backbone */
    //    public StsPoint[] profilePoints;
    /** fractional max error allowed in curve fitting wavelets */
    double maxError = 0.01f;
    /** times for semblance array (used in enhanced semblance calculations and display) */
    public double[][] semblanceTimes;
    //    public double[] semblanceVelocities;
    /** estimated maxSemblanceCurve based on stat analysis */
    StsPoint[] maxSemblanceCurve;//	double[] traceOffsets = null;
    /** volume min: line2d.zMin */
    public float seismicVolumeZMin;
    /** volume max: line2d.zMax */
    public float seismicVolumeZMax;
    /** volume zInc: line2d.zInc */
    public float seismicVolumeZInc;
    /** number of samples in volume or line (2d) */
    public int nVolumeSlices;
    /** time/depth of first semblanceBytes point; used by backbone plot where range is inside defined velocity profile range */
    public double semblanceZMin;
    /** time/depth of last semblanceBytes point; used by backbone plot where range is inside defined velocity profile range */
    public double semblanceZMax;
    /** time/depth increment of semblanceBytes points */
    public double semblanceZInc;
    /** index of first semblanceBytes point in vrms and vint arrays */
    public int semblanceZMinIndex;
    /** index of last semblanceBytes point in vrms and vint arrays */
    public int semblanceZMaxIndex;//	public int firstIntVelIndex = 0;
    /** Number of vertical samples on semblanceBytes; computed from range and increment information */
    public int nSemblanceSamples;
    /** temporarily use second-order (as when dragging NMO curve on gather display) */
    public boolean useSecondOrder;
    /** Number of resampled points for semblance computed by cubic curve-fitting */
    public int nEnhancedSemblancePoints = 0;//	public byte order;
    /** glPanel of parent of family this gather belongs to */
    public StsWin3dBase parentWindow = null;
    /** description used for title of gather display typically */
    public String gatherDescription = "null";
    /** number of total traces in superGather */
    public int nSuperGatherTraces = 0;
    /** lock object while computing velocity stack (CVS or VVS) */
    public boolean computingStacks = false;
    /** semblance data is being computed */
    public boolean computingSemblance = false;
    /** indicates files needed for this superGather are available */
    public boolean filesOk = false;
    /** current velocityProfile */
    public StsVelocityProfile velocityProfile = null;
    /** input velocityProfile (if available) */
    public StsVelocityProfile inputVelocityProfile = null;
    /** order of tzero calculation */
    public byte order = ORDER_2ND;
    /** datum shift to be applied */
    public double datumShift;
    /** maxAmplitude for this superGather used in normalizing CVS */
    public double superGatherMaxAmplitude;
    /** RMS normalized amplitude */
    public double rmsAmplitude;
    //    static private ComputeSemblanceStandardProcess semblanceStandardProcess = null;
    //    static private ComputeSemblanceEnhancedProcess semblanceEnhancedProcess = null;
    //    static private ComputeResidualSemblanceStandardProcess residualSemblanceStandardProcess = null;
    //    static private ComputeResidualSemblanceEnhancedProcess residualSemblanceEnhancedProcess = null;

    public boolean blockReprocessSemblance = false;
    /** Thread running CVS Panel construction if running */
    private Thread cvsThread;
    /** Thread running VVS Panel construction if running */
    private Thread vvsThread;

    public String errorMessage = "";
    /** interpolated points for each trace in the super gather */
    public StsSuperGatherTraceSet superGatherTraceSet;

    public float zMax;
    private double stretchFactor;
    private float stretchMuteMinOffset;
    /** currently loaded attribute array from all traces within supergather */
    transient private double[] attributes;
    /** name of currently loaded trace attribute */
    private String attributeName = StsWiggleDisplayProperties.ATTRIBUTE_NONE;

    /** approximate number of interpolated intervals from zInc to subSampled zInc */
    static final int nApproxSubIntervals = 5;
    // lineSets already have a user-selected minimum number of traces (lineSet.traceThreshold); use it instead.
    // TJL 4/23/09
    // public final static int minLiveTraces = 2;

    static final float colorSaturationFactor = 2.0f;

    static final boolean debug = false;
    static final boolean runTimer = false;
    static StsTimer timer = null;

    static final boolean debugVelTime = false;
    static final float debugVelocity = 6.5262f;
    static final float debugTime = 940.f;
    static int debugVelocityIndex;
    static int debugTimeIndex;
    static int debugSuperGatherRow;
    static int debugSuperGatherCol;
    static private boolean debugSemblance = false;

    static
    {
        if (runTimer) timer = new StsTimer();
    }

    public final static byte nullByte = StsParameters.nullByte;
    final static double roundOff = StsParameters.roundOff;

    /** minimum number of values to use when calculating coherence */
    static final int MIN_COHERENCE_VALS = 2;

    public final static byte POINT_ORIGINAL = 0;
    public final static byte POINT_MAXIMUM = 1;
    public final static byte POINT_MINIMUM = 2;
    public final static byte POINT_ZERO_CROSSING = 3;
    public final static byte POINT_FLAT_ZERO = 4;

    public final static int TOP_MUTE_INDEX = -1;
    public final static int BOT_MUTE_INDEX = -2;
    public final static int STRETCH_MUTE_INDEX = -3;
    public final static int NONE_INDEX = -99;

    public final static byte ORDER_2ND = StsSemblanceComputeProperties.ORDER_2ND;
    public final static byte ORDER_4TH = StsSemblanceComputeProperties.ORDER_4TH;
    public final static byte ORDER_6TH = StsSemblanceComputeProperties.ORDER_6TH;
    public final static byte ORDER_6TH_OPT = StsSemblanceComputeProperties.ORDER_6TH_OPT;

    static StsColor outsideMuteColor = StsColor.GRAY;
    //    transient StsColor muteColor = StsColor.GREEN;
    public static StsTimer StackTotalTimer;
    public static StsTimer StackingTimer;
    public static StsTimer LoadGathersTimer;
    public static StsTimer ApplyNMOTimer;
    public static StsTimer CalcVelocitiesTimer;
    public static StsTimer CalcTxTimer;
    public static StsTimer InterpolateTracesTimer;
    private static boolean debugTimer = false;
    public static StsTimer CalcMeanTimer;
    public static boolean debugStackTimer = false;

    public StsSuperGather()
    {

    }

    private StsSuperGather(StsModel model, StsPreStackLineSet lineSet, StsWin3dBase parentWindow)
    {
        this.model = model;
        this.lineSet = lineSet;
        lineSetClass = (StsPreStackLineSetClass) lineSet.getStsClass();
        this.parentWindow = parentWindow;
        wiggleProperties = lineSet.getWiggleDisplayProperties();
    }

    public StsSuperGather(StsModel model, StsPreStackLineSet lineSet)
    {
        this(model, lineSet, null);
    }

    public String getTitleString()
    {
        return getGatherDescription();
    }

    /** constructTraceAnalyzer */
    static public StsSuperGather constructor(StsModel model, StsPreStackLineSet lineSet, StsWin3dBase parentWindow)
    {
        try
        {
            StsSuperGather superGather = new StsSuperGather(model, lineSet, parentWindow);
            StsCursor3d cursor3d = parentWindow.getCursor3d();
            float x = cursor3d.getCurrentDirCoordinate(StsCursor3d.YDIR);
            float y = cursor3d.getCurrentDirCoordinate(StsCursor3d.XDIR);
            int[] rowCol = lineSet.getRowColFromCoors(x, y);
            int row = rowCol[0];
            int col = rowCol[1];
            superGather.initializeSuperGather(row, col);
            return superGather;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSuperGather.constructTraceAnalyzer(lineSet, glPanel3d) failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    static public StsSuperGather constructor(StsModel model, StsPreStackLineSet lineSet)
    {
        try
        {
            StsSuperGather superGather = new StsSuperGather(model, lineSet);
            return superGather;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSuperGather.constructTraceAnalyzer(lineSet) failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void setSemblanceRange()
    {
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        if (rangeProperties == null) return;
        semblanceZInc = rangeProperties.zInc;
        semblanceZMin = StsMath.intervalRoundUp(rangeProperties.zMin, seismicVolumeZMin, rangeProperties.zInc);
        semblanceZMax = StsMath.intervalRoundDown(rangeProperties.zMax, seismicVolumeZMin, rangeProperties.zInc);
        nSemblanceSamples = (int) Math.round((semblanceZMax - semblanceZMin) / semblanceZInc) + 1;
        semblanceZMinIndex = (int) Math.round((semblanceZMin - seismicVolumeZMin) / semblanceZInc);
        semblanceZMaxIndex = (int) Math.round((semblanceZMax - seismicVolumeZMin) / semblanceZInc);
    }

    public synchronized boolean computeSemblanceProcess(StsWin3dBase win3d)
    {
        if (blockReprocessSemblance)
        {
            if (debug) debugThread("computeSemblanceProcess running");
            return false;
        }
        blockReprocessSemblance = true;
        if (nSuperGatherTraces == 0) return false;

        byte semblanceType = lineSet.getSemblanceComputeProperties().semblanceType;
        ComputeSemblanceProcess computeSemblanceProcess = null;
        if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_STANDARD)
            computeSemblanceProcess = new ComputeSemblanceStandardProcess();
        else if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_STACKED_AMP)
            computeSemblanceProcess = new ComputeStackedAmplitudeProcess();
        else if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_ENHANCED)
            computeSemblanceProcess = new ComputeSemblanceEnhancedProcess();
        if (runTimer)
            timer.stopPrintReset("computeSemblance for super gather of semblance type " + semblanceType + " at row " + superGatherRow + " col " + superGatherCol);

        StsToolkit.runRunnable(computeSemblanceProcess);
        return true;
    }

    // The general ignore super gather is in superGatherProperties however, each view type may be
    // allowed to ignore it just for that type. (CVSView Only so far)
    //
    // Because we don't know who is ignoring and who isn't, we will setup superGather as if not ignored.
    // Actual data fetching and computing is done on an as needed basis, so there is little overhead/penalty for this.
    // TJL 4/10/07
    public void checkInitializeSuperGatherGeometry()
    {
        if (superGatherGeometryChanged())
            initializeGathers();
    }

    public boolean superGatherGeometryChanged()
    {
        boolean changed = (gathers == null);
        StsSuperGatherProperties superGatherProperties = lineSet.superGatherProperties;
        byte currentGatherType = superGatherProperties.getGatherType();
        if (gatherType != currentGatherType)
        {
            changed = true;
            gatherType = currentGatherType;
        }
        int currentNSuperGatherCols = superGatherProperties.getNSuperGatherCols();
        if (nSuperGatherCols != currentNSuperGatherCols)
        {
            changed = true;
            nSuperGatherCols = currentNSuperGatherCols;
        }
        int currentNSuperGatherRows = superGatherProperties.getNSuperGatherRows();
        if (nSuperGatherRows != currentNSuperGatherRows)
        {
            changed = true;
            nSuperGatherRows = currentNSuperGatherRows;
        }
        if (changed) initializeGeometry();
        return changed;
    }

    public void initializeGeometry()
    {
        StsSuperGatherProperties superGatherProperties = lineSet.superGatherProperties;
        nSuperGatherRows = superGatherProperties.getNSuperGatherRows();
        nSuperGatherCols = superGatherProperties.getNSuperGatherCols();
        gatherType = superGatherProperties.getGatherType();
        switch (gatherType)
        {
            case StsPreStackLineSetClass.SUPER_SINGLE:
                nGathersInSuperGather = 1;
                break;
            case StsPreStackLineSetClass.SUPER_CROSS:
                nGathersInSuperGather = nSuperGatherCols + nSuperGatherRows - 1;
                break;
            case StsPreStackLineSetClass.SUPER_RECT:
                nGathersInSuperGather = nSuperGatherCols * nSuperGatherRows;
                break;
            case StsPreStackLineSetClass.SUPER_INLINE:
                nGathersInSuperGather = nSuperGatherCols;
                break;
            case StsPreStackLineSetClass.SUPER_XLINE:
                nGathersInSuperGather = nSuperGatherRows;
                break;
        }
    }

    public int[] getSuperGatherRowColStart()
    {
        return new int[]{superGatherRow - nSuperGatherRows / 2, superGatherCol - nSuperGatherCols / 2};
    }

    private void debugThread(String message)
    {
        System.out.println(Thread.currentThread().getName() + " class " + StsToolkit.getSimpleClassname(this) + " " + message +
            " computingSemblance " + computingSemblance + " superGatherRow " + superGatherRow + " superGatherCol " + superGatherCol);
    }

    public synchronized boolean computeResidualSemblanceProcess(StsWin3dBase win3d)
    {
        if (nSuperGatherTraces == 0) return false;

        if (blockReprocessSemblance)
        {
            if (debug) debugThread("computeResidualSemblanceProcess running");
            return false;
        }
        if (velocityProfile == null || velocityProfile.getNProfilePoints() < 2) return false;
        byte semblanceType = lineSet.getSemblanceComputeProperties().semblanceType;
        ComputeSemblanceProcess computeSemblanceProcess = null;
        if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_STANDARD)
            computeSemblanceProcess = new ComputeResidualSemblanceStandardProcess();
            //else if(semblanceType == StsSemblanceComputeProperties.SEMBLANCE_ENHANCED)
            //    computeSemblanceProcess = new ComputeResidualSemblanceEnhancedProcess();
        else if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_STACKED_AMP)
            //computeSemblanceProcess = new ComputeResidualSemblanceStackedAmpProcess();
            System.out.println("ComputeResidualSemblanceStackedAmpProcess not implemented yet ");
        StsToolkit.runRunnable(computeSemblanceProcess);
        //backboneExecutionPipe.addCancelInProcess(computeProcess, this);
        return true;
    }

    protected class Semblance
    {
        protected double[] sumValue;
        protected double[] sumValueSq;
        protected int[] nVals;
        protected int nValues;

        protected Semblance(int nValues)
        {
            this.nValues = nValues;
            initialize();
        }

        protected void addSemblanceValue(int n, double v)
        {
            nVals[n]++;
            sumValue[n] += v;
            sumValueSq[n] += v * v;
        }

        public final float[] getSemblanceValues()
        {
            float[] semblanceValues = new float[nValues];
            for (int n = 0; n < nValues; n++)
            {
                if (nVals[n] < 4)
                    semblanceValues[n] = 0.0f;
                else
                    semblanceValues[n] = (float) (sumValue[n] * sumValue[n] / sumValueSq[n] / nVals[n]);
            }
            return semblanceValues;
        }

        public final double getSemblanceValue(int n)
        {
            if (nVals[n] < 4)
                    return 0.0;
                else
                    return (sumValue[n] * sumValue[n] / sumValueSq[n] / nVals[n]);
        }

        public void initialize()
        {
            sumValue = new double[nValues];
            sumValueSq = new double[nValues];
            nVals = new int[nValues];
        }

        protected void addStackValue(int i, double v)
        {
            nVals[i]++;
            sumValue[i] += v;
        }

        public float[] getStackValues()
        {
            float[] stackValues = new float[nValues];
            for (int i = 0; i < nValues; i++)
            {
                if (nVals[i] < MIN_COHERENCE_VALS)
                    stackValues[i] = 0.0f;
                else
                    stackValues[i] = (float) (sumValue[i] / nVals[i]);
            }
            return stackValues;
        }
    }

    public boolean scaleSemblanceTrace(float[] floatData, byte[] byteData)
    {
        if (floatData == null) return false;
        for (int n = 0; n < floatData.length; n++)
            byteData[n] = StsMath.unsignedIntToUnsignedByte254((int) (254 * floatData[n]));
        return true;
    }

    public synchronized void resetReprocessSemblanceBlockFlag()
    {
        computingSemblance = false;
        blockReprocessSemblance = false;
    }

    public synchronized void resetReprocessResidualBlockFlag()
    {
        computingSemblance = false;
        //System.out.println("resetReprocessBackboneBlockFlag");
        blockReprocessSemblance = false;
        //blockReprocessBackbone = false;
    }

    /*public synchronized void resetReprocessCVSBlockFlag()
    {
        blockReprocessCVS = false;
    }*/

    public void createNullSemblanceBytes()
    {
        semblanceBytes = new byte[4];
        for (int n = 0; n < 4; n++)
            semblanceBytes[n] = StsParameters.nullByte;
    }

    public int getNumberGathers()
    {
        return nGathersInSuperGather;
    }

    public boolean isInsideRowColRange(int row, int col)
    {
        return lineSet.isInsideRowColRange(row, col);
    }

    public float[][] computeSuperGatherCVSTraces(int nDisplaySlices, float velocity)
    {
        checkInitializeSuperGatherGeometry();
        float[][] cvsTraces = new float[nGathersInSuperGather][nDisplaySlices];
        superGatherMaxAmplitude = 0.0;
        for (int n = 0; n < nGathersInSuperGather; n++)
        {
            cvsTraces[n] = gathers[n].computeConstantVelocityStackTrace(nDisplaySlices, velocity);
            superGatherMaxAmplitude = Math.max(superGatherMaxAmplitude, gathers[n].maxAmplitude);
        }
        return cvsTraces;
    }

    public float[][] computeCVSTraces(int nDisplaySlices, float velocity)
    {
        checkInitializeSuperGatherGeometry();
        float[][] cvsTraces = new float[nGathersInSuperGather][nDisplaySlices];
        superGatherMaxAmplitude = 0.0;
        for (int n = 0; n < nGathersInSuperGather; n++)
        {
            cvsTraces[n] = gathers[n].computeConstantVelocityStackTrace(nDisplaySlices, velocity);
            superGatherMaxAmplitude = Math.max(superGatherMaxAmplitude, gathers[n].maxAmplitude);
        }
        return cvsTraces;
    }

    public byte getGatherType()
    {
        return gatherType;
    }

    /**
     * The single gather or superGather center row or col have changed.
     * Initialize the new gather/superGather.
     * Returns true if already initialized or if it succesfully reinitialized; otherwise returns false.
     */
    synchronized public boolean initializeSuperGather(int row, int col)
    {
        if (debug)
            debugThread("StsSuperGather.initializeSuperGather(row,col) called with row " + row + " col " + col + ".");
        boolean initialized = true;
        if (row != -1)
        {
            if (superGatherRow != row)
            {
                if (row >= lineSet.nRows) return false;
                superGatherRow = row;
                initialized = false;
            }
        }
        if (col != -1)
        {
            if (superGatherCol != col)
            {
                if (col >= lineSet.nCols) return false;
                superGatherCol = col;
                initialized = false;
            }
        }

        if (initialized) return true;
        return reinitializeSuperGather();
    }

    public void reinitialize()
    {
        initializeSuperGather();
    }

    /** called when we want to initialize the superGather location, but not the velocity.  This is used by the VVS panel. */
    private boolean initializeGeometry(int row, int col)
    {
        if (row == superGatherRow && col == superGatherCol) return true;
        if (row < 0 || row >= lineSet.nRows) return false;
        if (col < 0 || col >= lineSet.nCols) return false;
        filesOk = false;
        nSuperGatherTraces = 0;
        gathers = null;
        superGatherTraceSet = null;

        if (debug) debugThread("StsSuperGather.initializeSuperGather() called.");
        superGatherRow = row;
        superGatherCol = col;
        checkRowColLimits();

        if (lineSet == null) return false;
        if (!lineSet.hasGatherTraces(superGatherRow, superGatherCol)) return false;
        StsPreStackLine line = lineSet.getDataLine(superGatherRow, superGatherCol);
        if (line == null) return false;
        filesOk = line.checkFiles();
        if (!filesOk) return false;

        seismicVolumeZMin = line.zMin;
        seismicVolumeZMax = line.zMax;
        seismicVolumeZInc = line.zInc;
        nVolumeSlices = line.nSlices;
        initializeGeometry();
        initializeGathers();
        checkSetDatumShift();

        return true;
    }

    /** clear data associated with the superGather and return true if superGather is ok for this new row and col */
    private boolean reinitializeSuperGather()
    {
        clearData();
        return initializeSuperGather();
    }

    private void clearData()
    {
        semblanceBytes = null;
        residualSemblance = null;
        vvsTraces = null;
        cvsTraces = null;
    }

    /** clear data associated with the superGather and return true if superGather is ok for this new row and col */
    private boolean initializeSuperGather()
    {
        filesOk = false;
        nSuperGatherTraces = 0;
        gathers = null;
        superGatherTraceSet = null;

        if (debug) debugThread("StsSuperGather.initializeSuperGather() called.");
        if (superGatherRow == -1 || superGatherCol == -1) return false;
        checkRowColLimits();

        if (lineSet == null) return false;
        if (!lineSet.hasGatherTraces(superGatherRow, superGatherCol)) return false;
        StsPreStackLine line = lineSet.getDataLine(superGatherRow, superGatherCol);
        if (line == null) return false;
        filesOk = line.checkFiles();
        if (!filesOk) return false;

        seismicVolumeZMin = line.zMin;
        seismicVolumeZMax = line.zMax;
        seismicVolumeZInc = line.zInc;
        nVolumeSlices = line.nSlices;

        //        isNMOed = line.isNMOed;
        initializeCurrentVelocityProfile();
        // StsVelocityProfile.debugPrintProfile(this, "initializeSuperGather(after)", velocityProfile);

        setSemblanceRange();
        initializeGeometry();
        initializeGathers();
        setGatherDescription(); //needs to be here so that it can see how many traces!
        checkSetDatumShift();
        setStretchMuteFactor();
        setStretchMuteMinOffset();
        // StsVelocityProfile.debugPrintProfile(this, "initializeSuperGather(end)", velocityProfile);

        return true;
    }

    private void checkRowColLimits()
    {
        int[] rowCol = lineSet.adjustLimitRowCol(superGatherRow, superGatherCol);
        superGatherRow = rowCol[0];

        superGatherCol = rowCol[1];
    }

    private void initializeCurrentVelocityProfile()
    {
        if (lineSet.velocityModel == null) return;
        //        velocityProfile = lineSet.velocityModel.getVelocityProfile(superGatherRow, superGatherCol);
        velocityProfile = lineSet.velocityModel.getComputeVelocityProfile(superGatherRow, superGatherCol);
        if (lineSet.inputVelocityModel != null)
            inputVelocityProfile = lineSet.inputVelocityModel.getComputeVelocityProfile(superGatherRow, superGatherCol);
    }


    private void muteRangeChanged()
    {
        for (int n = 0; n < gathers.length; n++)
            gathers[n].muteRangeChanged();
    }

    /*
private void checkComputeMuteRange()
        {
            for(int n = 0; n < gathers.length; n++)
                gathers[n].checkComputeMuteRange(velocityProfile);
        }
    */
    private void computeDrawOffsets()
    {
        centerGather.computeDrawOffsets();
    }

    public void clearGathers()
    {
        gathers = null;
    }

    private void initializeGathers()
    {
        if (gathers == null || gathers.length != nGathersInSuperGather)
        {
            gathers = new StsGather[nGathersInSuperGather];
            for (int n = 0; n < nGathersInSuperGather; n++)
                gathers[n] = new StsGather(this);
        }
        int startRow, startCol, i;
        switch (gatherType)
        {
            case StsPreStackLineSetClass.SUPER_SINGLE:
                initializeGather(gathers[0], superGatherRow, superGatherCol);
                break;
            case StsPreStackLineSetClass.SUPER_CROSS:
                startRow = superGatherRow - nSuperGatherRows / 2;
                startCol = superGatherCol - nSuperGatherCols / 2;
                i = 0;
                for (int n = 0, col = startCol; n < nSuperGatherCols; n++, col++)
                    initializeGather(gathers[i++], superGatherRow, col);
                for (int n = 0, row = startRow; n < nSuperGatherRows; n++, row++)
                    if (row != superGatherRow)
                        initializeGather(gathers[i++], row, superGatherCol);
                break;
            case StsPreStackLineSetClass.SUPER_RECT:
                startRow = superGatherRow - nSuperGatherRows / 2;
                startCol = superGatherCol - nSuperGatherCols / 2;
                i = 0;
                for (int r = 0, row = startRow; r < nSuperGatherRows; r++, row++)
                    for (int c = 0, col = startCol; c < nSuperGatherCols; c++, col++)
                        initializeGather(gathers[i++], row, col);
                break;
            case StsPreStackLineSetClass.SUPER_INLINE:
                startCol = superGatherCol - nSuperGatherCols / 2;
                i = 0;
                for (int c = 0, col = startCol; c < nSuperGatherCols; c++, col++)
                    initializeGather(gathers[i++], superGatherRow, col);
                break;
            case StsPreStackLineSetClass.SUPER_XLINE:
                startRow = superGatherRow - nSuperGatherRows / 2;
                i = 0;
                for (int n = 0, row = startRow; n < nSuperGatherRows; n++, row++)
                    initializeGather(gathers[i++], row, superGatherCol);
                break;
        }
        computeDrawOffsets();
        //        checkComputeMuteRange(velocityProfile);
    }

    private void initializeGather(StsGather gather, int row, int col)
    {
        gather.initializeRowCol(row, col);
        gather.setSuperGather(this);
        if (row == superGatherRow && col == superGatherCol)
            centerGather = gather;
        nSuperGatherTraces += gather.nGatherTraces;
        gather.initializeTimeRange();
    }

    private void checkConstructSuperGatherTraceSet()
    {
        if (superGatherTraceSet != null) return;
        superGatherTraceSet = new StsSuperGatherTraceSet(this, nApproxSubIntervals);
    }

    private void setGatherDescription()
    {
        gatherDescription = lineSet.getGatherDescription(superGatherRow, superGatherCol, nSuperGatherTraces);
    }

    public boolean isComputingGather()
    {
        return computingSemblance;
    }

    /** the gather row or col have changed: null out all arrays */
    public void gatherDataChanged()
    {
        superGatherTraceSet = null;
        for (int g = 0; g < nGathersInSuperGather; g++)
            gathers[g].gatherDataChanged();
    }

    /** trace data at this gather has been changed: null out gather data */
    public void gatherTracesChanged()
    {
        superGatherTraceSet = null;
        stackPanelsChanged();
        if (gathers == null) return;
        for (int g = 0; g < nGathersInSuperGather; g++)
            gathers[g].gatherDataChanged();
    }

    public void semblanceChanged()
    {
        semblanceBytes = null;
    }

    public void residualSemblanceChanged()
    {
        residualSemblance = null;
    }

    public void vvsStackChanged()
    {
        vvsTraces = null;
    }

    public void cvsStackChanged()
    {
        cvsTraces = null;
    }

    public void stackPanelsChanged()
    {
        cvsTraces = null;
        vvsTraces = null;
    }

    public boolean velocityProfileChanged(StsVelocityProfile velocityProfile)
    {
        if (velocityProfile == null) return false;
        if (velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
            muteRangeChanged();
        else if (velocityProfile.changeType == StsVelocityProfile.CHANGE_POINT)
        {
            centerGather.velocities.clearVelocities();
            stretchMuteChanged();
        }
        return true;
    }

    private void stretchMuteChanged()
    {
        for (int i = 0; i < gathers.length; i++)
            gathers[i].stretchMuteChanged();
        setStretchMuteFactor();
        setStretchMuteMinOffset();
    }

    private void setStretchMuteFactor()
    {
        double percentStretch = lineSet.getWiggleDisplayProperties().getStretchMute();
        stretchFactor = percentStretch / 100 + 1;
    }

    private void setStretchMuteMinOffset()
    {
        stretchMuteMinOffset = lineSet.getWiggleDisplayProperties().getStretchMuteMinOffset();
    }

    public boolean semblanceChanged(StsVelocityProfile velocityProfile)
    {
        if (velocityProfile == null) return false;
        velocityProfileChanged(velocityProfile);
        if (velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
            semblanceChanged();
        return true;
    }

    public boolean residualSemblanceChanged(StsVelocityProfile velocityProfile)
    {
        if (velocityProfile == null) return false;
        velocityProfileChanged(velocityProfile);
        residualSemblanceChanged();
        return true;
    }

    /*
public void cvsChanged(StsVelocityProfile velocityProfile)
        {
            if(velocityProfile == null) return;
            velocityProfileChanged(velocityProfile);
            if(velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
                cvsChanged();
        }
    */
    public synchronized void displayWiggleTraces(GL gl, StsViewGather viewGather)
    {
        StsPreStackLine line = lineSet.getDataLine(superGatherRow, superGatherCol);
        boolean flatten = lineSet.lineSetClass.getFlatten();
        for (int i = 0; i < gathers.length; i++) gathers[i].flatten = flatten;
        boolean adjustNMO = lineSet.isNMOed || centerGather.flatten;
        for (int i = 0; i < gathers.length; i++)
        {
            gathers[i].checkComputeStretchMuteOffsets();
        }

        if (superGatherTraceSet == null)
            superGatherTraceSet = new StsSuperGatherTraceSet(this, nApproxSubIntervals);
        superGatherTraceSet.rmsNormalizeAmplitudes();

        if (adjustNMO) computeSuperGatherOffsetTimes();
        int offsetAxisType = wiggleProperties.getOffsetAxisType();
        double horizScaleFactor = computeHorizScale(offsetAxisType);
        superGatherTraceSet.initializeDraw(viewGather, horizScaleFactor, wiggleProperties);
        float[][] axisRanges = viewGather.axisRanges;
        StsGatherTrace[] traces = getSuperGatherTraces(offsetAxisType);
        int nTraces = traces.length;
        for (int i = 0; i < nTraces; i++)
        {
            StsGatherTrace trace = traces[i];
            double x0;
            if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
                x0 = trace.x;
            else
                x0 = i;

            if (x0 < axisRanges[0][0] || x0 > axisRanges[0][1]) continue;

            double[] traceMuteRange;
            // float[] traceBeforeNMO = superGatherData.data[index].clone();
            // if (adjustNMO) applyNMO(centerGather.velocityProfile, superGatherOffsets[index], superGatherData.data[index], zInc);
            if (centerGather.flatten)
                traceMuteRange = getFlattenedMuteRange(trace);
            else
                traceMuteRange = getUnFlattenedMuteRange(trace);
            superGatherTraceSet.displayInterpolatedPoints(gl, trace, (float) x0, traceMuteRange);
            // superGatherData.data[index] = traceBeforeNMO;
        }
    }

    public double[] getUnFlattenedMuteRange(StsGatherTrace trace)
    {
        double[] traceMuteRange;
        traceMuteRange = new double[]{semblanceZMin, semblanceZMax};
        double offset = trace.x;
        traceMuteRange[0] = Math.max(semblanceZMin, centerGather.getStretchMuteTZeroFromOffset(offset) * centerGather.stretchFactor);
        return traceMuteRange;
    }

    public double[] getFlattenedMuteRange(StsGatherTrace trace)
    {
        double[] traceMuteRange;
        traceMuteRange = new double[]{semblanceZMin, semblanceZMax};
        double offset = trace.x;
        traceMuteRange[0] = Math.max(semblanceZMin, centerGather.getStretchMuteTZeroFromOffset(offset));
        return traceMuteRange;
    }

    /*
private boolean applyNMO(int index, double[][] drawPoints)
    {
        int traceCount=0;
        for (int i=0; i < gathers.length; i++)
        {
            StsGather gather = gathers[i];
            int remainder = index - traceCount;
            if (remainder < gather.nGatherTraces)
            {
                if (gather.velocities == null || gather.velocities.offsetTimes == null || gather.velocities.offsetTimes.times == null) continue;
                if (gather.velocities.offsetTimes.times.length <= remainder)
                {
                    StsMessage.printMessage("StsSuperGather.applyNMO: remainder >= gather.velocities.offsetTimes.times.length");
                    continue;
                }
                gather.checkAdjustFlattening(drawPoints, remainder);
                return true;
            }
            traceCount += gather.nGatherTraces;
        }
        return false;
    }
    */

    private boolean applyNMO(StsVelocityProfile velProfile, double x, float[] traceData, double zInc)
    {
        for (int i = 0; i < traceData.length; i++)
        {
            double t0 = i * zInc;
            double vrms = velProfile.getVelocityFromTzero((float) t0);
            double tx = Math.sqrt(t0 * t0 + x * x / (vrms * vrms));
            double txData = StsTraceUtilities.getInterpolatedValue(traceData, tx, zInc);
            traceData[i] = (float) txData;
        }
        return true;
    }

    private double computeHorizScale(int offsetAxisType)
    {
        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
        {
            double firstTraceOffset = ((StsGatherTrace) superGatherTraceSet.getFirstTrace()).x;
            double lastTraceOffset = ((StsGatherTrace) superGatherTraceSet.getLastTrace()).x;
            int nTraceSpacings = superGatherTraceSet.nTraces - 1;
            return ((lastTraceOffset - firstTraceOffset) / nTraceSpacings);
        }
        else
            return 1.0;
    }

    /**
     * returns array of all offsets from gathers within supergather
     *
     * @return
     */
    private double[] computeSuperGatherOffsets()
    {
        nSuperGatherTraces = 0;
        for (int i = 0; i < gathers.length; i++) nSuperGatherTraces += gathers[i].nGatherTraces;
        if (nSuperGatherTraces == 0) return null;
        if (gathers == null || gathers.length == 0) return null;
        double[] array = new double[nSuperGatherTraces];
        int trace = 0;
        for (int i = 0; i < gathers.length; i++)
        {
            double[] offsets = gathers[i].getGatherOffsets();
            if (offsets == null) continue;
            for (int j = 0; j < offsets.length; j++)
            {
                if (trace >= array.length)
                {
                    StsMessage.printMessage("StsSuperGather.computeSuperGatherOffsets(): trace > array.length");
                    break;
                }
                array[trace] = offsets[j];
                trace++;
            }
        }
        return array;
    }

    /**
     * gathers all of the traces from the sub-gathers into the supergather
     *
     * @return super gather trace floats
     */
    private float[][] getSuperGatherTraceSet()
    {
        if (gathers == null || gathers.length == 0) return null;

        //...First, get planeData (floats for all traces in supergather
        if (centerGather.line == null) return null;
        float zInc = centerGather.line.getZInc();
        int nSlices = centerGather.line.nSlices;
        nSuperGatherTraces = 0;
        for (int i = 0; i < gathers.length; i++) nSuperGatherTraces += gathers[i].nGatherTraces;
        float[][] gatherData = new float[nSuperGatherTraces][nSlices];
        int nTrace = 0;
        for (int g = 0; g < gathers.length; g++)
        {
            StsGather gather = gathers[g];
            int nGatherTraces = gather.nGatherTraces;
            float[] dataF = gather.getGatherDataf();
            int nGatherSample = 0;
            if (dataF != null)
            {
                for (int t = 0; t < nGatherTraces; t++)
                {
                    System.arraycopy(dataF, nGatherSample, gatherData[nTrace++], 0, nSlices);
                    nGatherSample += nSlices;
                }
            }
            else
            {
                nTrace += nGatherTraces;
            }
        }
        double rmsAmplitude = StsMath.computeNormalizedRMSAmplitude(gatherData);
        for (int n = 0; n < nSuperGatherTraces; n++)
            StsMath.normalizeAmplitude(gatherData[n], rmsAmplitude);
        return gatherData;
    }

    private void computeSuperGatherOffsetTimes()
    {
        for (int i = 0; i < gathers.length; i++)
        {
            if (centerGather.velocityProfile == null) break;
            if (!gathers[i].velocities.computeVelocities(centerGather.velocityProfile)) continue;
            gathers[i].velocities.computeOffsetTimes();
        }
    }

    public String getGatherDescription()
    {
        int nTraces = 0;
        if (centerGather != null) nTraces = centerGather.nGatherTraces;
        int cdp = getCdp();
        if (lineSet instanceof StsPreStackLineSet2d)
        {
            return lineSet.getName() + " CDP " + cdp + " Line " + superGatherRow + " Gather " + superGatherCol + " Traces " + nTraces;
        }
        double inline = getInline();
        double xline = getXline();
        return lineSet.getName() + " Inline " + inline + " Xline " + xline + " CDP " + cdp + " Traces " + nTraces;
    }

    /**
     * gets xline of supergather for 3d
     *
     * @return xline or -1 for null linSet
     */
    public double getXline()
    {
        if (lineSet == null) return -1;
        if (lineSet.isInline)
            return lineSet.getColNumFromCol(superGatherCol);
        return lineSet.getRowNumFromRow(superGatherRow);
    }

    /**
     * gets inline of supergather for 3d
     *
     * @return inline or -1 for null linSet
     */
    public double getInline()
    {
        if (lineSet == null) return -1;
        if (lineSet.isInline)
            return lineSet.getRowNumFromRow(superGatherRow);
        return lineSet.getColNumFromCol(superGatherCol);
    }

    /**
     * CDP from header (expensive, but more trustworthy than calculated CDP
     *
     * @return CDP number or -1 if line of center gather is null
     */
    public int getCdp()
    {
        int cdp = -1;
        if (centerGather != null && centerGather.line != null)
        {
            cdp = centerGather.line.getCDP(centerGather.nLineGather);
        }
        return cdp;
    }

    public void adjustStretchMute(StsPoint pressedPick, StsPoint pick)
    {
        centerGather.adjustStretchMute(pressedPick, pick, velocityProfile, lineSet);
    }

    public void displayNMOCurve(StsGLPanel3d glPanel3d, float[] offsetAxisRange)
    {
        //velocityProfile = centerGather.velocityProfile;
        if (velocityProfile == null) return;
        if (velocityProfile.isInterpolated()) return;

        StsPoint[] points = velocityProfile.getProfilePoints();
        if (points.length == 0) return;

        if (!checkInitializeGather()) return;
        if (!checkComputeVelocitiesAndOffsetTimes()) return;
        if (wiggleProperties == null) return;
        int offsetAxisType = wiggleProperties.getOffsetAxisType();
        glPanel3d.getGL().glLineWidth(1);
        for (int i = 0; i < points.length; i++)
            displayPickNMOPoint(points[i], velocityProfile, StsGather.NONE_INDEX, offsetAxisType, true, false, glPanel3d.getGL(), StsColor.RED);
        glPanel3d.getGL().glLineWidth(2);
    }

    public void pickNMOCurve(StsGLPanel3d glPanel3d, float[] offsetAxisRange)
    {
        // velocityProfile = centerGather.velocityProfile;
        if (velocityProfile == null) return;
        if (velocityProfile.isInterpolated()) return;

        StsPoint[] points = velocityProfile.getProfilePoints();
        if (points.length == 0) return;

        checkComputeVelocitiesAndOffsetTimes();
        centerGather.checkSetOrder();  //getting order of nmo

        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        GL gl = glPanel3d.getGL();
        gl.glDisable(GL.GL_LIGHTING);
        gl.glLineWidth(2.0f);

        boolean pick = true;
        displayPickNMOPoint(velocityProfile.getTopMute(), velocityProfile, TOP_MUTE_INDEX, offsetAxisType, !pick, pick, gl, StsColor.GREEN);
        displayPickNMOPoint(velocityProfile.getBottomMute(), velocityProfile, BOT_MUTE_INDEX, offsetAxisType, !pick, pick, gl, StsColor.BLUE);
        displayPickStretchMute(offsetAxisType, pick, gl, StsColor.CYAN);

        int nPoints = points.length;
        for (int n = 0; n < nPoints; n++)
            displayPickNMOPoint(points[n], velocityProfile, n, offsetAxisType, !pick, pick, gl, StsColor.RED);

        gl.glEnable(GL.GL_LIGHTING);
        gl.glLineWidth(1.0f);

    }

    public void displayPickStretchMute(int offsetAxisType, boolean pick, GL gl, StsColor color)
    {
        try
        {
            if (!centerGather.checkComputeStretchMuteOffsets()) return;
            gl.glDisable(GL.GL_LIGHTING);
            color.setGLColor(gl);
            if (pick)
            {
                gl.glInitNames();
                gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_EDGE);
                gl.glPushName(STRETCH_MUTE_INDEX);
            }

            double unFlatten = 1.0;
            boolean flatten = lineSet.lineSetClass.getFlatten();
            if (!flatten)
            {
                unFlatten = centerGather.stretchFactor;
            }
            if (debug) StsException.systemDebug(this, "displayPickStretchMute", "stretchMute: " + unFlatten);

            double t0 = lineSet.getZMin();
            gl.glBegin(GL.GL_LINE_STRIP);
            StsGatherTrace[] traces = getSuperGatherTraces(offsetAxisType);
            int n = 0;
            for (StsGatherTrace trace : traces)
            {
                double offset = trace.x;
                t0 = centerGather.getStretchMuteTZeroFromOffset(Math.abs(offset));
                if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
                    gl.glVertex2d(offset, t0 * unFlatten);
                else
                    gl.glVertex2d((double) n, t0 * unFlatten);
            }
            if (pick)
            {
                gl.glPopName();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSuperGather.displayStretchMute() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public float getOffsetForTypeValue(float offsetValue)
    {
        return centerGather.getOffsetForTypeValue(offsetValue);
    }


    public float getTraceIndexFloatFromTypeValue(float offsetValue)
    {
        return centerGather.getTraceIndexFloatFromTypeValue(offsetValue);
    }

    /*
public float getTzeroForTypeValue(float offsetValue)
        {
            return centerGather.getTzeroForTypeValue(offsetValue);
        }
    */
    public void setUseSecondOrder(boolean b)
    {
        useSecondOrder = b;
        if (debug) System.out.println("Setting useSecondOrder " + b);
    }

    public boolean checkSetOrder(byte order)
    {
        centerGather.checkSetOrder(order);
        return true;
    }

    /*
public boolean checkComputeVelocities(StsVelocityProfile velocityProfile)
        {
            if(this.velocityProfile != velocityProfile)
            {
                StsException.systemError(this, "checkComputeVelocities", "VelocityProfiles not consistent between superGather and calling class.");
                return false;
            }
            return centerGather.velocities.checkComputeVelocities(velocityProfile);
        }
    */
    /**
     * Check to see if any changes have been made to velocityProfile: if so, recompute.  If velocity range has changed
     * recompute both velocityProfile and inputVelocityProfile if it exists.
     *
     * @return
     */
    public boolean checkComputeVelocities()
    {
        for (int i = 0; i < gathers.length; i++)
        {
            gathers[i].checkComputeVelocities();
        }
        //    velocityProfile = centerGather.velocityProfile;
        return centerGather.checkComputeVelocities();
    }

    public boolean checkInitializeGather()
    {
        for (int i = 0; i < gathers.length; i++)
        {
            gathers[i].checkInitializeGather();
        }
        return centerGather.checkInitializeGather();
    }

    public boolean checkComputeVelocitiesAndOffsetTimes()
    {
        for (int i = 0; i < gathers.length; i++)
        {
            gathers[i].checkComputeVelocitiesAndOffsetTimes();
        }
        return centerGather.checkComputeVelocitiesAndOffsetTimes();
    }

    // TODO check that datumShift is not being applied if no datumShift is available
    public boolean checkSetDatumShift()
    {
        double datumShift = 0.0f; // new datum shift
        double userDatum = 0.0f;
        double userVelocity = 0.0f;

        try
        {
            String datumAttribute = lineSet.datumProperties.getDatumAttribute();
            if (datumAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.NONE]))
                datumShift = 0.0;
            else
            {
                if (!datumAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.USER_SPECIFIED]))
                {
                    double[] datumValues = centerGather.line.getAttributeArray(datumAttribute);
                    if (datumValues != null)
                        userDatum = (float) datumValues[0]; // Told it would be the same on every trace by Bill
                }
                else
                    userDatum = lineSet.datumProperties.getDatum();

                int datumDomain = lineSet.datumProperties.getDatumDomainIndex();
                if ((datumDomain == StsDatumProperties.TIME_DOMAIN && !StsModel.getCurrentModel().getIsDepth()) ||
                    (datumDomain == StsDatumProperties.DEPTH_DOMAIN && StsModel.getCurrentModel().getIsDepth()))
                    datumShift = userDatum;
                else
                {
                    String velocityAttribute = lineSet.datumProperties.getVelocityAttribute();
                    if (!velocityAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.NONE]))
                    {
                        if (!velocityAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.USER_SPECIFIED]))
                        {
                            double[] velocityValues = centerGather.line.getAttributeArray(velocityAttribute);
                            if (velocityValues != null)
                                userVelocity = (float) velocityValues[0]; // Told it would be the same on every trace by Bill
                        }
                        else
                            userVelocity = lineSet.datumProperties.getVelocity();
                        if (userVelocity != 0.0f)
                        {
                            if (StsModel.getCurrentModel().getIsDepth())
                                datumShift = userVelocity * (userDatum / 1000.0f); // userDatum in msec
                            else
                                datumShift = 2.0 * (userDatum / (userVelocity / 1000.0f));  // two-way time msec
                        }
                    }
                    else
                        datumShift = 0.0; // depthDatum, but no velocity so can't shift time data
                }
            }
            if (datumShift == this.datumShift) return false;
            this.datumShift = datumShift;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("Error accessing datum information", e, StsException.WARNING);
            return false;
        }
    }

    public StsGatherTrace[] getCenterGatherTraces()
    {
        return centerGather.getConstructGatherTraces();
    }

    public StsGatherTrace[] getSuperGatherTraces(int offsetAxisType)
    {
        if (superGatherTraceSet == null)
            superGatherTraceSet = new StsSuperGatherTraceSet(this, nApproxSubIntervals);
        return superGatherTraceSet.getGatherTraces(offsetAxisType);
    }

    /*
ComputeResidualSemblanceStandardProcess getResidualSemblanceStandardProcessInstance(StsWin3dBase win3d)
        {
            if(residualSemblanceStandardProcess != null)
                residualSemblanceStandardProcess.initialize(win3d);
            else
                residualSemblanceStandardProcess = new ComputeResidualSemblanceStandardProcess(win3d);
            return residualSemblanceStandardProcess;
        }
    */

    /*
ComputeResidualSemblanceEnhancedProcess getResidualSemblanceEnhancedProcessInstance(StsWin3dBase win3d)
        {
            if(residualSemblanceEnhancedProcess != null)
                residualSemblanceEnhancedProcess.initialize(win3d);
            else
                residualSemblanceEnhancedProcess = new ComputeResidualSemblanceEnhancedProcess(win3d);
            return residualSemblanceEnhancedProcess;
        }
    */
    /*
public float[] computeSemblanceTrace(int nDisplayValues)
        {
            ComputeSemblanceStandardProcess semblanceProcess = new ComputeSemblanceStandardProcess();
            return semblanceProcess.computeSemblanceTrace(nDisplayValues);
        }

        public float[] computeStackedTrace(int nDisplayValues)
        {
            ComputeStackedAmplitudeProcess stackProcess = new ComputeStackedAmplitudeProcess();
            return stackProcess.computeStackedTrace(nDisplayValues);
        }

        public float[] computeStackedTraceNormalized(int nDisplayValues)
        {
            ComputeStackedAmplitudeProcess stackProcess = new ComputeStackedAmplitudeProcess();
            float[] stackedData = stackProcess.computeStackedTrace(nDisplayValues);
            if(stackedData == null) return stackedData;
            StsMath.normalizeAmplitude(stackedData);
            return stackedData;
        }
    */

    public float[] computeSemblanceTrace(int nDisplaySlices)
    {
        try
        {
            //			 profilePoints = lineSet.getVelocityModel().getVelocityProfilePoints(superGatherRow, StsSuperGather.this.superGatherCol);
            //			 if (profilePoints == null) return null;
            //			 profilePoints = (StsPoint[])StsToolkit.deepCopy(profilePoints);


            checkInitializeSuperGatherGeometry();

            Semblance semblance = new Semblance(nDisplaySlices);

            for (int g = 0; g < nGathersInSuperGather; g++)
            {
                StsGather gather = gathers[g];
                if (!gather.checkInitializeGather()) continue;
                if (lineSet.isNMOed && !checkComputeVelocitiesAndOffsetTimes()) continue;
                gather.computeSemblanceTrace(semblance);
            }
            return semblance.getSemblanceValues();
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeSemblanceStandard", e);
            return null;
        }
    }


    abstract private class ComputeSemblanceProcess implements Runnable
    {
        StsProgressBarDialog progressDialog;
        //double[][] sumValues;
        //double[][] sumValuesSq;
        //int[][] nSumValues;
        Semblance semblance;
        int nVels;
        float vMin;
        float vInc;
        int nDataSamples;
        float zMin;
        float zMax;
        float zInc;
        int nFirstGatherLineTrace;
        int nGatherTraces;
        double[] traceOffsets;
        float maxSemblance;
        float minSemblance;
        double[] semblances;
        double[] phase;
        byte[] localSemblanceBytes;
        boolean autoscale = false;
        float maxScale = 1.0f;
        StsGather velocityGather;
        boolean useWindow = true;
        double semblanceTotalSum = 0.0;
        int nBytesPerSample = 1;
        // double[][] gatherData;

        abstract protected boolean initializeSemblance();
        abstract protected void computeSemblanceTrace(double vf, int nInc, int offset, int nTraceValues);
        abstract protected void semblanceSuccessful();

        public ComputeSemblanceProcess()
        {
        }

        protected void initialize()
        {
            if (debug) System.out.println(StsToolkit.getSimpleClassname(this) + " initialized.");
            computingSemblance = true;
            constructDialog();
            stretchMuteChanged(); //force stretch mute to update
        }

        protected void constructDialog()
        {
            Runnable runnable = new Runnable()
            {
                public void run()
                {
                    progressDialog = StsProgressBarDialog.constructor(parentWindow, "Calculating...");
                }
            };
            StsToolkit.runWaitOnEventThread(runnable);
        }

        public boolean isCanceled()
        {
            return progressDialog.isCanceled();
        }

        public void run()
        {
            computeSemblance();
        }

        public boolean computeSemblance()
        {
            try
            {
                if (debug) debugThread("computing semblance standard");
                if (!filesOk) return false;
                checkInitializeSuperGatherGeometry();
                if(lineSet.isNMOed) checkComputeVelocitiesAndOffsetTimes();
                initialize();
                initializeSemblance();
                StsSemblanceComputeProperties computeProperties = lineSet.semblanceComputeProperties;
                int windowWidth = computeProperties.windowWidth;

                if (progressDialog != null) progressDialog.setProgressMax(nVels);

                if (superGatherTraceSet.nTraces == 0)
                {
                    semblanceBytes = new byte[]{nullByte};
                    return false;
                }
                initializeSemblanceOutput();
                float v = vMin;
                int offset = 0;
                for (int n = 0; n < nVels; n++, v += vInc, offset += nSemblanceSamples)
                {
                    if (progressDialog != null && progressDialog.isCanceled()) return false;
                    computeSemblanceTrace(v, n, offset, nSemblanceSamples);
                    if (useWindow && windowWidth > 0)
                        StsMath.windowAverages(semblances, windowWidth, offset, nSemblanceSamples);
                    if (progressDialog != null) progressDialog.setProgress(n);
                }
                // maxSemblance = StsMath.max(semblances);
                semblanceSuccessful();
                return true;
            }
            catch (Exception e)
            {
                semblanceBytes = new byte[]{nullByte};
                errorMessage = e.getMessage();
                StsException.outputWarningException(this, "computeSemblanceStandardProcess", e);
                return false;
            }
            finally
            {
                completeProcess();
            }
        }

        protected void initializeSemblanceOutput()
        {
            semblances = new double[nVels * nSemblanceSamples];
            localSemblanceBytes = new byte[nVels * nSemblanceSamples * nBytesPerSample];
        }

        /*
        protected boolean initializeSemblanceGather(StsGather gather)
        {
            this.gather = gather;
            nDataSamples = gather.line.nSlices;
            zMin = gather.line.zMin;
            zMax = gather.line.zMin;
            zInc = gather.line.zInc;
            nFirstGatherLineTrace = gather.nFirstGatherLineTrace;
            nGatherTraces = gather.nGatherTraces;
            traceOffsets = gather.line.traceOffsets;

            float[] data = gather.getGatherDataf();
            if(data == null) return false;

            gatherData = new double[nGatherTraces][nDataSamples];
            int ii = 0;
            for(int t = 0; t < nGatherTraces; t++)
                for(int s = 0; s < nDataSamples; s++)
                    gatherData[t][s] = data[ii++];
            return true;
        }
        */
        protected void completeProcess()
        {
            computingSemblance = false;
            blockReprocessSemblance = false;
            if (progressDialog != null)
            {
                progressDialog.finished();
                progressDialog.dispose();
                progressDialog = null;
            }
            if (computeSemblanceFailed())
                new StsMessage(model.win3d, StsMessage.WARNING, "Gather compute process has failed with error: " + errorMessage);
            if (debug) System.out.println(StsToolkit.getSimpleClassname(this) + " completed.");
            model.viewObjectRepaintFamily(parentWindow, lineSet);
        }

        protected boolean computeSemblanceFailed()
        {
            return semblanceBytes != null && semblanceBytes.length == 1;
        }
    }

    public boolean computeSemblanceStandardProcess(StsPreStackLineSet lineSet, StsProgressBarDialog dialog)
    {
        ComputeSemblanceStandardProcess computeSemblanceProcess = new ComputeSemblanceStandardProcess();
        return computeSemblanceProcess.computeSemblance();
    }


    class ComputeSemblanceStandardProcess extends ComputeSemblanceProcess
    {
        public ComputeSemblanceStandardProcess()
        {
        }

        protected void computeSemblanceTrace(double v, int nVel, int offset, int nTraceValues)
        {
            double rvsq = 1.0 / (v * v);
            double t0 = semblanceZMin;
            double tInc = semblanceZInc;
            StsGatherTrace[] gatherTraces = superGatherTraceSet.getGatherTraces();
            for (int s = 0; s < nSemblanceSamples; s++, t0 += tInc, offset++)
            {
                double sumValue = 0.0;
                double sumValueSq = 0.0;
                int nVals = 0;
                for (StsGatherTrace trace : gatherTraces)
                {
                    double xsq = trace.xsq;
                    double tOffset = Math.sqrt(t0 * t0 + xsq * rvsq);
                    // if(!centerGather.isInsideMuteRange(trace, tOffset, t0)) continue;
                    double tValue;
                    if (!lineSet.isNMOed)
                        tValue = tOffset;
                    else
                        tValue = velocityGather.inputVelocities.computeTzero(trace, tOffset);
                    double value = superGatherTraceSet.getTOffsetValue(tValue, trace);
                    nVals++;
                    if (value == StsParameters.doubleNullValue) continue;
                    sumValue += value;
                    sumValueSq += value * value;
                    if (debugSemblance)
                    {
                        if (t0 == 1500 && (v > 2.03 && v < 2.05))
                        {
                            StsException.systemDebug(this, "computeSemblanceStandardProcess", "t0: " + t0 + " vel: " + v + " xsq: " + xsq + " t: " + tOffset +
                                "\n    val: " + value + " sumValue: " + sumValue + " sumValueSq: " + sumValueSq + " nVel:" + nVel + " s: " + s);
                        }
                    }
                }

                double semblance = (nVals < MIN_COHERENCE_VALS || sumValueSq == 0.0) ? 0.0 : sumValue * sumValue / sumValueSq / nVals;
                semblances[offset] = semblance;
                semblanceTotalSum += semblance;
                maxSemblance = Math.max(maxSemblance, (float)semblance);
            }
        }

        protected boolean initializeSemblance()
        {
            checkConstructSuperGatherTraceSet();
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            nBytesPerSample = 1;
            StsJOGLShader.nBytesPerTexel = 1;
            vMin = rangeProperties.velocityMin;
            float vMax = rangeProperties.velocityMax;
            vInc = rangeProperties.velocityStep;
            nVels = 1 + StsMath.ceiling((vMax - vMin) / vInc);
            nSemblanceSamples = (int) Math.round((semblanceZMax - semblanceZMin) / semblanceZInc) + 1;
            semblance = new Semblance(nSemblanceSamples);
            //sumValues = new double[nVels][nSemblanceSamples];
            //sumValuesSq = new double[nVels][nSemblanceSamples];
            //nSumValues = new int[nVels][nSemblanceSamples];
            StsSemblanceComputeProperties computeProperties = lineSet.semblanceComputeProperties;
            autoscale = computeProperties.normalizeCoherence;
            maxScale = (float)computeProperties.clip;
            maxSemblance = 0.0f;
            velocityGather = centerGather;
            return true;
        }

        protected void semblanceSuccessful()
        {
            float[] scaleData = computeSemblanceScaleAndOffset();
            float scale= scaleData[0];
            float scaleOffset = scaleData[1];

            for (int i = 0; i < semblances.length; i++)
                localSemblanceBytes[i] = StsMath.floatToUnsignedByte254WithScale((float) semblances[i], scale, scaleOffset);

            semblanceBytes = localSemblanceBytes;
            lineSet.semblanceColorList.colorscale.setRange(0, maxSemblance);
        }

        protected float[] computeSemblanceScaleAndOffset()
        {
            float scale, scaleOffset;
            if (autoscale)
            {
                minSemblance = (float)(semblanceTotalSum/(nVels*nSemblanceSamples));
                if(minSemblance < 0.001f) minSemblance = 0.001f;
                scale = 254*colorSaturationFactor*5/(maxSemblance - minSemblance);
                scaleOffset = -(float)minSemblance*scale;
            }
            else
            {
                scale = 254*colorSaturationFactor*maxScale/maxSemblance;
                scaleOffset = 0.0f;
                maxSemblance = maxScale;
            }
            return new float[] { scale, scaleOffset };
        }

        /*
        public float[] computeSemblanceTrace()
        {
            try
            {
                initializeSemblance();

                checkInitializeSuperGatherGeometry();
                if (lineSet.isNMOed && !checkComputeVelocitiesAndOffsetTimes()) return null;
                Semblance semblance = new Semblance(nDataSamples);
                for (int g = 0; g < nGathersInSuperGather; g++)
                {
                    StsGather gather = gathers[g];
                    if (!gather.checkInitializeGather()) continue;
                    initializeSemblanceGather(gather);
                    computeSemblanceTrace(semblance);
                }
                return semblance.getSemblanceValues();
            }
            catch (Exception e)
            {
                residualSemblance = null;
                StsException.outputWarningException(this, "computeSemblanceStandard", e);
                return null;
            }
        }
        */
        /*
        public float[] computeSemblanceTrace(int nDisplayValues)
        {
            try
            {
                initializeSemblance();

                checkInitializeSuperGatherGeometry();
                if (lineSet.isNMOed && !checkComputeVelocitiesAndOffsetTimes()) return null;
                Semblance semblance = new Semblance(nDisplayValues);
                for (int g = 0; g < nGathersInSuperGather; g++)
                {
                    StsGather gather = gathers[g];
                    if (!gather.checkInitializeGather()) continue;
                    initializeSemblanceGather(gather);
                    computeSemblanceTrace(semblance);
                }
                return semblance.getSemblanceValues();
            }
            catch (Exception e)
            {
                residualSemblance = null;
                StsException.outputWarningException(this, "computeSemblanceTrace", e);
                return null;
            }
        }
        */
    }


    public boolean computeSemblanceEnhancedProcess(StsPreStackLineSet lineSet, StsProgressBarDialog dialog)
    {
        ComputeSemblanceEnhancedProcess computeSemblanceProcess = new ComputeSemblanceEnhancedProcess();
        return computeSemblanceProcess.computeSemblance();
    }

    static public final boolean useHilbertTransform = false;
    static public final boolean phaseDisplay = true;
    static public final boolean phaseSuppressAmplitude = false;

    class ComputeSemblanceEnhancedProcess extends ComputeSemblanceStandardProcess
    {
        boolean enhancedAdjust;
        double phaseAdjust = 30.0;

        public ComputeSemblanceEnhancedProcess()
        {
        }

        protected void computeSemblanceTrace(double v, int nVel, int offset, int nTraceValues)
        {
            try
            {
                double rvsq = 1.0 / (v * v);
                double t0 = semblanceZMin;
                double tInc = semblanceZInc;
                StsGatherTrace[] gatherTraces = superGatherTraceSet.getGatherTraces();
                int nTraces = gatherTraces.length;
                for (int s = 0; s < nSemblanceSamples; s++, t0 += tInc, offset++)
                {
                    int nVals = 0;
                    double phaseRAvg = 0.0;
                    double phaseIAvg = 0.0;
                    for (int n = 0; n < nTraces; n++)
                    {
                        StsGatherTrace trace = gatherTraces[n];
                        double xsq = trace.xsq;
                        double tOffset = Math.sqrt(t0 * t0 + xsq * rvsq);
                        // if(!centerGather.isInsideMuteRange(trace, tOffset, t0)) continue;
                        double tValue;
                        if (!lineSet.isNMOed)
                            tValue = tOffset;
                        else
                            tValue = velocityGather.inputVelocities.computeTzero(trace, tOffset);
                        double[] complex = superGatherTraceSet.getComplexComponents(tValue, trace); // returns x and y complex components
                        if(complex == null) continue;
                        if(phaseSuppressAmplitude) StsMath.normalize(complex);
                        nVals++;
                        phaseRAvg += complex[0];
                        phaseIAvg += complex[1];
                    }
                    if(nVals == 0) continue;
                    phaseRAvg /= nVals;
                    phaseIAvg /= nVals;
                    double semblance = StsMath.length(phaseRAvg, phaseIAvg);
                    semblances[offset] = semblance;
                    semblanceTotalSum += semblance;
                    maxSemblance = Math.max(maxSemblance, (float)semblance);
                    if(phaseDisplay) phase[offset] = phaseRAvg/semblance;
                }
            }
            catch (Exception e)
            {
                StsException.outputWarningException(this, "computeSemblanceTrace", e);
            }
        }

        protected boolean initializeSemblance()
        {
            checkConstructSuperGatherTraceSet();
            if(!useHilbertTransform)
                superGatherTraceSet.computeInstantaneousAmpAndPhaseGatherData();
            useWindow = false;
            StsSemblanceComputeProperties computeProperties = lineSet.semblanceComputeProperties;
            enhancedAdjust = computeProperties.enhancedAdjust;
            boolean ok = super.initializeSemblance();
            minSemblance = 0.0f;
            if(phaseDisplay)
            {
                nBytesPerSample = 2;
                StsJOGLShader.nBytesPerTexel = 2;
            }
            return ok;
        }

        protected void initializeSemblanceOutput()
        {
            super.initializeSemblanceOutput();
            phase = new double[nVels * nSemblanceSamples];
            if(useHilbertTransform)
                computeHilbertTransform(3);
        }

        private void computeHilbertTransform(int windowHalfSize)
        {
            StsGatherTrace[] gatherTraces = superGatherTraceSet.getGatherTraces();
            int nTraces = gatherTraces.length;
            for (int n = 0; n < nTraces; n++)
            {
                StsGatherTrace trace = gatherTraces[n];
                trace.computeHilbertTransform(windowHalfSize);
            }
        }

        /** Enhanced semblance is displayed as phase with the amplitude controlling the transparency */
        protected void semblanceSuccessful()
        {
            if(!phaseDisplay)
            {
                super.semblanceSuccessful();
                return;
            }
            float[] scaleData = computeSemblanceScaleAndOffset();
            float semblanceScale= scaleData[0];
            float semblanceScaleOffset = scaleData[1];
            int nValues = nVels * nSemblanceSamples;
            int i = 0;
            for (int n = 0; n < nValues; n++)
            {
                localSemblanceBytes[i++] = StsMath.floatMinusOneToOneToUnsignedByte254((float) phase[n]);
                if(semblances[n] < minSemblance)
                    localSemblanceBytes[i++] = (byte)0;
                else
                    localSemblanceBytes[i++] = StsMath.floatToUnsignedByte254WithScale((float)semblances[n], semblanceScale, semblanceScaleOffset);
                // localSemblanceBytes[i++] = (byte)0;
            }
            semblanceBytes = localSemblanceBytes;
            lineSet.semblanceColorList.colorscale.setRange(-1, 1);
        }

        /** suppress semblance below average */
        protected float[] computeSemblanceScaleAndOffset()
        {
            float scale, scaleOffset;
            /* suppress semblance below average */
            minSemblance = (float)(semblanceTotalSum/(nVels*nSemblanceSamples));
            scale = 254/(maxSemblance - minSemblance);
            if (autoscale) scale *= colorSaturationFactor;
            scaleOffset = -minSemblance*scale;
            return new float[] { scale, scaleOffset };
        }
    }

    /** Semblance inner class is used here, but only to compute summed amplitudes and not semblance. */
    private class ComputeStackedAmplitudeProcess extends ComputeSemblanceStandardProcess
    {
        public ComputeStackedAmplitudeProcess()
        {
        }

        protected void computeSemblanceTrace(double v, int nVel, int offset, int nTraceValues)
         {
             double rvsq = 1.0 / (v * v);
             double t0 = semblanceZMin;
             double tInc = semblanceZInc;
             StsGatherTrace[] gatherTraces = superGatherTraceSet.getGatherTraces();
             for (int s = 0; s < nSemblanceSamples; s++, t0 += tInc, offset++)
             {
                 double sumValue = 0.0;
                 int nVals = 0;
                 for (StsGatherTrace trace : gatherTraces)
                 {
                     double xsq = trace.xsq;
                     double tOffset = Math.sqrt(t0 * t0 + xsq * rvsq);
                     // if(!centerGather.isInsideMuteRange(trace, tOffset, t0)) continue;
                     double tValue;
                     if (!lineSet.isNMOed)
                         tValue = tOffset;
                     else
                         tValue = velocityGather.inputVelocities.computeTzero(trace, tOffset);
                     double value = superGatherTraceSet.getTOffsetValue(tValue, trace);
                     nVals++;
                     if (value == StsParameters.doubleNullValue) continue;
                     sumValue += value;
                     if (debugSemblance)
                     {
                         if (t0 == 1500 && (v > 2.03 && v < 2.05))
                         {
                             StsException.systemDebug(this, "ComputeStackedAmplitudeProcess", "t0: " + t0 + " vel: " + v + " xsq: " + xsq + " t: " + tOffset +
                                 "\n    val: " + value + " sumValue: " + sumValue + " nVel:" + nVel + " s: " + s);
                         }
                     }
                 }
                 double semblance = (nVals < MIN_COHERENCE_VALS) ? 0.0 : Math.abs(sumValue) / nVals;
                 semblances[offset] = semblance;
                 semblanceTotalSum += semblance;
                 maxSemblance = Math.max(maxSemblance, (float)semblance);
             }
         }
    }

    private class ComputeResidualSemblanceStandardProcess extends ComputeSemblanceStandardProcess
    {

        ComputeResidualSemblanceStandardProcess()
        {
        }

        protected void computeSemblanceTrace(double vFactor, int nVel)
        {
            velocityGather.velocities.nmoCoefficients.setFactor(vFactor);
            double t0 = zMin;
            //double[] stretchMuteOffsets = gather.compute2ndOrderStretchMuteOffsets(velocityGather.velocityProfile, vFactor);
            StsGatherTrace[] gatherTraces = superGatherTraceSet.getGatherTraces();
            for (int m = 0; m < nSemblanceSamples; m++, t0 += zInc)
            {
                double sumValue = 0.0;
                double sumValueSq = 0.0;
                int nVals = 0;
                int nTrace = nFirstGatherLineTrace;
                velocityGather.velocities.setVelocityFromTZero(t0, vFactor);
                //System.out.println("BB: rvsq=" + velocities.nmoCoefficients.c2);
                for (StsGatherTrace trace : gatherTraces)
                {
                    double x = traceOffsets[nTrace];
                    float[] data = trace.data;
                    double tOffset = velocityGather.computeOffsetTime(x);
                    if (!centerGather.isInsideMuteRange(trace, tOffset, t0)) continue;
                    double tf;
                    if (!lineSet.isNMOed)
                        tf = (tOffset - zMin) / zInc;
                    else
                    {
                        double tZeroInput = velocityGather.inputVelocities.computeTzero(trace, tOffset);
                        tf = (tZeroInput - zMin) / zInc;
                    }
                    if (tf < 0.0 || tf > nDataSamples - 1)
                    {
                        nVals++;
                        continue;
                    }
                    double value;
                    int index = StsMath.floor(tf);
                    if (index == nDataSamples - 1)
                    {
                        value = data[index];
                    }
                    else
                    {
                        tf = tf - index;
                        double v0 = data[index];
                        double v1 = data[index + 1];
                        value = v0 + tf * (v1 - v0);
                    }
                    sumValue += value;
                    sumValueSq += value * value;
                    nVals++;
                }

                semblance.sumValue[m] += sumValue;
                semblance.sumValueSq[m] += sumValueSq;
                semblance.nVals[m] += nVals;
            }
        }

        protected boolean initializeSemblance()
        {
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            float percentRange = rangeProperties.percentRange;
            float percentRangeInc = rangeProperties.percentRangeInc;
            int nHalfIncs = Math.round(percentRange / percentRangeInc);
            nVels = 1 + 2 * nHalfIncs;
            vInc = percentRangeInc / 100;
            vMin = 1 - nHalfIncs * vInc;

            semblance = new Semblance(nSemblanceSamples);
            //sumValues = new double[nVels][nSemblanceSamples];
            //sumValuesSq = new double[nVels][nSemblanceSamples];
            //nSumValues = new int[nVels][nSemblanceSamples];

            StsSemblanceComputeProperties computeProperties = lineSet.semblanceComputeProperties;
            autoscale = computeProperties.normalizeCoherence;
            maxScale = (float)computeProperties.clip;

            velocityGather = centerGather;
            return true;
        }

        protected void semblanceSuccessful()
        {
            if (autoscale)
                maxSemblance = maxSemblance * maxScale;
            else
                maxSemblance = maxScale;

            for (int i = 0; i < semblances.length; i++)
                localSemblanceBytes[i] = StsMath.floatToUnsignedByte254((float) semblances[i], 0, maxSemblance);

            residualSemblance = localSemblanceBytes;
            lineSet.semblanceColorList.colorscale.setRange(0, maxSemblance);
        }
    }

    public synchronized boolean computeCVSPanelsProcess(StsWin3dBase win3d)
    {
        if (blockReprocessSemblance)
        {
            if (debug) debugThread("computeCVSPanelsProcess running");
            return false;
        }
        blockReprocessSemblance = true;
        ComputeCVSPanelsProcess computeProcess = new ComputeCVSPanelsProcess();
        cvsThread = StsToolkit.checkRunRunnable(computeProcess);
        return true;
    }


    public synchronized boolean computeVVSPanelsProcess()
    {
        if (blockReprocessSemblance)
        {
            if (debug) debugThread("computeVVSPanelsProcess running");
            return false;
        }
        blockReprocessSemblance = true;
        ComputeVVSPanelsProcess computeProcess = new ComputeVVSPanelsProcess();
        vvsThread = StsToolkit.checkRunRunnable(computeProcess);
        return true;
    }

    abstract private class ComputeStackPanelsProcess implements Runnable
    {
        StsProgressBarDialog progressDialog;
        public double[][][] superGatherData;
        float[][][] localStackedTraces;
        int nPanels;
        double vMin;
        double vInc;
        int nTracesPerPanel;
        boolean ignoreSuperGather;
        boolean isInline = true;
        float[] zeroTrace;
        int row, col;
        int rowInc, colInc;
        int orientation;
        StsSuperGather superGatherComputer;
        StsGather velocityGather;
        double[] traceOffsets;
        int nFirstGatherLineTrace;
        int nGatherTraces;
        int nTraceSamples;
        float zMin;
        float zInc;

        abstract protected void initializePanel();

        abstract protected float[] computeStackTrace(StsGather gather, double v, int nGather);

        abstract protected void completeProcess(boolean successful);

        public ComputeStackPanelsProcess()
        {
        }

        protected void initializeProcess()
        {
            if (debug) System.out.println(StsToolkit.getSimpleClassname(this) + " initialized.");
            computingStacks = true;
            constructDialog();
        }

        protected void constructDialog()
        {
            Runnable runnable = new Runnable()
            {
                public void run()
                {
                    progressDialog = StsProgressBarDialog.constructor(parentWindow, "Calculating...");
                }
            };
            StsToolkit.runWaitOnEventThread(runnable);
        }

        public void run()
        {
            byte order = lineSet.semblanceComputeProperties.order;
            if (order == StsSemblanceComputeProperties.ORDER_FAST)
                computePanelFast();
            else
                computePanel();
        }

        public void computePanel()
        {
            try
            {
                if (!filesOk) return;
                if (nSuperGatherTraces == 0) return;
                initializeProcess();
                initializePanel();
                initializeOrientation();
                if (StsSuperGather.debugStackTimer)
                {
                    StsSuperGather.StackTotalTimer = StsTimer.constructor("Total Stacking Time");
                    StsSuperGather.StackingTimer = StsTimer.constructor("Inside Stack Loop");
                    StsSuperGather.LoadGathersTimer = StsTimer.constructor("Loading Data From Gathers");
                    StsSuperGather.ApplyNMOTimer = StsTimer.constructor("Applying NMO");
                    StsSuperGather.CalcVelocitiesTimer = StsTimer.constructor("Calculating Velocities");
                    StsSuperGather.CalcMeanTimer = StsTimer.constructor("Calculating Mean of NMO Curve");
                    StsSuperGather.CalcTxTimer = StsTimer.constructor("Calculate Tx");
                    StsSuperGather.InterpolateTracesTimer = StsTimer.constructor("Interpolate To Get NMO Curve");

                    StsSuperGather.StackTotalTimer.start();
                }
                for (int t = 0; t < nTracesPerPanel; t++, row += rowInc, col += colInc)
                {
                    if (!superGatherComputer.initializeGeometry(row, col))
                        zeroTraceOnPanels(localStackedTraces, t, nPanels, zeroTrace);
                    else
                    {
                        if (ignoreSuperGather)
                        {
                            superGatherComputer.nGathersInSuperGather = 1;
                            superGatherComputer.gatherType = StsPreStackLineSetClass.SUPER_SINGLE;
                            superGatherComputer.initializeGathers();
                        }
                        double velocity = vMin;
                        if (StsSuperGather.debugTimer) StsSuperGather.LoadGathersTimer.start();
                        initializeSuperGatherDataArrays(superGatherComputer.gathers);
                        if (StsSuperGather.debugTimer) StsSuperGather.LoadGathersTimer.stopAccumulate();
                        for (int p = 0; p < nPanels; p++, velocity += vInc)
                            localStackedTraces[p][t] = computeStackTrace(velocity, superGatherComputer.gathers);
                    }
                    progressDialog.setLabelText("Computing superGather " + (t + 1));
                    progressDialog.setProgress(t + 1);
                }
                if (StsSuperGather.debugStackTimer)
                {
                    StsSuperGather.StackTotalTimer.stopAccumulate();

                    StsSuperGather.StackTotalTimer.printElapsedTime();
                    StsSuperGather.StackingTimer.printElapsedTime();
                    StsSuperGather.LoadGathersTimer.printElapsedTime();
                    StsSuperGather.ApplyNMOTimer.printElapsedTime();
                    StsSuperGather.CalcVelocitiesTimer.printElapsedTime();
                    StsSuperGather.CalcMeanTimer.printElapsedTime();
                    StsSuperGather.CalcTxTimer.printElapsedTime();
                    StsSuperGather.InterpolateTracesTimer.printElapsedTime();

                    StsSuperGather.StackTotalTimer.reset();
                    StsSuperGather.StackingTimer.reset();
                    StsSuperGather.LoadGathersTimer.reset();
                    StsSuperGather.ApplyNMOTimer.reset();
                    StsSuperGather.CalcVelocitiesTimer.reset();
                    StsSuperGather.CalcMeanTimer.reset();
                    StsSuperGather.CalcTxTimer.reset();
                    StsSuperGather.InterpolateTracesTimer.reset();
                }
                completeProcess(true);
            }
            catch (Exception e)
            {
                StsException.outputWarningException(this, "run", e);
                completeProcess(false);
            }
        }

        public void computePanelFast()
        {
            try
            {
                if (!filesOk) return;
                if (nSuperGatherTraces == 0) return;
                initializeProcess();
                initializePanel();          //this is where localStackedTraces is allocated - nTracesPerPanel is left alone for this
                if (!ignoreSuperGather)
                    nTracesPerPanel += Math.max(nSuperGatherCols, nSuperGatherRows) - 1; //we will need additional traces for the tracemix if using supergather
                initializeOrientation();
                if (StsSuperGather.debugStackTimer)
                {
                    StsSuperGather.StackTotalTimer = StsTimer.constructor("Total Stacking Time");
                    StsSuperGather.StackingTimer = StsTimer.constructor("Inside Stack Loop");
                    StsSuperGather.LoadGathersTimer = StsTimer.constructor("Loading Data From Gathers");
                    StsSuperGather.ApplyNMOTimer = StsTimer.constructor("Applying NMO");
                    StsSuperGather.CalcVelocitiesTimer = StsTimer.constructor("Calculating Velocities");
                    StsSuperGather.CalcMeanTimer = StsTimer.constructor("Calculating Mean of NMO Curve");
                    StsSuperGather.CalcTxTimer = StsTimer.constructor("Calculate Tx");
                    StsSuperGather.InterpolateTracesTimer = StsTimer.constructor("Interpolate To Get NMO Curve");

                    StsSuperGather.StackTotalTimer.start();
                }
                if (StsSuperGather.debugTimer) StsSuperGather.LoadGathersTimer.start();
                float[][] gatherTraces = new float[nTracesPerPanel][];
                double[][] gatherOffsets = new double[nTracesPerPanel][];
                for (int t = 0; t < nTracesPerPanel; t++, row += rowInc, col += colInc)      //nTracesPerPanel includes supergather traces at this point
                {
                    if (!superGatherComputer.initializeGeometry(row, col)) continue;
                    if (superGatherComputer.centerGather != null && superGatherComputer.centerGather.initializeGather())
                    {
                        gatherTraces[t] = superGatherComputer.centerGather.getGatherDataf().clone();
                        gatherOffsets[t] = superGatherComputer.centerGather.getGatherOffsets();
                    }
                }
                if (StsSuperGather.debugTimer) StsSuperGather.LoadGathersTimer.stopAccumulate();
                double vFactor = vMin;
                StsPoint[] velPoints1 = new StsPoint[velocityProfile.getProfilePoints().length];
                for (int i = 0; i < velPoints1.length; i++) velPoints1[i] = new StsPoint(0, 0);
                float[][] panelTraces = new float[nTracesPerPanel][];
                for (int p = 0; p < nPanels; p++, vFactor += vInc)
                {
                    if (StsSuperGather.debugTimer) StsSuperGather.CalcVelocitiesTimer.start();
                    calcVelPoints(velPoints1, vFactor);
                    StsPoint[] velPoints2 = expandVels(velPoints1);
                    if (StsSuperGather.debugTimer) StsSuperGather.CalcVelocitiesTimer.stopAccumulate();
                    for (int t = 0; t < nTracesPerPanel; t++)
                    {
                        progressDialog.setLabelText("Computing panel " + (p + 1));
                        progressDialog.setProgress(p + 1);
                        panelTraces[t] = computeStackTraceFast(gatherTraces[t], gatherOffsets[t], velPoints2);
                    }
                    for (int i = 0; i < localStackedTraces[p].length; i++)
                    {
                        int nSlices = centerGather.line.nSlices;
                        float[] trace = new float[nSlices];
                        if (!ignoreSuperGather)
                        {
                            int nStack = nTracesPerPanel - localStackedTraces[p].length + 1;
                            StsTraceUtilities.stack(panelTraces, i, i + nStack, trace); //trace mix of centerGather stacks is same as stacking superGather
                        }
                        else
                            trace = panelTraces[i];
                        if (trace != null)
                            localStackedTraces[p][i] = StsSeismicFilter.filter(StsFilterProperties.POSTSTACK, trace, 1, nSlices, zInc, lineSet.filterProperties, lineSet.agcProperties, lineSet.dataMin, lineSet.dataMax);
                        else
                        {
                            zeroTraceOnPanels(localStackedTraces, i, nPanels, zeroTrace);
                        }
                    }

                }
                if (StsSuperGather.debugStackTimer)
                {
                    StsSuperGather.StackTotalTimer.stopAccumulate();

                    StsSuperGather.StackTotalTimer.printElapsedTime();
                    StsSuperGather.StackingTimer.printElapsedTime();
                    StsSuperGather.LoadGathersTimer.printElapsedTime();
                    StsSuperGather.ApplyNMOTimer.printElapsedTime();
                    StsSuperGather.CalcVelocitiesTimer.printElapsedTime();
                    StsSuperGather.CalcMeanTimer.printElapsedTime();
                    StsSuperGather.CalcTxTimer.printElapsedTime();
                    StsSuperGather.InterpolateTracesTimer.printElapsedTime();

                    StsSuperGather.StackTotalTimer.reset();
                    StsSuperGather.StackingTimer.reset();
                    StsSuperGather.LoadGathersTimer.reset();
                    StsSuperGather.ApplyNMOTimer.reset();
                    StsSuperGather.CalcVelocitiesTimer.reset();
                    StsSuperGather.CalcMeanTimer.reset();
                    StsSuperGather.CalcTxTimer.reset();
                    StsSuperGather.InterpolateTracesTimer.reset();
                }
                completeProcess(true);
            }
            catch (Exception e)
            {
                StsException.outputWarningException(this, "run", e);
                completeProcess(false);
            }
        }

        public abstract void calcVelPoints(StsPoint[] velPoints1, double factor);

        /**
         * non-linearly add more vels to velocity profile
         * (exponentially more velocities at earlier times)
         * <p/>
         * new points are interpolated from user picks (same as first pick before first pick, extrapolated after last pick)
         *
         * @param velPoints
         * @return
         */
        private StsPoint[] expandVels(StsPoint[] velPoints)
        {
            if (velPoints == null || velPoints.length == 0) return null;
            ArrayList<StsPoint> newVels = new ArrayList<StsPoint>();
            double zMin = lineSet.zMin;
            double zMax = lineSet.zMax;
            double zInc = lineSet.zInc;
            double t1 = 0, t2 = 0, t = zMin;
            double v1 = 0, v2 = 0, v = 0;
            double delV = 0;
            double tInc = zInc;
            int inc = 2;

            t1 = t2 = velPoints[0].v[1];
            v1 = v2 = v = velPoints[0].v[0];

            //before picks
            while (t < t1)
            {
                newVels.add(new StsPoint((float) v, (float) t));
                t += tInc;
                tInc += inc;
            }

            //in the middle of the picks
            for (int i = 1; i < velPoints.length; i++)
            {
                t = t1;  //we want to be sure to include elbows of actual picks too!!!
                t2 = velPoints[i].v[1];
                v2 = velPoints[i].v[0];
                delV = (v2 - v1) / (t2 - t1);
                while (t < t2)
                {
                    v = v1 + (t - t1) * delV;
                    newVels.add(new StsPoint((float) v, (float) t));
                    t += tInc;
                    tInc += inc;
                }
                t1 = t2;
                v1 = v2;
            }

            //after picks
            t = t2;
            while (t < zMax)
            {
                v = v2 + (t - t2) * delV;
                newVels.add(new StsPoint((float) v, (float) t));
                t += tInc;
                tInc += inc;
            }

            //last point (in case t overshot zMax)
            t = zMax;
            v = v1 + (t - t1) * delV;
            newVels.add(new StsPoint((float) v, (float) t));

            //put into new array
            StsPoint[] newPoints = new StsPoint[newVels.size()];
            for (int i = 0; i < newPoints.length; i++)
            {
                newPoints[i] = newVels.get(i);
            }

            return newPoints;
        }

        private void initializeOrientation()
        {
            row = 0;
            col = 0;
            rowInc = 0;
            colInc = 0;
            StsCVSProperties cvsProps = lineSet.cvsProperties;
            int orientation = cvsProps.getOrientation();

            switch (orientation)
            {
                case StsCVSProperties.ORIENT_INLINE:
                    col = superGatherCol - nTracesPerPanel / 2;
                    row = superGatherRow;
                    rowInc = 0;
                    colInc = 1;
                    break;
                case StsCVSProperties.ORIENT_XLINE:
                    row = superGatherRow - nTracesPerPanel / 2;
                    col = superGatherCol;
                    rowInc = 1;
                    colInc = 0;
                    break;
                case StsCVSProperties.ORIENT_ANGLE_UP:
                    row = superGatherRow - nTracesPerPanel / 2;
                    col = superGatherCol - nTracesPerPanel / 2;
                    rowInc = 1;
                    colInc = 1;
                    break;
                case StsCVSProperties.ORIENT_ANGLE_DOWN:
                    row = superGatherRow + nTracesPerPanel / 2;
                    col = superGatherCol - nTracesPerPanel / 2;
                    rowInc = -1;
                    colInc = 1;
                    break;
            }
        }

        private void initializeSuperGatherDataArrays(StsGather[] gathers)
        {
            int nGathers = gathers.length;
            superGatherData = new double[nGathers][][];
            for (int n = 0; n < nGathers; n++)
            {
                StsGather gather = gathers[n];
                float[] gatherDataf = gather.getGatherDataf();
                if (gatherDataf == null) continue;
                int nGatherTraces = gather.nGatherTraces;
                int nTraceSamples = gather.line.nSlices;
                double[][] gatherData = new double[nGatherTraces][nTraceSamples];
                int i = 0;
                for (int t = 0; t < nGatherTraces; t++)
                    for (int s = 0; s < nTraceSamples; s++)
                        gatherData[t][s] = gatherDataf[i++];
                superGatherData[n] = gatherData;
            }
        }

        /*
                protected float[] computeStackTrace(double velocity)
                {
                    if(!ignoreSuperGather)
                        return computeStackTrace(velocity, superGatherComputer.gathers);
                    else
                        return computeStackTrace(velocity, new StsGather[] { superGatherComputer.centerGather } );
                }
        */
        private float[] computeStackTrace(double v, StsGather[] gathers)
        {
            int nGathers = gathers.length;
            float[][] traces = new float[nGathers][lineSet.nSlices]; //by giving full dimensions, array will be initialized to zeroes which handles dead gathers correctly
            int i = -1, n = -1;
            try
            {
                for (n = 0; n < nGathers; n++)
                {
                    if (!gathers[n].checkInitializeGather()) continue;
                    initializePanelGather(gathers[n]);
                    float trace[] = computeStackTrace(gathers[n], v, n);
                    if (trace == null) continue;
                    traces[n] = trace;
                }

                float[] trace = new float[nTraceSamples];
                System.arraycopy(traces[0], 0, trace, 0, nTraceSamples);
                for (i = 0; i < nTraceSamples; i++)
                {
                    for (n = 1; n < nGathers; n++)
                        trace[i] += traces[n][i];
                }
                trace = StsSeismicFilter.filter(StsFilterProperties.POSTSTACK, trace, 1, nTraceSamples, zInc, lineSet.filterProperties, lineSet.agcProperties, lineSet.dataMin, lineSet.dataMax);
                return trace;
            }
            catch (Exception e)
            {
                StsException.outputWarningException(this, "computeConstantVelocitySuperGatherStackTrace", " slice " + i + " gather " + n, e);
                return null;
            }
        }

        private float[] computeStackTraceFast(float[] traces, double[] offsets, StsPoint[] velPoints)
        {
            if (traces == null || offsets == null || velPoints == null) return null;
            int nSamples = centerGather.line.nSlices;
            if (nSamples == 0) return null;
            float[] stackedData = new float[nSamples];
            int nTraces = offsets.length;

            //... Now we can actually stack data
            if (StsSuperGather.debugTimer) StsSuperGather.StackingTimer.start();

            if (StsSuperGather.debugTimer) StsSuperGather.ApplyNMOTimer.start();
            float[] nmoedTraces = traces.clone();
            applyNMO2(traces, nmoedTraces, offsets, velPoints);
            //applyNMO(nmoedTraces, offsets, vels);
            if (StsSuperGather.debugTimer) StsSuperGather.ApplyNMOTimer.stopAccumulate();
            if (StsSuperGather.debugTimer) StsSuperGather.CalcMeanTimer.start();
            StsTraceUtilities.stack(nmoedTraces, nTraces, stackedData);
            if (StsSuperGather.debugTimer) StsSuperGather.CalcMeanTimer.stopAccumulate();

            if (StsSuperGather.debugTimer) StsSuperGather.StackingTimer.stopAccumulate();

            return stackedData;
        }

        protected void initializePanelGather(StsGather gather)
        {
            nTraceSamples = gather.line.nSlices;
            zMin = gather.line.zMin;
            zMax = gather.line.zMax;
            zInc = gather.line.zInc;
            nFirstGatherLineTrace = gather.nFirstGatherLineTrace;
            nGatherTraces = gather.nGatherTraces;
            traceOffsets = gather.line.traceOffsets;
        }

        /** panelTraces is the float data for [nPanels][nPanelTraces][nTraceSamples]. If a trace is dead, zero it on each of the panels */
        protected void zeroTraceOnPanels(float[][][] panelTraces, int nPanelTrace, int nPanels, float[] zeroTrace)
        {
            for (int p = 0; p < nPanels; p++)
                panelTraces[p][nPanelTrace] = zeroTrace;
        }
    }

    private class ComputeCVSPanelsProcess extends ComputeStackPanelsProcess
    {
        public ComputeCVSPanelsProcess()
        {
        }

        protected void initializePanel()
        {
            StsCVSProperties cvsProps = lineSet.cvsProperties;
            StsSemblanceRangeProperties rangeProps = lineSet.semblanceRangeProperties;
            superGatherComputer = new StsSuperGather(model, lineSet);
            nPanels = cvsProps.getNumberPanels();
            vMin = rangeProps.velocityMin;
            vInc = cvsProps.getCvsVelocityStep(lineSet);
            nTracesPerPanel = cvsProps.getTracesPerPanel();
            ignoreSuperGather = cvsProps.getIgnoreSuperGather();
            localStackedTraces = new float[nPanels][nTracesPerPanel][];
            zeroTrace = new float[nVolumeSlices];
            progressDialog.setProgressMax(nPanels);
            velocityGather = centerGather;
            velocityGather.checkComputeVelocitiesAndOffsetTimes();
        }

        /**
         * For a CVS panel, stacks are at a constant velocity.
         * For the stretchMute to be applied:
         * From Oz p.439, stretch mute factor smf = dT/dt = T/t where T is offset time and t is tzero time.
         * The stretch mute is a percent stretch with 0% as no stretch (which is at x=0).
         * So stretch mute = (s - 1)*100.
         * For each offset, compute the minimum tzero muteTime tzm. Times less than this are muted.
         * Since tsq = tzsq + xsq/vsq
         * muteTime = F*x where F = 1/(Math.sqrt(s*s - 1)*velocity)
         */
        public float[] computeStackTrace(StsGather gather, double v, int nGather)
        {
            int index;
            double f;
            try
            {
                if (StsPreStackLineSet3d.debugTimer) StsSeismicTimer.stackTimer.start();

                double[][] gatherData = superGatherData[nGather];
                if (gatherData == null) return null;
                if (StsPreStackLineSet3d.debugTimer) StsVolumeConstructorTimer.stackTimer.start();
                float t0;
                double rvsq = 1.0 / (v * v);

                double[] sumValues = new double[nTraceSamples];
                int[] nValues = new int[nTraceSamples];
                float[] stackedData = new float[nTraceSamples];

                int nTrace = nFirstGatherLineTrace;
                boolean isNMOed = lineSet.getIsNMOed();
                StsGatherTrace[] gatherTraces = getCenterGatherTraces();
                for (StsGatherTrace trace : gatherTraces)
                {
                    double x = Math.abs(traceOffsets[nTrace++]);
                    double xovsq = x * x * rvsq;
                    t0 = zMin;
                    double value;
                    double v0, v1;
                    float[] data = trace.data;
                    for (int s = 0; s < nTraceSamples; s++, t0 += zInc)
                    {
                        double tOffset = Math.sqrt(t0 * t0 + xovsq);
                        if (!gather.isInsideMuteRange(trace, tOffset, t0)) continue;
                        if (!isNMOed)
                            f = (tOffset - zMin) / zInc;
                        else
                        {
                            float t0i = (float) velocityGather.inputVelocities.computeTZero(x, tOffset, t0);
                            f = (t0i - zMin) / zInc;
                        }

                        if (f < 0.0 || f > nTraceSamples - 1) continue;

                        index = StsMath.floor(f);
                        if (index < 0)
                        {
                            System.out.println("System error! StsSuperGather.computeStackedTrace.index < 0");
                        }
                        if (index == nTraceSamples - 1)
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
                        nValues[s]++;
                    }
                }
                for (int s = 0; s < nTraceSamples; s++)
                {
                    if (nValues[s] > 0) stackedData[s] = ((float) sumValues[s] / nValues[s]);
                    // maxAmplitude = Math.max(maxAmplitude, Math.abs(stackedData[s]));
                }

                //StsMath.normalizeAmplitude(stackedData, (float) maxAmplitude);
                if (StsPreStackLineSet3d.debugTimer) StsVolumeConstructorTimer.stackTimer.stopAccumulate();
                return stackedData;
            }
            catch (Exception e)
            {
                StsException.outputException("StsSuperGather.computeConstantVelocityStackTrace() failed.", e, StsException.WARNING);
                return null;
            }
        }

        protected void completeProcess(boolean successful)
        {
            if (successful)
                //cvsTraces = StsMath.normalizeAmplitude(localStackedTraces); //without this screen is all black
                cvsTraces = StsMath.normalizeAmplitudeRMS(localStackedTraces);
            cvsThread = null;
            computingStacks = false;
            blockReprocessSemblance = false;
            if (progressDialog != null)
            {
                progressDialog.finished();
                progressDialog.dispose();
                progressDialog = null;
            }
            if (debug) System.out.println(StsToolkit.getSimpleClassname(this) + " completed.");
            model.viewObjectRepaintFamily(parentWindow, lineSet);
        }

        @Override
        public void calcVelPoints(StsPoint[] velPoints1, double vel)
        {
            for (int i = 0; i < velPoints1.length; i++)
            {
                velPoints1[i].v[0] = (float) vel;
                velPoints1[i].v[1] = velocityProfile.getProfilePoints()[i].v[1];
            }
        }
    }

    private class ComputeVVSPanelsProcess extends ComputeStackPanelsProcess
    {

        public ComputeVVSPanelsProcess()
        {
        }

        protected void initializePanel()
        {
            StsCVSProperties cvsProps = lineSet.cvsProperties;

            // checkInitializeSuperGatherGeometry();
            // copy the current superGather; this will be used as the compute engine for all gathers in the velocity stack panel//
            superGatherComputer = new StsSuperGather(model, lineSet);
            velocityGather = centerGather;
            if (velocityGather.velocityProfile == null) return;
            velocityGather.velocityProfile.clearVvsInitialProfilePoints();
            velocityGather.velocityProfile.getSetVvsInitialProfilePoints(); //set initialVVS to previously picked velocity, if exists
            velocityGather.checkComputeVelocities();
            nPanels = cvsProps.numberPanels;
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            vInc = rangeProperties.vvsPercentRangeInc / 100;
            vMin = 1 - (nPanels / 2) * vInc;

            nTracesPerPanel = cvsProps.tracesPerPanel;
            ignoreSuperGather = cvsProps.getIgnoreSuperGather();
            localStackedTraces = new float[nPanels][nTracesPerPanel][];
            zeroTrace = new float[nVolumeSlices];
            progressDialog.setProgressMax(nPanels);
            velocityGather = centerGather;
            velocityGather.checkComputeVelocitiesAndOffsetTimes();
        }

        public float[] computeStackTrace(StsGather gather, double vFactor, int nGather)
        {
            int index;
            double f;
            float[] stackedData = null;

            try
            {
                stackedData = new float[nTraceSamples];
                velocityGather.velocities.nmoCoefficients.setFactor(vFactor);
                gather.checkComputeMuteRange();
                //double[] stretchMuteOffsets = gather.compute2ndOrderStretchMuteOffsets(velocityGather.velocityProfile, vFactor);
                //gather.checkComputeStretchMuteTimes();
                //double[] tableTimes = gather.stretchMuteTimes;
                //double[] focusTimes = gather.computeFocusStretchMuteTimes(velocityGather.velocityProfile, vFactor);
                double[][] gatherData = superGatherData[nGather];
                if (gatherData == null) return stackedData;
                if (StsSuperGather.debugTimer) StsSuperGather.StackingTimer.start();

                float t0 = zMin;
                double[] sumValues = new double[nTraceSamples];
                int[] nValues = new int[nTraceSamples];
                StsGather.Velocities velocities = velocityGather.velocities;
                StsGatherTrace[] gatherTraces = getCenterGatherTraces();
                for (int m = 0; m < nTraceSamples; m++, t0 += zInc)
                {
                    int nTrace = nFirstGatherLineTrace;
                    if (StsSuperGather.debugTimer) StsSuperGather.CalcVelocitiesTimer.start();
                    velocities.setVelocityFromTZero(t0, vFactor);
                    if (StsSuperGather.debugTimer) StsSuperGather.CalcVelocitiesTimer.stopAccumulate();
                    if (StsSuperGather.debugTimer) StsSuperGather.ApplyNMOTimer.start();
                    for (StsGatherTrace trace : gatherTraces)
                    {
                        double x = traceOffsets[nTrace];
                        float[] data = trace.data;
                        // if (StsSuperGather.debugTimer) StsSuperGather.CalcTxTimer.start();
                        double tOffset = velocityGather.computeOffsetTime(x);
                        // if (StsSuperGather.debugTimer) StsSuperGather.CalcTxTimer.stopAccumulate();
                        if (!gather.isInsideMuteRange(trace, tOffset, t0)) continue;
                        double tf;
                        if (!lineSet.isNMOed)
                            tf = (tOffset - zMin) / zInc;
                        else
                        {
                            double tZeroInput = velocityGather.inputVelocities.computeTzero(trace, tOffset);
                            tf = (tZeroInput - zMin) / zInc;
                        }
                        if (tf < 0.0 || tf > nTraceSamples - 1) continue;

                        double value;
                        index = StsMath.floor(tf);
                        if (index == nTraceSamples - 1)
                        {
                            value = data[index];
                        }
                        else
                        {
                            tf = tf - index;
                            double v0 = data[index];
                            double v1 = data[index + 1];
                            value = v0 + tf * (v1 - v0);
                        }
                        sumValues[m] += value;
                        nValues[m]++;
                    }
                    if (StsSuperGather.debugTimer) StsSuperGather.ApplyNMOTimer.stopAccumulate();
                }
                if (StsSuperGather.debugTimer) StsSuperGather.CalcMeanTimer.start();
                for (int m = 0; m < nTraceSamples; m++, t0 += zInc)
                {
                    if (nValues[m] == 0) continue;
                    stackedData[m] = (float) (sumValues[m] / nValues[m]);
                }
                if (StsSuperGather.debugTimer) StsSuperGather.CalcMeanTimer.stopAccumulate();
                //			float debugValue = stackedData[100];
                //			System.out.println("row " + row + " col " + col + " value " + debugValue);
                if (StsSuperGather.debugTimer) StsSuperGather.StackingTimer.stopAccumulate();
                return stackedData;
            }
            catch (Exception e)
            {
                //StsException.outputWarningException(this, "computeStackedTrace", e); //this prints out many statements when velocity profile not defined yet
                return stackedData;
            }
        }

        protected void completeProcess(boolean successful)
        {
            if (successful)
            {
                //vvsTraces = StsMath.normalizeAmplitude(localStackedTraces, lineSet.dataMax); //had to set to one normalization per dataset, otherwise unexpected results as effective gain changed from shot-to-shot
                vvsTraces = StsMath.normalizeAmplitudeRMS(localStackedTraces);
            }
            vvsThread = null;
            computingStacks = false;
            blockReprocessSemblance = false;
            if (progressDialog != null)
            {
                progressDialog.finished();
                progressDialog.dispose();
                progressDialog = null;
            }
            if (debug) System.out.println(StsToolkit.getSimpleClassname(this) + " completed.");
            model.viewObjectRepaintFamily(parentWindow, lineSet);
        }

        @Override
        public void calcVelPoints(StsPoint[] velPoints1, double factor)
        {
            for (int i = 0; i < velPoints1.length; i++)
            {
                velPoints1[i].v[0] = (float) (factor * velocityProfile.getProfilePoints()[i].v[0]);
                velPoints1[i].v[1] = velocityProfile.getProfilePoints()[i].v[1];
            }
        }
    }

    /**
     * returns sorted index (not original sort order index) of trace with minimum absolute offset (unsigned offset)
     *
     * @return
     */
    public int getAbsMinOffsetTrace()
    {
        return this.superGatherTraceSet.minOffsetTraceIndex;
    }

    /**
     * gets values of supergather along NMO curve determined by velocity vrms and zero-offset time t0
     * <p/>
     * if curve goes outside gather bounds or into mute-zone, value is zero
     *
     * @return double array
     */
    public void applyNMO(float[] traces, double[] offsets, StsPoint[] velPoints)
    // public void applyNMO(float[] traces, double[] offsets, double[] vels)
    {
        float zMax = lineSet.zMax;
        float zMin = lineSet.zMin;
        double zInc = lineSet.zInc;
        double x = 0;
        double tx = 0;
        double txData = 0;
        double t0 = 0;
        int index = 0;
        double r = 0;
        double vrms = 0;
        double v0 = 0;
        double v1 = 0;
        int nSamples = centerGather.line.nSlices;
        double nmoShift = 0;
        double shift0 = 0;
        double shift1 = 0;
        double tshift0 = 0;
        double tshift1 = 0;
        int nVelPoints = velPoints.length;

        double[][] nmoShifts = new double[nVelPoints][offsets.length];

        for (int i = 0; i < nVelPoints; i++)
        {
            for (int j = 0; j < offsets.length; j++)
            {
                x = offsets[j];
                vrms = velPoints[i].v[0];
                t0 = velPoints[i].v[1];
                tx = Math.sqrt(t0 * t0 + x * x / (vrms * vrms));
                nmoShifts[i][j] = tx;// - t0;
            }
        }

        int iNmo = 0;

        for (int i = 0; i < nSamples; i++)
        {
            t0 = i * zInc;
            while (iNmo < nVelPoints && t0 > velPoints[iNmo].v[1]) iNmo++;
            for (int j = 0; j < offsets.length; j++)
            {

                if (iNmo == 0) nmoShift = nmoShifts[iNmo][j];
                else if (iNmo < nVelPoints && iNmo > 0)
                {
                    shift0 = nmoShifts[iNmo - 1][j];
                    shift1 = nmoShifts[iNmo][j];
                    tshift0 = velPoints[iNmo - 1].v[1];
                    tshift1 = velPoints[iNmo].v[1];
                    nmoShift = shift0 + (t0 - tshift0) * (shift1 - shift0) / (tshift1 - tshift0);
                }
                else
                {
                    shift0 = nmoShifts[nVelPoints - 2][j];
                    shift1 = nmoShifts[nVelPoints - 1][j];
                    tshift0 = velPoints[nVelPoints - 2].v[1];
                    tshift1 = velPoints[nVelPoints - 1].v[1];
                    nmoShift = shift1 + (t0 - tshift1) * (shift1 - shift0) / (tshift1 - tshift0);
                }
                tx = nmoShift;//t0 + nmoShift;

                // tx = Math.sqrt(t0*t0 + offsets[j]*offsets[j]/(vels[i]*vels[i]));
                //  if (StsSuperGather.debugTimer = true) StsSuperGather.CalcTxTimer.stopAccumulate();

                if (tx < zMin || tx > zMax) txData = 0;
                else if (!isInsideStretchMuteRange(tx, t0, x)) txData = 0;
                else
                {
                    //   if (StsSuperGather.debugTimer = true) StsSuperGather.InterpolateTracesTimer.start();
                    //txData = StsTraceUtilities.getInterpolatedValue(traces[i], tx, zInc);
                    index = (int) (tx / zInc); //Math.floor(tx/zInc);
                    if (index < nSamples - 1)
                    {
                        r = tx - index * zInc;
                        v0 = traces[j * nSamples + index];
                        v1 = traces[j * nSamples + index + 1];
                        txData = v0 + r * (v1 - v0) / zInc;
                    }
                    if (index >= nSamples - 1) txData = traces[j * nSamples + nSamples - 1];
                }
                //   if (StsSuperGather.debugTimer = true) StsSuperGather.InterpolateTracesTimer.stopAccumulate();

                traces[j * nSamples + i] = (float) txData;
            }
        }

        //return curveVals;
    }

    /**
     * gets values of supergather along NMO curve determined by velocity vrms and zero-offset time t0
     * <p/>
     * uses super-fast algorithm that interpolates tOffset between user's velocity picks instead of velocity.
     * This greatly simplifies the mathematics for calculating tOffset at each sample along the trace.
     * It is merely a matter of adding deltaTx to the current tx (instead of squaring numbers and calculating square root).
     * Some error is involved in this (particularly at earlier times) - which is why extra picks are added (by interpolation)
     * into the users velocity profile before running this. (run expandVels() on velPoints first!!).
     * <p/>
     * if curve goes outside gather bounds or into mute-zone, value is zero
     *
     * @return double array
     */
    public void applyNMO2(float[] traces, float[] nmoedTraces, double[] offsets, StsPoint[] velPoints)
    {
        float zMax = lineSet.zMax;
        float zMin = lineSet.zMin;
        double zInc = lineSet.zInc;
        double x = 0; //offset of current trace
        double tx = 0; //time at offset (in samples)
        double txData = 0; //value of data at tx (along nmo curve)
        double vp1 = 0; //velocity of 1st pick
        int nSamples = centerGather.line.nSlices;
        double deltaTx = 0; //change in tx per sample between velocity picks
        double t0p1 = 0; //zero-offset time for 1st pick
        double tx1 = 0; //time at offset for 1st pick
        double vp2 = 0; //velocity of 2nd pick
        double t0p2 = 0; //zero-offset time for 2nd pick
        double tx2 = 0; //time at offset for 2nd pick
        double iT0p2 = 0; //zero-offset time for 1st pick (in samples)

        for (int j = 0; j < offsets.length; j++)
        {
            int i = 1;
            x = offsets[j];
            vp1 = velPoints[0].v[0];
            t0p1 = velPoints[0].v[1];  //time of 1st pick in milliseconds;
            tx1 = Math.sqrt(t0p1 * t0p1 + x * x / (vp1 * vp1)) / zInc; //offset time in samples of 1st pick
            for (int t = 0; t < nSamples; t++)
            {
                vp2 = velPoints[i].v[0];
                t0p2 = velPoints[i].v[1]; //time of2nd pick in milliseconds;
                tx2 = Math.sqrt(t0p2 * t0p2 + x * x / (vp2 * vp2)) / zInc; //offset time in samples of 2nd pick;

                deltaTx = (tx2 - tx1) / ((t0p2 - t0p1) / zInc);  //change in tx per sample between pick 1 and 2
                tx = tx1 + (t - t0p1 / zInc) * deltaTx;      //offset time in samples
                iT0p2 = t0p2 / zInc;
                // StsMessage.printMessage("\n i="+i+"\n");
                while (t < iT0p2) //loop through samples between velocity picks
                {
                    // double t0 = t*zInc;
                    // double vel = vp1 + (t0-t0p1)*(vp2-vp1)/(t0p2-t0p1);
                    // double txTrue = Math.sqrt(t0*t0 + x*x/(vel*vel));
                    // StsMessage.printMessage("tx="+tx*zInc+" txTrue="+txTrue+" diff="+(txTrue-tx*zInc));
                    if (tx < zMin || tx > zMax) txData = 0;
                    else if (!isInsideStretchMuteRange(tx, t, x)) txData = 0;
                    else
                    {
                        txData = StsTraceUtilities.getInterpolatedValue(traces, tx, j, nSamples);//, zInc);
                    }
                    nmoedTraces[j * nSamples + t] = (float) txData;
                    tx += deltaTx;
                    t++;
                }
                t--; //we haven't done t=iT0p2 yet
                if (i < velPoints.length - 1)
                {
                    i++;
                    tx1 = tx2;
                    t0p1 = t0p2;
                    vp1 = vp2;
                }
            }
        }
    }

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

    /**
     * Pick/display point on NMO curve
     * <p/>
     * Picking on the gather is this (see StsVelocityProfile.pickOnGather):
     * <p/>
     * picked offset time:                           tp = pick.getY();
     * initial velocity (before editing):            vi = initialPickedPoint.getX();
     * initial t0 (constant while bending residual): t0 = initialPickedPoint.getY();
     * offset distance:                              x = gather.getOffsetForTypeValue(pick.getX());
     * offset time at picked point:                  t = Math.sqrt(t0*t0 + x*x/(vi*vi)) + tp - t0;         equation 1
     * computed velocity for picked point:           v = (float)(Math.abs(x)/Math.sqrt(t*t - t0*t0));     equation 2
     * store new velocity in vertex:                 vertex.setX(v);
     * <p/>
     * Drawing the NMOed curve (see StsGather.displayPickNMOPoint) is this:
     * <p/>
     * velocity at picked point:                    v = vertex.point.getX();
     * offset distance at draw point:               x = drawOffsets[i];
     * initial velocity before picking:             vi = initialPickedPoint.getX();
     * offset time:                                 t = Math.sqrt(t0*t0 + x*x/(v*v));  equation 3
     * flattened offset time:                       tf = Math.sqrt(t0*t0 + x*x/(v*v)) - Math.sqrt(t0*t0 + x*x/(vi*vi) + t0;
     * <p/>
     * equations 1 and 2 should be equivalent to 3, but aren't exactly so.
     */
    public void displayPickNMOPoint(StsPoint point, StsVelocityProfile velocityProfile, int pointIndex, int offsetAxisType, boolean checkMute, boolean pick, GL gl, StsColor color)
    {
        double v = point.v[0]; // current velocity at this point (which possibly has been moved/edited)
        double t0 = point.v[1]; // current tzero at this point (possibly moved/edited)
        double rvsq = 1.0 / (v * v);
        double x = 0.0, t = 0.0;

        boolean flatten = lineSet.lineSetClass.getFlatten();
        boolean isMutePoint = StsVelocityProfile.isMuteIndex(pointIndex);
        boolean nonFlatPoint = velocityProfile.indexPicked == pointIndex || isMutePoint;

        // draw tzero point
        if (pick)
        {
            gl.glInitNames();
            gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX);
            gl.glPushName(pointIndex);
        }

        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_INDEX)
            StsGLDraw.drawPoint2d((float) getAbsMinOffsetTrace(), (float) t0, color, gl, 6);
        else
            StsGLDraw.drawPoint2d(0.0f, (float) t0, color, gl, 6);

        if (pick)
        {
            gl.glPopName();
            gl.glInitNames();
            gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_EDGE);
            gl.glPushName(pointIndex);
        }
        /*
                if(isMutePoint)
                {
                    double v0 = velocityProfile.getVelocityFromTzero((float)t0);
                    rvsq0 = 1.0 / (v0 * v0);
        //                System.out.println("displayPickNMOPoint v: " + v + " initialV: " + v0 + " index: " + pointIndex + " t0: " + t0);
                }
        */

        //      NMOCoefficients nmoCoefficients = getNMOCoefficients(order);
        //        velocities.nmoCoefficients.setFactor(1.0);
        //        velocities.nmoCoefficients.setTZero(t0);
        //       velocities.nmoCoefficients.setIndex(v);
        if (debug) System.out.println("drawing gather point with velocity " + v + " at t0 " + t0);

        if (!checkMute || velocityProfile == null) color.setGLColor(gl);
        gl.glDisable(GL.GL_LIGHTING);
        StsPoint[] profilePoints = velocityProfile.getProfilePoints();
        float velocity = StsMath.interpolateValue(profilePoints, (float) t0, 1, 0);
        double rvsq2 = 1.0 / (velocity * velocity);
        StsGatherTrace[] traces = getSuperGatherTraces(offsetAxisType);
        if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
        {
            gl.glBegin(GL.GL_LINE_STRIP);
            for (StsGatherTrace trace : traces)
            {
                x = trace.x;

                double tOffset = Math.sqrt(t0 * t0 + x * x * rvsq);
                double tOffset2 = Math.sqrt(t0 * t0 + x * x * rvsq2);

                if (checkMute)
                {
                    if (centerGather.isInsideMuteRange(trace, tOffset, t0))
                        color.setGLColor(gl);
                    else
                        outsideMuteColor.setGLColor(gl);
                }
                else
                    color.setGLColor(gl);

                if (flatten)
                {
                    if (nonFlatPoint)
                        //t = velocities.offsetTimes.computeTZeroFromT(i, tOffset, x, t0);
                        t = t0 + (tOffset - tOffset2);
                    else
                        t = t0;
                }
                else
                    t = tOffset;

                gl.glVertex2d(x, t);
            }
            if (debug)
                System.out.println("t0 " + t0 + " x " + x + " c2 " + centerGather.velocities.nmoCoefficients.c2 + " t at offset " + t);
            gl.glEnd();
        }
        else
        {
            //if (absOffsetSortedLineIndices == null) computeSortedTraceIndicesByAbsOffset();
            // double[] traceOffsets = line.getTraceOffsets();
            // if (traceOffsets == null) return;
            gl.glBegin(GL.GL_LINE_STRIP);
            int i = 0;
            for (StsGatherTrace trace : traces)
            {
                x = trace.x;
                //ii = absOffsetSortedLineIndices[i] - nFirstGatherLineTrace;
                double tOffset = Math.sqrt(t0 * t0 + x * x * rvsq);
                double tOffset2 = Math.sqrt(t0 * t0 + x * x * rvsq2);

                if (checkMute)
                {
                    if (centerGather.isInsideStretchMuteRange(tOffset, t0, x))
                        color.setGLColor(gl);
                    else
                        outsideMuteColor.setGLColor(gl);
                }
                else
                    color.setGLColor(gl);

                if (flatten)
                {
                    if (nonFlatPoint)
                        //t = velocities.offsetTimes.computeTZeroFromT(ii, tOffset, x, t0);
                        t = t0 + (tOffset - tOffset2);
                    else
                        t = t0;
                }
                else
                    t = tOffset;

                gl.glVertex2d(i++, t);
            }
            gl.glEnd();
        }
        if (pick) gl.glPopName();
        gl.glEnable(GL.GL_LIGHTING);
    }

    public void displayPickStretchMute(GL gl, StsColor color)
    {
        StsWiggleDisplayProperties wiggleDisplayProperties = lineSet.getWiggleDisplayProperties();
        if (wiggleDisplayProperties == null) return;
        int offsetAxisType = wiggleDisplayProperties.getOffsetAxisType();
        displayPickStretchMute(offsetAxisType, false, gl, color);
    }

    public double getTraceAttribute(String attributeName, int traceNum)
    {
        attributes = getTraceAttributes(attributeName);
        if (attributes == null) return traceNum;
        return attributes[traceNum];
    }

    public double[] getTraceAttributes(String attributeName)
    {
        if (this.attributeName == attributeName) return attributes;
        if (attributeName.equals(StsWiggleDisplayProperties.ATTRIBUTE_NONE)) return null;
        double[] attributes = new double[nSuperGatherTraces];
        int pos = 0;
        for (int i = 0; i < nGathersInSuperGather; i++)
        {
            StsGather gather = gathers[i];
            if (gather == null) continue;
            try
            {
                double[] gatherAtts = gather.line.getAttributeArray(attributeName, gather.gatherRow, gather.gatherCol);
                if (gatherAtts != null)
                    System.arraycopy(gatherAtts, 0, attributes, pos, gatherAtts.length);
                pos += gather.nGatherTraces;
            }
            catch (FileNotFoundException e) { return null; }
        }

        return attributes;
    }

    /**
     * Displays any attribute on ViewGather.
     * If currently selected attribute is not measured in time and
     * therefore is outside time bounds (axisRange), it is rescaled to fit within
     * the trace length.
     *
     * @param glPanel3d
     * @param color
     * @return
     */
    public boolean displayAttributeOnGather(StsGLPanel3d glPanel3d, StsColor color)
    {
        if (attributes == null) return false;

        int offsetAxisType = wiggleProperties.getOffsetAxisType();

        double attScale = 1; //used to scale non-time attributes to the time scale
        double attShift = 0; //used to shift non-time attributes into the time scale range;

        GL gl = glPanel3d.getGL();

        gl.glDisable(GL.GL_LIGHTING);
        color.setGLColor(gl);

        // Check attribute range against axis range
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < attributes.length; i++)
        {
            min = Math.min(min, attributes[i]);
            max = Math.max(max, attributes[i]);
        }

        // If selected attribute outside time range (non-time attribute), rescale
        double zMax = lineSet.zMax;
        double zMin = lineSet.zMin;
        if (max == min)
        {
            attScale = 1;
            attShift = (zMin + zMax) / 2;
        }
        else if ((max > zMax) || (min < zMin))
        {
            //StsMessage.printMessage("rescaling attribute: " + attributeName + " to fit within time axis.");
            double attRange = max - min;
            double timeRange = zMax - zMin;
            attScale = timeRange / attRange;
            attShift = zMin - min * attScale;
        }

        double y = 0;
        gl.glBegin(GL.GL_LINE_STRIP);
        int nTraces = superGatherTraceSet.nTraces;
        for (int n = 0; n < nTraces; n++)
        {
            StsGatherTrace trace = (StsGatherTrace) superGatherTraceSet.traces[n];
            int index = trace.index;
            y = attributes[index] * attScale + attShift;
            if (offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
                gl.glVertex2d(trace.x, y);
            else
                gl.glVertex2d((double) n, y);
        }

        gl.glEnd();
        gl.glEnable(GL.GL_LIGHTING);
        return true;
    }

    public double getSortedTraceOffset(int index)
    {
        StsGatherTrace trace = (StsGatherTrace) superGatherTraceSet.traces[index];
        return trace.x;
    }
}
