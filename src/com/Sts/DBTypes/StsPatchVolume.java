package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SurfaceCurvature.*;
import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.nio.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 14, 2009
 * Time: 10:35:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPatchVolume extends StsSeismicBoundingBox implements StsTreeObjectI, StsSerializable
{
    /** patches created are sorted by row, col, and z in this array */
    public StsPatchGrid[] rowSortedPatchGrids;
    protected StsColorscale colorscale;
    protected StsHistogram histogram;
    protected StsSeismicVolume seismicVolume;
    protected String seismicName;
    transient public int nInterpolatedSlices;
    transient public float interpolatedZInc;
    transient public int interpolatedSliceMin;
    transient public int interpolatedSliceMax;
    transient StsPatchVolumeClass patchVolumeClass;
    transient public StsCroppedBoundingBox croppedBoundingBox;
    /** gridList contains patches which have been completed. gridList patches are added from prevRowGrids when they are disconnected
     *  or when the process is completed and remaining grids are added.  Because an original rowGrid patch maybe multilayered,
     *  a rowGrid patch may generate several gridList patches (one for each layer).
     */
    transient ArrayList<StsPatchGrid> gridList;
    /** rowGrids contains new patches and existing patches from previous row connected to the latest row;
     *  at the completion of the row, these become the previousRowGrids.
     */
    transient HashMap<Integer, StsPatchGrid> rowGrids = null;
    /** prevRowGrids are the active patches in the previousRow; when making a connection to a point in the previous row, we look here for a patch */
    transient HashMap<Integer, StsPatchGrid> prevRowGrids = null;
    /** a temporary array built at initialization with final patches sorted by col, row, and z */
    transient public StsPatchGrid[] colSortedPatchGrids;
    /** total number of points on all patches on this volume */
    transient int nPointsTotal;
    /** window size in wavelengths. If zero, window size is 1/2 wavelength. */
    transient int windowSize;
    /** half window size in wavelengths. */
    transient float halfWindowSize;
    /**
     * Minimum data amplitude that will be used for a minimum or maximum event.
     * Max data amplitude is Min(-dataMin, dataMax). So to be used, the abs(amplitude) >= minDataAmplitude
     */
    transient double minDataAmplitude;
    /** pick point on adjoining trace cannot be more than this many wavelengths away from picked point */
    transient float pickDifWavelengths;
    /** Type(s) of pick events to be correlated: min, max, min&max, all */
    transient byte pickType;
    /** number of points to find in half of window */
    transient int nHalfSamples;
    /** Window center is either max or min; window ends are either the same or a zero crossing;  This flag indicates window ends with zero-crossing */
    transient boolean windowEndIsZeroCrossing;
    /** multiplier of half window size yielding the max allowed pick difference */
    transient double halfWindowPickDifFactor;
    /** number of interpolation intervals between samples */
    transient int nInterpolationIntervals;
    /** indicates iterative stretchCorrel is to be applied from max down to min by inc */
    transient boolean isIterative;
    /** max value for stretchCorrelation */
    transient float autoCorMax;
    /** min value for stretchCorrelation */
    transient float autoCorMin;
    /** increment of stretch correlation in interative pick */
    transient float autoCorInc;
    /** manual picking minimum acceptable cross-correlation */
    transient float manualCorMin;
    /** sequence of stretchCorrelations: 1 if not iterative, max to min by inc if iterative */
    transient float[] stretchCorrelations;
    /** number of stretchCorrelations in sequence */
    transient int nIterations;

    /** row currently being computed: used for debug print out only */
    transient int row;
    transient int nPatchPointsMin;
    transient public byte curveType = CURVPos;
    transient public int filterSize = 0;
    transient protected boolean displaySurfs = false;
    transient protected boolean displayVoxels = false;
    transient int nSmallGridsRemoved;
    transient int nParentGrids;
    transient StsPoint currentCursorPoint;
    transient StsPatchGrid cursorPointPatch;
    transient private StsPatchGrid[] selectedPatchGrids;

    /** if the row or col correl is < this value, that link and the associated edge is not drawn */
    transient float minLinkCorrel = 0.0f;

    transient float[] histogramValues;
    transient int nHistogramValues = 10000;

    static protected StsObjectPanel objectPanel = null;

    public static final byte PICK_MAX = 0;
    public static final byte PICK_MIN = 1;
    public static final byte PICK_MIN_MAX = 2;
    public static final byte PICK_ALL = 3;

    public static final byte WINDOW_CENTERED = 0;
    public static final byte WINDOW_ABOVE = -1;
    public static final byte WINDOW_BELOW = 1;

    public static final String[] pickTypeNames = new String[]{"Maximum", "Minimum", "Min+Max", "All"}; //, "Zero-crossing+", "Zero-crossing-", "All"};
    public static final String[] stopCriteriaNames = new String[]{"Stop", "Replace", "Stop if same Z"};
    //static StsComboBoxFieldBean displayAttributeBean = new StsComboBoxFieldBean(StsPatchVolume.class, "displayAttribute", "Attribute");
    static public final float badCurvature = StsPatchGrid.badCurvature;

    public String getSeismicName()
    {
        return seismicName;
    }

    public void setSeismicName(String seismicName)
    {
        this.seismicName = seismicName;
    }

    static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsPatchVolume.class, "colorscale");


    static public final StsFieldBean[] displayFields =
    {
        new StsBooleanFieldBean(StsPatchVolume.class, "isVisible", "Display on Cursors"),
        new StsBooleanFieldBean(StsPatchVolume.class, "displaySurfs", "Display as Surfaces"),
        new StsBooleanFieldBean(StsPatchVolume.class, "displayVoxels", "Display as Voxel cloud"),
    };

    static public final StsFieldBean[] propertyFields = new StsFieldBean[]
    {
        new StsStringFieldBean(StsPatchVolume.class, "name", true, "Name"),
        colorscaleBean
    };

    final public float getZ(int slice)
    {
        return zMin + this.interpolatedZInc*slice;
    }

    public float getDataMin()
    {
        return dataMin;
    }

    public float getDataMax()
    {
        return dataMax;
    }

    static public final int defaultTestWindow = 21;
    static public final float defaultTestWavelengths = 1.0f;

    //Curvature Attribute Types
    static public final byte CURVDip = StsSurfaceCurvatureAttribute.CURVDip;
    static public final byte CURVStrike = StsSurfaceCurvatureAttribute.CURVStrike;
    static public final byte CURVMean = StsSurfaceCurvatureAttribute.CURVMean;
    static public final byte CURVGauss = StsSurfaceCurvatureAttribute.CURVGauss;
    static public final byte CURVPos = StsSurfaceCurvatureAttribute.CURVPos;
    static public final byte CURVNeg = StsSurfaceCurvatureAttribute.CURVNeg;
    static public final byte CURVMin = StsSurfaceCurvatureAttribute.CURVMin;
    static public final byte CURVMax = StsSurfaceCurvatureAttribute.CURVMax;


    static final float largeFloat = StsParameters.largeFloat;

    static final boolean debug = true;
    static final boolean runTimer = false;
    static StsTimer timer;

    private static final long serialVersionUID = 1L;

    public StsPatchVolume()
    {
    }

    public StsPatchVolume(StsSeismicVolume seismicVolume)
    {
        super(false);
        StsToolkit.copySubToSuperclass(seismicVolume, this, StsRotatedGridBoundingBox.class, StsBoundingBox.class, true);
        this.seismicVolume = seismicVolume;
        seismicName = seismicVolume.getName();
        zDomain = seismicVolume.zDomain;
        stsDirectory = seismicVolume.stsDirectory;
        initialize(currentModel);
    }

    public boolean initialize(StsModel model)
    {
        initializeColorscale();
        initializeSliceInterpolation();
        initializePatchPointTotal();
        return true;
    }

    private void initializeSliceInterpolation()
    {
        nInterpolationIntervals = StsTraceUtilities.computeInterpolationInterval(zInc, 5);
        interpolatedZInc = zInc / nInterpolationIntervals;
        interpolatedSliceMin = 0;
        interpolatedSliceMax = (nSlices - 1) * nInterpolationIntervals;
        nInterpolatedSlices = interpolatedSliceMax + 1;
    }

    public void initialize()
    {
        clearSelectedPatches();
    }

    public void initializeColorscale()
    {
        try
        {
            if(colorscale == null)
            {
                StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
                colorscale = new StsColorscale("Curvature", spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW), dataMin, dataMax);
                colorscale.setEditRange(dataMin, dataMax);
                colorscale.addActionListener(this);
            }
            colorscale.setRange(dataMin, dataMax);
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatchcVolume.initializeColorscale() failed.", e, StsException.WARNING);
        }
    }

    public void actionPerformed(ActionEvent e)
    {
        if(e.getSource() instanceof StsColorscale)
        {
            int modifiers = e.getModifiers();
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
        else
        {
            //fireActionPerformed(e);
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
        return;
    }

    public void constructPatchVolume(StsPatchPickPanel pickPanel)
    {
        StsProgressBar progressPanel = pickPanel.progressBar;
        windowSize = pickPanel.corWavelength;
        pickDifWavelengths = pickPanel.maxPickDif;
        float minAmpFraction = pickPanel.minAmpFraction;
        float maxStretch = pickPanel.maxStretch;
        pickType = pickPanel.pickType;
        nPatchPointsMin = pickPanel.minPatchSize;
        isIterative = pickPanel.isIterative;
        autoCorMax = pickPanel.autoCorMax;
        autoCorMin = pickPanel.autoCorMin;
        autoCorInc = pickPanel.autoCorInc;
        manualCorMin = pickPanel.manualCorMin;
        gridList = new ArrayList<StsPatchGrid>(100);
        rowGrids = new HashMap<Integer, StsPatchGrid>();
        prevRowGrids = null;

        if(!isIterative)
        {
            nIterations = 1;
            stretchCorrelations = new float[]{manualCorMin};
        }
        else
        {
            nIterations = StsMath.ceiling(1 + (autoCorMax - autoCorMin) / autoCorInc);
            stretchCorrelations = new float[nIterations];
            float stretchCorrelation = autoCorMax;
            for(int n = 0; n < nIterations; n++, stretchCorrelation -= autoCorInc)
                stretchCorrelations[n] = stretchCorrelation;
        }

        rowSortedPatchGrids = new StsPatchGrid[0];
        colSortedPatchGrids = null;
        initialize();
        StsPatchGrid.staticInitialize(seismicVolume.nRows, seismicVolume.nCols);

        if(progressPanel != null)
            progressPanel.initialize(croppedBoundingBox.nRows);

        TracePoints[] rowTracePoints = null; //new TracePoints[nCols];
        float absSeismicDataMax = Math.min(Math.abs(seismicVolume.dataMin), Math.abs(seismicVolume.dataMax));
        minDataAmplitude = minAmpFraction * absSeismicDataMax;
        initializeParameters();
        initializeToBoundingBox(croppedBoundingBox);
        initializeSliceInterpolation();

        int croppedRowMin = croppedBoundingBox.rowMin;
        int croppedColMin = croppedBoundingBox.colMin;
        int croppedSliceMin = croppedBoundingBox.sliceMin;
        int nCroppedSlices = croppedBoundingBox.sliceMax - croppedSliceMin;
        int volRow = croppedRowMin;
        int nVolSlices = seismicVolume.nSlices;
        float[] traceValues = new float[nVolSlices];
        for(row = 0; row < nRows; row++, volRow++)
        {
            //statusArea.setProgress(row*40.f/nRows);
            TracePoints[] prevRowTracePoints = rowTracePoints;
            rowTracePoints = new TracePoints[nCols];
            TracePoints prevTracePoints = null;
            prevRowGrids = rowGrids;
            rowGrids = new HashMap<>();
            FloatBuffer rowFloatBuffer = seismicVolume.getRowPlaneFloatBuffer(volRow, croppedColMin);
            if(rowFloatBuffer == null) return;
            // if(croppedColMin > 0) rowFloatBuffer.position(croppedColMin * nVolSlices);
            int volCol = croppedColMin;
            for(int col = 0; col < nCols; col++, volCol++)
            {
                rowFloatBuffer.position(col*nSlices + croppedSliceMin);
                rowFloatBuffer.get(traceValues);

                TracePoints tracePoints = new TracePoints(row, col, traceValues, pickType, nSlices, croppedSliceMin, nVolSlices);
                rowTracePoints[col] = tracePoints;
                //checkAddTracePatches(tracePoints, prevRowTracePoints[col], prevTracePoints);
                addTracePatches(tracePoints, prevTracePoints, prevRowTracePoints);
                prevTracePoints = tracePoints;
            }
            processPrevRowGrids(row);
            if(progressPanel == null) continue;
            if(progressPanel.isCanceled())
            {
                progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
                return;
            }
            progressPanel.setValue(row + 1);

        }
        addGridSet(rowGrids);
        finish();
        getPatchVolumeClass().setDisplayCurvature(true);
        progressPanel.finished();
    }

    private void initializeParameters()
    {
        //maxStretchLimit = 1 + maxStretch;
        //minStretchLimit = 1/(1 + maxStretch);

        if(windowSize == 0) // windowSize is 0.5 wavelengths, half-window is 0.25 wavelengths
        {
            nHalfSamples = 1;
            windowEndIsZeroCrossing = true;
            halfWindowSize = 0.25f;
            halfWindowPickDifFactor = pickDifWavelengths / 0.25f;
        }
        else
        {
            boolean isHalfWave = !StsMath.isEven(windowSize);
            if(isHalfWave)
            {
                // window size is odd, so window ends with zero-crossing; half-window size is windowSize/2.
                // we need to find (windowSize +1)/2 zero-crossings above and below window center (which is a max or min).
                nHalfSamples = (windowSize + 1) / 2;
                windowEndIsZeroCrossing = true;
                halfWindowPickDifFactor = pickDifWavelengths * 2 / windowSize;
            }
            else
            {
                // window size is even, so window ends with same point type as center; half-window size is windowSize/2.
                // we need to find windowSize/2 points above and below with same point type as window center (which is a max or min).
                nHalfSamples = windowSize / 2;
                halfWindowPickDifFactor = pickDifWavelengths / nHalfSamples;
                windowEndIsZeroCrossing = false;
            }
        }
    }

    private void addTracePatches(TracePoints trace, TracePoints prevColTrace, TracePoints[] prevRowTraces)
    {
        StsPatchPoint[] newPatchPoints = trace.tracePatchPoints;
        int nNewPatchPoints = newPatchPoints.length;
        if(nNewPatchPoints == 0) return;

        ArrayList<TracePoints> otherTraces = getOtherTraces(trace, prevColTrace, prevRowTraces);
        int nOtherTraces = otherTraces.size();
        if(nOtherTraces == 0) return;
        // For each patchPoint, check previous traces for points which are same type and just above or below this new point
        // and find which of these two possible prev trace points has the best correlation.
        // If this correlation is above the minCorrelation, add this new point to the prev point patch.
        // If the new point already has a patch (because it correlated with one of the other of the 4 otherTraces, then the addPatchPointToPatch
        // will merge the two patches.
        for(int i = 0; i < nIterations; i++)
        {
            initializeTraceIndices(trace, otherTraces);
            float minStretchCorrelation = stretchCorrelations[i];
            for(int centerPointIndex = 0; centerPointIndex < nNewPatchPoints; centerPointIndex++)
            {
                StsPatchPoint centerPoint = newPatchPoints[centerPointIndex];
                if(centerPoint.patchID != -1) continue;
                CorrelationWindow window;
                window = new CorrelationWindow(trace, centerPoint, centerPointIndex);
                if(!window.initialize()) continue;
                int centerSlice = centerPoint.slice;
                //iterate through adjacent traces & find best patch (if exists)
                StsPatchConnection[] connections = new StsPatchConnection[0];
                for(TracePoints otherTrace : otherTraces)
                {
                    // int[] possibleMatchIndices = otherTrace.getPossibleEvents(z0, z1);
                    CorrelationWindow matchingWindow = window.findMatchingWindows(otherTrace, minStretchCorrelation);
                    if(matchingWindow == null) continue;
                    StsPatchPoint otherCenterPoint = matchingWindow.centerPoint;

                    double distance = Math.abs(otherCenterPoint.slice - centerSlice);
                    StsPatchConnection connection = new StsPatchConnection(centerPoint, otherCenterPoint, matchingWindow.stretchCorrelation, distance);
                    connections = (StsPatchConnection[])StsMath.arrayAddElement(connections, connection);
                    otherTrace.centerPointIndex = matchingWindow.centerPointIndex;
                    // otherTrace.centerOffset = otherTrace.centerPointIndex - centerPointIndex;
                }
                addConnections(centerPoint, connections);
                // checkAddConnections(newPoint, newConnections);
            }
        }
    }

    private void checkForMultipleConnections(StsPatchConnection connection, StsPatchConnection[] newConnections)
    {
        if(newConnections.length == 0) return;

    }

    private void initializeTraceIndices(TracePoints prevColTrace, ArrayList<TracePoints> otherTraces)
    {
        if(prevColTrace != null)
            prevColTrace.centerPointIndex = -1;
        if(otherTraces == null) return;
        for(TracePoints otherTrace : otherTraces)
            if(otherTrace != null) otherTrace.centerPointIndex = -1;
    }

    void addConnections(StsPatchPoint newPoint, StsPatchConnection[] newConnections)
    {
        int nNewConnections = newConnections.length;
        if(nNewConnections == 0) return;
        for(StsPatchConnection connection : newConnections)
            addPatchConnection(connection);
    }

    public StsPatchVolumeClass getPatchVolumeClass()
    {
        if(patchVolumeClass != null) return patchVolumeClass;
        patchVolumeClass = (StsPatchVolumeClass)getCreateStsClass();
        return patchVolumeClass;
    }

    void addPatchConnections(StsPatchConnection[] newConnections)
    {
        for(StsPatchConnection connection : newConnections)
            addPatchConnection(connection);
    }

    boolean valueOK(double v)
    {
        return Math.abs(v) >= minDataAmplitude;
    }

    private void finish()
    {
        // build patchGrid arrays. Overlapped grids will be split out into new grids
        // So construct an array of existing patchGrids and construct the various 2D arrays for each.
        // New grids generated will be added to the gridList
        int nGrids = gridList.size();
        rowSortedPatchGrids = new StsPatchGrid[nGrids];
        gridList.toArray(rowSortedPatchGrids);
        debugCheckEmptyFraction(rowSortedPatchGrids);
        StsException.systemDebug(this, "finish", " Number of parent grids: " + nParentGrids + "Number of child grids: " + nGrids + " too small: " + nSmallGridsRemoved);
        StsException.systemDebug(this, "finish", "max grid dimension: " + StsPatchGrid.maxGridSize);

        StsPatchGrid.sortRowFirst = true;
        Arrays.sort(rowSortedPatchGrids);
        // compute the total number of points
        initializePatchPointTotal();
        StsException.systemDebug(this, "finish", "Total number of points: " + nPointsTotal);
       // reset the index of each patch to its sequence in row-ordered array
        nGrids = 0;
        for(StsPatchGrid grid : rowSortedPatchGrids)
            grid.resetIndex(nGrids++);

        int num = getPatchVolumeClass().getSize();
        String newname = seismicName + ".patchVolume" + (num + 1);
        setName(newname);
        clearConstructionArrays();
        if(!isPersistent()) addToModel();
        getPatchVolumeClass().setIsVisible(true);
        setIsVisible(true);
        if(StsPatchGrid.nOverlappedPoints > 0) StsException.systemError("Error: number of overlapping points for all grids > 0: " + StsPatchGrid.nOverlappedPoints);
        if(StsPatchGrid.debug) StsPatchGrid.printOverlapPoints();
    }

    private void debugCheckEmptyFraction(StsPatchGrid[] grids)
    {
        float nTotal = 0;
        int[] used = null;
        float nUsed = 0;
        float nActualUsed = 0;

        for(StsPatchGrid grid : grids)
        {
            nTotal += grid.getGridSize();
            used  = grid.getGridPointsUsed();
            nUsed += used[0];
            nActualUsed += used[1];
        }
        float fractionUsed = nUsed/nTotal;
        float fractionActualUsed = nActualUsed/nTotal;
        StsException.systemDebug(this, "debugCheckEmptyFraction", "Fraction used: " + fractionUsed + "Fraction actualUsed: " + fractionActualUsed);
    }

    private void clearConstructionArrays()
    {
        gridList = null;
        rowGrids = null;
        prevRowGrids = null;
    }

    private void initializePatchPointTotal()
    {
        nPointsTotal = 0;
        if(rowSortedPatchGrids == null) return;
        for(StsPatchGrid patchGrid : rowSortedPatchGrids)
            nPointsTotal += patchGrid.nPatchPoints;
    }

    private boolean checkColSortedPatchGrids()
    {
        if(colSortedPatchGrids != null) return true;
        if(rowSortedPatchGrids == null) return false;
        int nGrids = rowSortedPatchGrids.length;
        if(nGrids == 0) return false;
        colSortedPatchGrids = new StsPatchGrid[nGrids];
        System.arraycopy(rowSortedPatchGrids, 0, colSortedPatchGrids, 0, nGrids);
        StsPatchGrid.sortRowFirst = false;
        Arrays.sort(colSortedPatchGrids);
        return true;
    }

    private boolean checkRowSortedPatchGrids()
    {
        return rowSortedPatchGrids != null;
    }

    void addGrid(StsPatchGrid patchGrid)
    {
        gridList.add(patchGrid);
    }

    StsPatchGrid getGrid(int target)
    {
        int number = gridList.size();
        int high = number, low = -1, probe;
        while(high - low > 1)
        {
            probe = (high + low) / 2;
            int id = gridList.get(probe).idFinal;
            if(id > target)
                high = probe;
            else
                low = probe;
        }
        if(low == -1 || gridList.get(low).idFinal != target)
            return null;
        else
            return gridList.get(low);
    }

    /** Called for the last row only; unconditionally add these grids to volume unless they are too small. */
    void addGridSet(HashMap<Integer, StsPatchGrid> gridSetHashMap)
    {
        StsPatchGrid[] patchGrids = gridSetHashMap.values().toArray(new StsPatchGrid[0]);
        for(int  n = 0; n < patchGrids.length; n++)
        {
            StsPatchGrid patchGrid = patchGrids[n];
            patchGrid.initializeGrid();
            if(!patchGrid.isTooSmall(nPatchPointsMin))
                patchGrid.constructGridArray(gridList);
        }
        if(StsPatchVolume.debug)
            StsException.systemDebug(this, "addToGrids", "added " + patchGrids.length + " grids from the last row.");
    }

    public void setCroppedBoundingBox(StsCroppedBoundingBox croppedBoundingBox)
    {
        this.croppedBoundingBox = croppedBoundingBox;
        croppedBoundingBox.setCroppedBoxRange();
    }

    /*
    * run the curvature calculation on the patches
    */
    public void runCurvature(StsProgressPanel progressPanel, int filterSize, byte curveType, boolean runAllPatches)
    {
        this.filterSize = filterSize;

        if(runTimer)
        {
            timer = new StsTimer();
            timer.start();
        }
        StsPatchGrid[] runPatches = getRunCurvaturePatches(runAllPatches);
        int numPatches = runPatches.length;
        if(progressPanel != null)
        progressPanel.initialize(numPatches);
        int nValuePoints = 0;
        int nPoints = 0;
        float[] values = new float[nPointsTotal];
        double sum = 0;
        int progressUpdateInterval = Math.max(numPatches/200, 1);
        int minNPoints = Math.min(filterSize*filterSize, StsQuadraticCurvature.minNPoints);
        for(int i = 0; i < numPatches; i++)
        {
            StsPatchGrid patch = runPatches[i];
            if(patch == null) continue;
            nPoints += patch.nPatchPoints;
            if(patch.computeCurvature(xInc, yInc, curveType, filterSize, minNPoints))
            {
                for(int row = patch.rowMin; row <= patch.rowMax; row++)
                {
                    for(int col = patch.colMin; col <= patch.colMax; col++)
                    {
                        int ptCol = col - patch.colMin;
                        int ptRow = row - patch.rowMin;
                        float value = patch.curvature[ptRow][ptCol];
                        if (value == StsPatchVolume.nullValue) continue;
                        if(value == badCurvature || value == -badCurvature) continue;
                        values[nValuePoints++] = value;
                        sum += value;
                    }
                }
            }

            if(progressPanel == null) continue;
            if(progressPanel.isCanceled())
            {
                progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
                clearPatches();
                return;
            }
            if(i%progressUpdateInterval == 0) progressPanel.setValue(i);
        }

        values = (float[])StsMath.trimArray(values, nValuePoints);

        StsMessageFiles.infoMessage("Number curvature points: " + nValuePoints + " number of patch points " + nPointsTotal);
        if(debug)
            StsException.systemDebug(this, "runCurvature", "Number curvature points: " + nValuePoints + " number of patch points " + nPointsTotal);


        double mean = sum/nValuePoints;
        int histogramDataInc = Math.max(1, nValuePoints/nHistogramValues);
        nHistogramValues = nValuePoints/histogramDataInc;
        histogramValues = new float[nHistogramValues+1];
        double avgDev = 0;
        double sumSqr = 0;
        nHistogramValues = 0;
        for(int n = 0; n < nValuePoints; n++)
        {
            float value = values[n];
            double dif = value - mean;
            avgDev += Math.abs(dif);
            sumSqr += dif*dif;
            if(n%histogramDataInc == 0)
                histogramValues[nHistogramValues++] = value;
        }
        avgDev = avgDev/nValuePoints;
        double variance = (sumSqr) / nValuePoints;
        double stdDev = Math.sqrt(variance);
        StsMessageFiles.infoMessage("mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
        if(debug) StsException.systemDebug(this, "runCurvature", "mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
        dataMin = (float)(mean - 2.0*avgDev);
        dataMax = (float)(mean + 2.0*avgDev);
        colorscale.setRange(dataMin, dataMax);
        StsMessageFiles.infoMessage("range set to +- 2.0*avg dev: " + dataMin + " to " + dataMax);
        if(debug) StsException.systemDebug(this, "runCurvature", "range set to +- 2.0*std dev: " + dataMin + " to " + dataMax);

//        // reset outliers to dataMin/dataMax
//        for (StsPatchGrid patch : rowSortedPatchGrids)
//        {
//            if (patch == null) continue;
//            for (int row = patch.rowMin; row <= patch.rowMax; row++)
//            {
//                for (int col = patch.colMin; col <= patch.colMax; col++)
//                {
//                    float value = patch.getCurvature(row, col, dataMin, dataMax);
//                    if (value == StsPatchVolume.nullValue) continue;
//                }
//            }
//        }
        calculateHistogram(histogramValues, nHistogramValues);
        progressPanel.finished();

        if (runTimer) timer.stopPrint("Time to compute curvature for " + numPatches + " patches.");
    }

    private void clearPatches()
    {
        for(StsPatchGrid patch : rowSortedPatchGrids)
        {
            if(patch == null) patch.clear();
        }
    }

    private StsPatchGrid[] getRunCurvaturePatches(boolean runAllPatches)
    {
        if(runAllPatches || this.selectedPatchGrids == null) return rowSortedPatchGrids;
        else return selectedPatchGrids;
    }

    private ArrayList<TracePoints> getOtherTraces(TracePoints newTrace, TracePoints prevColTrace, TracePoints[] prevRowTraces)
    {
        ArrayList<TracePoints> otherTraces = new ArrayList<>();
        if(prevColTrace != null)
        {
            addOtherTrace(otherTraces, prevColTrace);
        }
        int col = newTrace.col;
        if(prevRowTraces != null)
        {
            // if (col > colMin)
            //    addOtherTrace(otherTraces, prevRowTraces[col - 1]);
            addOtherTrace(otherTraces, prevRowTraces[col]);
            // if (col < colMax)
            //    addOtherTrace(otherTraces, prevRowTraces[col + 1]);
        }
        return otherTraces;
    }

    private void addOtherTrace(ArrayList<TracePoints> otherTraces, TracePoints otherTrace)
    {
        if(otherTrace.nTracePatchPoints == 0) return;
        otherTraces.add(otherTrace);
    }

    //if(debug) StsException.systemDebug(this, "addCorrelatedTracePatches", "Cor: " + bestCor + " with " + minDistance + " distance");

    // Would avoid an iterator if iterator is simple; here instead use: for(PatchCheck patchCheck : possiblePatches) instead/  TJL 5/23/09
    // Loop below doesn't eliminate two previous traces who disagree on which patch new event should connect to.
    // I suspect a "cracked" patch is very common as we discussed and we need to handle accordingly.
    // There are countless possibilities here, so we need to be very shrewd in selecting the one we think is optimal.
    // I think our distance rejection should be based on wavelength.  As we have discussed, if we are chasing a max,
    // the only possible candidates on previous traces would be the one above and the one below.  If our new trace event
    // is exactly between 2 maxes on the prevTrace, then the distance would be 90 degrees.  I think this is being generous
    // and in practice we might find 30 to 45 to be a better range. Regardless, it would be user definable.
    // When examining the 4 previous traces, if the traces that have a candidate all agree that it's the same patch we can
    // then move to the correlation calculation.  As long as one prevTrace passes the correlation test, we should add the new
    // event to this patch.  If the previous traces don't agree on the patch, then we apply the correlation test.  If the ones
    // that pass still don't agree on the patch we can merge the patches if they don't have any overlap.  It doesn't appear
    // that the algorithm currently would do any merging which I believe is an essential operation.
    // If there is overlap,  we then pick the best one; I would opt for the closest since I think any that
    // pass the correlation test are "equal" in that regard.
    // In order to decide if two or more patches are mergeable, the StsPatchGrid should use a TreeMap for the patchPoint entries
    // instead of an ArrayList.  The key would be the surface index: row*nCols + col.  For two patches, we need to run thru them
    // using the two TreeMaps to determine if there is a common entry.  If there is, then they overlap and can't be merged.
    // A more efficient structure for this would be a simplification of the StsEdgeLoopRadialLinkGrid class which uses row and col links to
    // define a sparse collection of points on a regular grid.
    // TJL 6/23/09

/*
    private float computeCrossCorrelation(CorrelationWindow newWindow, CorrelationWindow otherWindow, byte pointType)
    {
        if(newWindow == null || otherWindow == null) return 0.0f;

        TracePoints traceNew = newWindow.trace;
        int centerNew = newWindow.centerSlice;
        int minNew = newWindow.minSlice;
        int maxNew = newWindow.maxSlice;

        TracePoints traceOther = otherWindow.trace;
        int centerOther = otherWindow.centerSlice;
        int minOther = otherWindow.minSlice;
        int maxOther = otherWindow.maxSlice;

        int nPointsNew = maxNew - minNew + 1;
        StsPatchPoint[] tracePointsNew = traceNew.tracePatchPoints;
        int nTracePointsNew = tracePointsNew.length;
        int[] centersNew = new int[nPointsNew];
        for(int n = minNew, i = 0; n <= maxNew; n++, i++)
            centersNew[i] = tracePointsNew[n].slice;
        int centerMinNew = centersNew[0];
        int centerMaxNew = centersNew[nPointsNew - 1];

        int nPointsOther = maxOther - minOther + 1;
        StsPatchPoint[] tracePointsOther = traceOther.tracePatchPoints;
        int nTracePointsOther = tracePointsOther.length;
        int centerMinOther = tracePointsOther[minOther].slice;
        int centerMaxOther = tracePointsOther[maxOther].slice;

        int[] centerPointsOther = new int[nPointsOther];
        for(int n = minOther, i = 0; n <= maxOther; n++, i++)
            centerPointsOther[i] = tracePointsOther[n].slice;
        // translate and stretch/shrink pointsOther z values so they line up with pointsNew

        int dzMinusOther = centerOther - centerMinOther;
        int dzMinusNew = centerNew - centerMinNew;
        float dzMinusOtherScalar = (float)dzMinusNew / dzMinusOther;

        int dzPlusOther = centerMaxOther - centerOther;
        int dzPlusNew = centerMaxNew - centerNew;
        double dzPlusOtherScalar = (float)dzPlusNew / dzPlusOther;

        for(int i = 0; i < nPointsOther; i++)
        {
            centerOther = centerPointsOther[i];
            double dz = zOther - centerOther;
            if(dz < 0.0)
                pointsOther[i][0] = centerNew + dz * dzMinusOtherScalar;
            else
                pointsOther[i][0] = centerNew + dz * dzPlusOtherScalar;
        }
        int startIndex = minNew;
        double[][] interpolatedNewPoints = new double[nPointsOther][];
        CorrelationPointPair[] otherPointPairs = new CorrelationPointPair[nPointsOther];
        for(int n = 0; n < nPointsOther; n++)
        {
            double zOther = pointsOther[n][0];
            double[] point = new double[]{zOther, 0, 0};
            double[] pointNextNew = tracePointsNew[startIndex].point;
            while(startIndex < nTracePointsNew - 1)
            {
                double[] pointPrevNew = pointNextNew;
                pointNextNew = tracePointsNew[startIndex + 1].point;

                if(zOther <= pointNextNew[0])
                {
                    point = StsTraceUtilities.interpolateCubicValueAndSlope(pointPrevNew, pointNextNew, zOther, zInc);
                    break;
                }
                startIndex++;
            }
            otherPointPairs[n] = new CorrelationPointPair(point, pointsOther[n]);
        }
        startIndex = minOther;
        double[][] interpolatedOtherPoints = new double[nPointsNew][];
        CorrelationPointPair[] newPointPairs = new CorrelationPointPair[nPointsNew];
        for(int n = 0; n < nPointsNew; n++)
        {
            double zNew = pointsNew[n][0];
            double dzNew = zNew - zCenterNew;
            double zOther;
            if(dzNew < 0.0)
                zOther = zCenterOther + dzNew / dzMinusOtherScalar;
            else
                zOther = zCenterOther + dzNew / dzPlusOtherScalar;
            double[] point = new double[]{zNew, 0, 0};
            double[] pointNextOther = tracePointsOther[startIndex].point;
            while(startIndex < nTracePointsOther - 1)
            {
                double[] pointPrevOther = pointNextOther;
                pointNextOther = tracePointsOther[startIndex + 1].point;
                if(zOther <= pointNextOther[0])
                {
                    point = StsTraceUtilities.interpolateCubicValueAndSlope(pointPrevOther, pointNextOther, zOther, zInc);
                    break;
                }
                startIndex++;
            }
            newPointPairs[n] = new CorrelationPointPair(pointsNew[n], point);
        }
        // first, center, and last otherPointPairs are repeats of first, center, and last newPointPairs, so remove them
        CorrelationPointPair[] newOtherPointPairs = new CorrelationPointPair[nPointsOther - 3];
        int windowCenterOther = centerOther - minOther;
        int nMinusPoints = windowCenterOther - 1;
        int nPlusPoints = maxOther - centerOther - 1;
        System.arraycopy(otherPointPairs, 1, newOtherPointPairs, 0, nMinusPoints);
        System.arraycopy(otherPointPairs, windowCenterOther + 1, newOtherPointPairs, nMinusPoints, nPlusPoints);
        otherPointPairs = newOtherPointPairs;
        CorrelationPointPair[] pointPairs = (CorrelationPointPair[])StsMath.arrayAddArray(newPointPairs, otherPointPairs);
        Arrays.sort(pointPairs);

        int nTotalValues = pointPairs.length;
        double valueNew;
        double valueOther;
        double dz;
        double cor = 0;
        double newCovar = 0;
        double otherCovar = 0;

        for(int n = 1; n < nTotalValues; n++)
        {
            dz = pointPairs[n].z - pointPairs[n - 1].z;
            valueNew = dz * (pointPairs[n - 1].newValue + pointPairs[n].newValue) / 2;
            valueOther = dz * (pointPairs[n - 1].otherValue + pointPairs[n].otherValue) / 2;

            cor += valueNew * valueOther;
            newCovar += valueNew * valueNew;
            otherCovar += valueOther * valueOther;
        }
        cor /= Math.sqrt(newCovar * otherCovar);
        return (float)cor;
    }
*/
    private float computeStretchCorrelation(CorrelationWindow newWindow, CorrelationWindow otherWindow)
    {
        if(newWindow == null || otherWindow == null) return 0.0f;

        TracePoints traceNew = newWindow.trace;
        int centerNew = newWindow.centerSlice;
        int minNew = newWindow.minSlice;
        int maxNew = newWindow.maxSlice;

        TracePoints traceOther = otherWindow.trace;
        int centerOther = otherWindow.centerSlice;
        int minOther = otherWindow.minSlice;
        int maxOther = otherWindow.maxSlice;

        // translate and stretch/shrink pointsOther z values so they line up with pointsNew

        int dSliceMinusOther = centerOther - minOther;
        int dSliceMinusNew = centerNew - minNew;
        float dSliceMinusOtherScalar = (float)dSliceMinusNew / dSliceMinusOther;
        // if(dzMinusOtherScalar < minStretchLimit || dzMinusOtherScalar > maxStretchLimit) return 0.0f;

        float minusStretchFactor = dSliceMinusOtherScalar;
        if(minusStretchFactor > 1.0f)
            minusStretchFactor = 1 / minusStretchFactor;

        int dSlicePlusOther = maxOther - centerOther;
        int dSlicePlusNew = maxNew - centerNew;
        float dSlicePlusOtherScalar = (float)dSlicePlusNew / dSlicePlusOther;
        // if(dzPlusOtherScalar < minStretchLimit || dzPlusOtherScalar > maxStretchLimit) return 0.0f;

        float plusStretchFactor = dSlicePlusOtherScalar;
        if(plusStretchFactor > 1.0f)
            plusStretchFactor = 1 / plusStretchFactor;

        return (float)Math.min(minusStretchFactor, plusStretchFactor);
    }

    class CorrelationPointPair implements Comparable<CorrelationPointPair>
    {
        double z;
        double newValue, otherValue;

        CorrelationPointPair(double[] pointNew, double[] pointOther)
        {
            z = pointNew[0];
            newValue = pointNew[1];
            otherValue = pointOther[1];
        }

        public int compareTo(CorrelationPointPair other)
        {
            double otherZ = other.z;
            if(z < otherZ)
                return -1;
            else if(z > otherZ)
                return 1;
            else
                return 0;
        }
    }
/*
    private double getTraceInterpolatedValue(StsPatchPoint[] patchPoints, double z, int startIndex)
    {
        int nPoints = patchPoints.length;
        while(startIndex < nPoints - 1)
        {
            if(patchPoints[startIndex].point[0] <= z)
                return StsTraceUtilities.interpolateCubic(patchPoints[startIndex].point, patchPoints[startIndex + 1].point, z, zInc);
            startIndex++;
        }
        return StsParameters.doubleNullValue;
    }
*/
/*
    private int[] getOtherCorrelationWindowPoints(TracePoints otherTrace, int nOtherPatchIndex, double zMin, double zMax)
    {
        StsPatchPoint otherPatchPoint = otherTrace.tracePickTypePoints[nOtherPatchIndex];
        byte pointType = otherPatchPoint.pointType;
        int traceIndex = otherPatchPoint.traceIndex;
        StsPatchPoint[] tracePatchPoints = otherTrace.tracePatchPoints;
        int nPoints = tracePatchPoints.length;
        int firstTraceIndex = traceIndex;
        StsPatchPoint patchPoint = otherPatchPoint;
        while(firstTraceIndex > 0 && zMin <= tracePatchPoints[firstTraceIndex-1].point[0])
            firstTraceIndex--;

        int lastTraceIndex = firstTraceIndex;
        while(lastTraceIndex < nPoints-1 && zMax >= tracePatchPoints[lastTraceIndex+1].point[0])
            lastTraceIndex++;

        return new int[] { firstTraceIndex, lastTraceIndex };
    }
*/

    /**
     * Given a newPatchPoint at newRow-newCol, which correlates with a prevPatchPoint at prevRow-prevCol which is possibly part of a patchGrid in the prevPatchGridsSet,
     * combine these two points in the same patch.  The prevPatchPoint may be on the previous col (same row), or previous row (same col).
     * If the previousPatchPoint is not part of an existing patchGrid (prevID == -1), then we will create a new patchGrid and add both points to it.
     * If the previousPatchPoint is part of a patchGrid we will add the newPatchPoint to this patchGrid, unless the newPatchPoint already belongs to another patchGrid
     * (this occurs when we first correlate with the previous column and find one patchGrid and also correlate with the previous row and find a different patchGrid).
     */
    private void addPatchConnection(StsPatchConnection connection)
    {
        StsPatchPoint newPatchPoint = connection.newPoint;
        StsPatchPoint prevPatchPoint = connection.otherPoint;

        StsPatchGrid patchGrid;

        int prevID = prevPatchPoint.patchID;
        int id = newPatchPoint.patchID;
        if(id == -1)
        {
            if(prevID == -1) // prevPatchGrid doesn't exist, so create it
            {
                patchGrid = StsPatchGrid.construct(this, newPatchPoint.pointType);
            }
            else // prevPatchGrid does exist, so get it
            {
                patchGrid = getPatchGrid(prevID, prevPatchPoint);
                if(patchGrid == null) return;
            }
            // add points and connection to this patchGrid
            patchGrid.addInitialPatchConnection(connection);
            // this patchGrid may or may not have already been added; add it, HashMap won't allow value with same key
            putPatchGridInRowList(patchGrid);

            // add newPatchPoint to prevPatchGrid
            // this might be replacing an existing patchGrid in rowGrids with the same one
        }
        else // id != -1 means this point was just added to a patch from prev row and the patchGrid would have been added to the rowGrids array
        {
            if(prevID == -1) // the prevPatchPoint from the prev col (same row) is not assigned to a patch; assign it to this one
            {
                patchGrid = getPatchGrid(id, newPatchPoint);
                if(patchGrid == null) return; // patchGrid not found; systemDebugError was printed
                patchGrid.addInitialPatchConnection(connection);
            }
            else if(prevID == id) // prevPatchPoint on prev row (same col) is already assigned to the same patch: add this connection
            {
                patchGrid = getPatchGrid(id, newPatchPoint);
                if(patchGrid == null) return; // patchGrid not found; systemDebugError was printed
                patchGrid.addInitialPatchConnection(connection);
            }

            else // prevPoint and this point belong to different patches: merge newPatchGrid into prevPatchGrid and add connection
            {
                StsPatchGrid prevPatchGrid = getPatchGrid(prevID, prevPatchPoint);
                if(prevPatchGrid == null) return;
                StsPatchGrid newPatchGrid = getPatchGrid(id, newPatchPoint);
                if(newPatchGrid == null) return;
                mergePatchGrid(prevPatchGrid, newPatchGrid, connection);
            }
        }
    }

    private void mergePatchGrid(StsPatchGrid patchGrid, StsPatchGrid mergedPatchGrid, StsPatchConnection newConnection)
    {
        // just to keep numbering compact, merge larger id grid into smaller id (generally the newPatchGrid into the prevPatchGrid)
        if(mergedPatchGrid.id < patchGrid.id)
        {
            mergePatchGrid(mergedPatchGrid, patchGrid, newConnection);
            return;
        }
        patchGrid.mergePatchGrid(mergedPatchGrid, newConnection);
        removePatchGridFromRowList(mergedPatchGrid);
        putPatchGridInRowList(patchGrid);
    }


    private void putPatchGridInRowList(StsPatchGrid patchGrid)
    {
        int patchID = patchGrid.id;
        if(StsPatchGrid.debugPatchID != -1 && patchID == StsPatchGrid.debugPatchID)
            StsException.systemDebug(this, "putPatchGridInGridList", "patch " + patchID + " added to rowGridsHashMap for row: " + row);
        rowGrids.put(patchID, patchGrid);
    }

    private void removePatchGridFromRowList(StsPatchGrid patchGrid)
    {
        int patchID = patchGrid.id;
        if(StsPatchGrid.debugPatchID != -1 && patchID == StsPatchGrid.debugPatchID)
            StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from rowGridsHashMap for row: " + row);

        StsPatchGrid removedPatchGrid;
        removedPatchGrid = rowGrids.remove(patchID);
//        if(removedPatchGrid == null || removedPatchGrid  != patchGrid)
//            StsException.systemError(this, "removePatchGridFromRowList", "Failed to find patchGrid " + patchID + " to remove from rowGrids");
        removedPatchGrid = prevRowGrids.remove(patchID);
        // if(removedPatchGrid == null || removedPatchGrid  != patchGrid)
        //    StsException.systemError(this, "removePatchGridFromRowList", "Failed to find patchGrid " + patchID + " to remove from prevRowGrids");
    }

    private StsPatchGrid getPatchGrid(StsPatchPoint patchPoint)
    {
        return getPatchGrid(patchPoint.patchID, patchPoint);
    }

    private StsPatchGrid getPatchGrid(int id, StsPatchPoint patchPoint)
    {
        StsPatchGrid patchGrid = rowGrids.get(id);
        if(patchGrid != null)
        {
            if(StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
                StsException.systemDebug(this, "getPatchGrid", "patch grid " + id + " gotten from rowGridsHashMap.");
            return patchGrid;
        }

        if(prevRowGrids != null)
        {
            patchGrid = prevRowGrids.get(id);
            if(patchGrid != null)
            {
                if(StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
                    StsException.systemDebug(this, "getPatchGrid", "patch grid " + id + " gotten from prevRowsGridsHashMap.");
                return patchGrid;
            }
        }
        StsException.systemError(this, "getPatchGrid",  "Couldn't get patchGrid for id " + id +
                " from patchPoint at row " + patchPoint.row + " col " + patchPoint.col);
        return null;
    }

    /**
     * This PatchGridSet is for the row before row just finished.
     * If a grid in this prev row is disconnected (doesn't have same patch in row just finished),
     * then either delete it if it is a small point, or write it out (disconnected so we won't refer to it again).
     * If not disconnected, add it to volume set.
     */
    void processPrevRowGrids(int row)
    {
        StsPatchGrid[] prevRowPatchGrids = prevRowGrids.values().toArray(new StsPatchGrid[0]);
        //int nPrevRowGrids = prevRowPatchGrids.length;
        int nDisconnectedGrids = 0;
        for(StsPatchGrid patchGrid : prevRowPatchGrids)
        {
            boolean disconnected = patchGrid.isDisconnected(row);
            if(!disconnected) continue;
            patchGrid.initializeGrid();
            // if disconnected but too small, don't add
            if(patchGrid.isTooSmall(nPatchPointsMin))
            {
                nSmallGridsRemoved++;
                continue;
            }
            else
                patchGrid.constructGridArray(gridList);
            prevRowGrids.remove(patchGrid.id);
            nDisconnectedGrids++;
        }
        if(debug) StsException.systemDebug(this, "processPrevRowGrids", "row: " + row + " added " + nDisconnectedGrids + " disconnected grids");
    }

    class TracePoints
    {
        byte pickType;
        int row;
        int col;
        int slice;
        StsPatchPoint[] tracePatchPoints = new StsPatchPoint[0];
        int nTracePatchPoints;
        // int centerOffset = 0;  // this is dif between centers for this trace and the master trace in slice increments
        int centerPointIndex = -1; // index of the last correlated connection for this trace with master trace
        // int indexAbove = 0;
        // int indexBelow = 1;

        TracePoints(int row, int col, float[] traceValues, byte pickType, int nSlices, int volSliceMin, int nVolSlices)
        {
            this.pickType = pickType;        // to match StsTraceUtilities.POINT...
            this.row = row;
            this.col = col;
            float[] tracePoints = StsTraceUtilities.computeCubicInterpolatedPoints(traceValues, nInterpolationIntervals);
            float z = zMin;
            if(tracePoints == null) return;
            nTracePatchPoints = tracePoints.length;
            tracePatchPoints = new StsPatchPoint[nTracePatchPoints];
            int i = 0;
            for(int n = 0; n < nTracePatchPoints; n++, z += interpolatedZInc)
            {
                byte tracePointType = StsTraceUtilities.getPointType(tracePoints, n);
                if(StsTraceUtilities.isMaxMinOrZero(tracePointType))
                    tracePatchPoints[i++] = new StsPatchPoint(row, col, n, z, tracePoints[n], tracePointType);
            }
            nTracePatchPoints = i;
            tracePatchPoints = (StsPatchPoint[])StsMath.trimArray(tracePatchPoints, nTracePatchPoints);
        }
/*
        private int[] findNearestPicks(double pickZ, byte pickType)
        {
            boolean foundPickAbove = getPrevTracePoint(pickType, pickZ);
            indexBelow = Math.min(indexAbove + 1, nTracePatchPoints - 1);
            boolean foundPickBelow = getNextTracePoint(pickType, pickZ);
            if(foundPickAbove && foundPickBelow)
                return new int[]{indexAbove, indexBelow};
            else if(foundPickAbove)
                return new int[]{indexAbove};
            else if(foundPickBelow)
                return new int[]{indexBelow};
            else
                return null;
        }
*/
        /*
        private boolean getNextTracePoint(byte pickType, double pickZ)
        {
            int index = indexBelow;
            StsPatchPoint tracePoint = tracePatchPoints[index];
            // current z is above or at pickZ; increase increment until below pickZ;
            // then search for first pickType match
            if(tracePoint.point[0] <= pickZ)
            {
                while(tracePoint.point[0] <= pickZ)
                {
                    index++;
                    if(index >= nTracePatchPoints) return false;
                    tracePoint = tracePatchPoints[index];
                }
                while(true)
                {
                    if(tracePoint.pointType == pickType)
                    {
                        indexBelow = index;
                        return true;
                    }
                    index++;
                    if(index >= nTracePatchPoints) return false;
                    tracePoint = tracePatchPoints[index];
                }
            }
            else // current z is below pickZ; decrement index until current z is <= pickZ searching for pickType match
            {
                boolean found = false;
                while(tracePoint.point[0] > pickZ)
                {
                    if(tracePoint.pointType == pickType)
                    {
                        indexBelow = index;
                        found = true;
                    }
                    index--;
                    if(index < 0) return false;
                    tracePoint = tracePatchPoints[index];
                }
                return found;
            }
        }

        private boolean getPrevTracePoint(byte pickType, double pickZ)
        {
            int index = indexAbove;
            StsPatchPoint tracePoint = tracePatchPoints[index];
            // if current z is above pickZ, increase index  looking for a pickType match until current z is at or below pickZ
            if(tracePoint.point[0] <= pickZ)
            {
                boolean found = false;
                while(tracePoint.point[0] <= pickZ)
                {
                    if(tracePoint.pointType == pickType)
                    {
                        indexAbove = index;
                        found = true;
                    }
                    index++;
                    if(index >= nTracePatchPoints) return false;
                    tracePoint = tracePatchPoints[index];
                }
                return found;
            }
            else // current z is at or below pickZ; decrement index until we current z is just at or above pickZ
            {
                while(tracePoint.point[0] > pickZ)
                {
                    index--;
                    if(index < 0) return false;
                    tracePoint = tracePatchPoints[index];
                }
                // now continue decreasing until we find the first matching pointType
                while(true)
                {
                    if(tracePoint.pointType == pickType)
                    {
                        indexAbove = index;
                        return true;
                    }
                    index--;
                    if(index < 0) return false;
                    tracePoint = tracePatchPoints[index];
                }
            }
        }
        */
        /*
                private CorrelationWindow getCorrelationWindowRange(int nPatchCenterPoint, byte windowType)
                {
                    StsPatchPoint newPatchPoint = tracePatchPoints[nPatchCenterPoint];
                    byte endPointType;
                    if(windowEndIsZeroCrossing)
                        endPointType = StsTraceUtilities.POINT_ZERO_CROSSING;
                    else
                        endPointType = newPatchPoint.pointType;
                    int nAbovePoints, nBelowPoints;
                    if(windowType == WINDOW_CENTERED)
                    {
                        nAbovePoints = nHalfSamples;
                        nBelowPoints = nHalfSamples;
                    }
                    else if(windowType == WINDOW_ABOVE)
                    {
                        nAbovePoints = 2*nHalfSamples;
                        nBelowPoints = 0;
                    }
                    else // WINDOW_BELOW
                    {
                        nAbovePoints = 0;
                        nBelowPoints = 2*nHalfSamples;
                    }
                    double zCenter = newPatchPoint.point[0];
                    double corSliceMin = zCenter;
                    double corSliceMax = zCenter;

                    int nSamePointType = 0;
                    int nFirstPickTypePoint = nPatchCenterPoint;
                    while (nSamePointType < nAbovePoints)
                    {
                        nFirstPickTypePoint--;
                        if (nFirstPickTypePoint < 0) return null;
                        if (matchesEndPointType(tracePatchPoints[nFirstPickTypePoint].pointType, windowEndIsZeroCrossing, endPointType))
                        {
                            nSamePointType++;
                            if (nSamePointType == 1) corSliceMin = tracePatchPoints[nFirstPickTypePoint].point[0];
                        }
                    }
                    nSamePointType = 0;
                    int nLastPickTypePoint = nPatchCenterPoint;
                    while (nSamePointType < nBelowPoints)
                    {
                        nLastPickTypePoint++;
                        if (nLastPickTypePoint >= nTracePatchPoints) return null;
                        if (matchesEndPointType(tracePatchPoints[nLastPickTypePoint].pointType, windowEndIsZeroCrossing, endPointType))
                        {
                            nSamePointType++;
                            if (nSamePointType == 1) corSliceMax = tracePatchPoints[nLastPickTypePoint].point[0];
                        }
                    }
                    int nFirstTracePoint = tracePatchPoints[nFirstPickTypePoint].traceIndex;
                    int nLastTracePoint = tracePatchPoints[nLastPickTypePoint].traceIndex;
                    int nTraceCenterPoint = tracePatchPoints[nPatchCenterPoint].traceIndex;
                    corSliceMin = halfWindowPickDifFactor * (corSliceMin - zCenter) + zCenter;
                    corSliceMax = halfWindowPickDifFactor * (corSliceMax - zCenter) + zCenter;
                    return new CorrelationWindow(this, nFirstTracePoint, nLastTracePoint, nTraceCenterPoint, corSliceMin, corSliceMax);
                }
        */
        boolean matchesEndPointType(byte pointType, boolean windowEndIsZeroCrossing, byte endPointType)
        {
            if(windowEndIsZeroCrossing)
                return pointType == StsTraceUtilities.POINT_PLUS_ZERO_CROSSING || pointType == StsTraceUtilities.POINT_MINUS_ZERO_CROSSING;
            else
                return pointType == endPointType;
        }
    }

    class CorrelationWindow
    {
        TracePoints trace;
        StsPatchPoint centerPoint;
        int centerPointIndex;
        StsPatchPoint pointAbove;
        StsPatchPoint pointBelow;
        int centerSlice;
        int minSlice;
        int maxSlice;
        int dSliceMinus;
        int dSlicePlus;
        double stretchCorrelation;
        byte centerPointType;
        byte abovePointType;
        byte belowPointType;
        byte windowType;

        CorrelationWindow(TracePoints trace, StsPatchPoint centerPoint, int centerPointIndex)
        {
            this.trace = trace;
            this.centerPoint = centerPoint;
            this.centerPointIndex = centerPointIndex;
            this.centerSlice = centerPoint.slice;
            centerPointType = centerPoint.pointType;
            abovePointType = StsTraceUtilities.pointTypesBefore[centerPointType];
            belowPointType = StsTraceUtilities.pointTypesAfter[centerPointType];
        }

        boolean initialize()
        {
            try
            {
                StsPatchPoint[] tracePatchPoints = trace.tracePatchPoints;
                int nTracePatchPoints = tracePatchPoints.length;
                if(centerPointIndex <= 0)
                {
                    if(nTracePatchPoints < 2) return false;
                    windowType = WINDOW_BELOW;
                    pointAbove = centerPoint;
                    pointBelow = getTracePatchPointBelow();
                    maxSlice = pointBelow.slice;
                    dSlicePlus = maxSlice - centerSlice;
                    dSliceMinus = dSlicePlus;
                    minSlice = centerSlice - dSliceMinus;
                }
                else if(centerPointIndex == nTracePatchPoints - 1)
                {
                    if(nTracePatchPoints < 2) return false;
                    windowType = WINDOW_ABOVE;
                    pointBelow = centerPoint;
                    pointAbove = getTracePatchPointAbove();
                    minSlice = pointAbove.slice;
                    dSliceMinus = centerSlice - minSlice;
                    dSlicePlus = dSliceMinus;
                    maxSlice = centerSlice + dSlicePlus;
                }
                else
                {
                    windowType = WINDOW_CENTERED;
                    pointAbove = getTracePatchPointAbove();
                    pointBelow = getTracePatchPointBelow();
                    minSlice = pointAbove.slice;
                    maxSlice = pointBelow.slice;
                    dSliceMinus = centerSlice - minSlice;
                    dSlicePlus = maxSlice - centerSlice;
                }
                return true;
            }
            catch(Exception e)
            {
                StsException.outputFatalException(this, "initialize", e);
                return false;
            }
        }

        private StsPatchPoint getTracePatchPointBelow()
        {
            StsPatchPoint[] tracePatchPoints = trace.tracePatchPoints;
            int nTracePatchPoints = tracePatchPoints.length;

            for(int index = centerPointIndex + 1; index < nTracePatchPoints; index++)
                if(tracePatchPoints[index].pointType == belowPointType)
                    return tracePatchPoints[index];
            return tracePatchPoints[centerPointIndex + 1];
        }

        private StsPatchPoint getTracePatchPointAbove()
        {
            StsPatchPoint[] tracePatchPoints = trace.tracePatchPoints;

            for(int index = centerPointIndex - 1; index >= 0; index--)
                if(tracePatchPoints[index].pointType == abovePointType)
                    return tracePatchPoints[index];
            return tracePatchPoints[centerPointIndex - 1];
        }

        /**
         * On the other trace, starting from the point following the last correlated point,
         * find any and all picks of the center pick type on the other trace which match.
         * check that the above and below pick types match as well
         */
        CorrelationWindow findMatchingWindows(TracePoints otherTrace, float stretchCorrelation)
        {
            int otherTraceCenterPointIndex = otherTrace.centerPointIndex + 1;
            int nOtherTracePoints = otherTrace.nTracePatchPoints;
            CorrelationWindow matchingWindow = null;
            double bestCorrelation = 0.0;
            int centerSlice = centerPoint.slice;
            while(otherTraceCenterPointIndex < nOtherTracePoints)
            {
                StsPatchPoint otherTraceCenterPoint = otherTrace.tracePatchPoints[otherTraceCenterPointIndex];
                int otherCenterSlice = otherTraceCenterPoint.slice;
                if(otherCenterSlice > maxSlice) break;
                if(otherCenterSlice >= minSlice && centerPointType == otherTraceCenterPoint.pointType)
                {
                    CorrelationWindow otherWindow = checkCreateMatchingWindow(otherTrace, otherTraceCenterPoint, otherCenterSlice, otherTraceCenterPointIndex, centerSlice, stretchCorrelation);
                    if(otherWindow != null && otherWindow.stretchCorrelation > bestCorrelation)
                    {
                        matchingWindow = otherWindow;
                        bestCorrelation = matchingWindow.stretchCorrelation;
                    }
                }
                otherTraceCenterPointIndex++;
            }
            return matchingWindow;
        }

        CorrelationWindow checkCreateMatchingWindow(TracePoints otherTrace, StsPatchPoint otherTraceCenterPoint, int otherCenterSlice, int otherTraceCenterPointIndex, int centerSlice, float stretchCorrelation)
        {
            CorrelationWindow otherWindow = new CorrelationWindow(otherTrace, otherTraceCenterPoint, otherTraceCenterPointIndex);
            if(!otherWindow.initialize()) return null;
            if(otherWindow.isZCenterOutsideWindow(centerSlice) || !matches(otherWindow, stretchCorrelation))
                return null;
            else
                return otherWindow;
        }

        /** check the various correlation measures and return the otherWindow if it matches; otherwise return null */
        boolean matches(CorrelationWindow otherWindow, float stretchCorrelation)
        {
            if(otherWindow.windowType != windowType) return false;
            if(otherWindow.pointAbove.pointType != pointAbove.pointType) return false;
            if(otherWindow.pointBelow.pointType != pointBelow.pointType) return false;
            // check correlation stretch
            double minusStretchFactor = 1.0;
            double plusStretchFactor = 1.0;
            if(windowType != WINDOW_BELOW)
            {
                double dzMinusOtherScalar = ((double)dSliceMinus) / otherWindow.dSliceMinus;
                minusStretchFactor = dzMinusOtherScalar;
                if(minusStretchFactor > 1.0f)
                    minusStretchFactor = 1 / minusStretchFactor;
            }
            else if(windowType != WINDOW_ABOVE)
            {
                double dzPlusOtherScalar = (double)dSlicePlus / otherWindow.dSlicePlus;
                plusStretchFactor = dzPlusOtherScalar;
                if(plusStretchFactor > 1.0f)
                    plusStretchFactor = 1 / plusStretchFactor;
            }
            double stretchFactor = Math.min(minusStretchFactor, plusStretchFactor);
            otherWindow.stretchCorrelation = stretchFactor;
            return stretchFactor >= stretchCorrelation;
        }

        boolean isZCenterOutsideWindow(int centerSlice)
        {
            return centerSlice < minSlice || centerSlice > maxSlice;
        }
    }

    public StsColorscale getCurvatureColorscale()
    {
        return colorscale;
    }

    /* Draw any map edges on all 2d sections */
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped,
                               boolean xAxisReversed, boolean yAxisReversed)
    {
        if(!getIsVisible()) return;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;
        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
        StsColor drawColor = StsColor.BLACK; //getStsColor();
        gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
        drawColor.setGLColor(gl);

        if(dirNo == StsCursor3d.XDIR) /* map edge is along a col	*/
        {
            if(!checkColSortedPatchGrids()) return;
            int col = getNearestColCoor(dirCoordinate);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glShadeModel(GL.GL_SMOOTH);

            float x = dirCoordinate;
            int nFirst = -1;
            int n = -1;

            for(StsPatchGrid patchGrid : colSortedPatchGrids)
            {
                if(patchGrid.colMin > col) break;
                n++;
                if(patchGrid.colMax < col) continue;
                patchGrid.drawCol(gl, col, x, yMin, yInc, colorscale, false, displayCurvature);
                if(nFirst == -1) nFirst = n;
            }
            gl.glEnable(GL.GL_LIGHTING);
        }
        else if(dirNo == StsCursor3d.YDIR)
        {
            if(!checkRowSortedPatchGrids()) return;
            int row = getNearestRowCoor(dirCoordinate);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
            gl.glShadeModel(GL.GL_SMOOTH);

            float y = dirCoordinate;
            int nFirst = -1;
            int n = -1;
            for(StsPatchGrid patchGrid : rowSortedPatchGrids)
            {
                if(patchGrid.rowMin > row) break;
                n++;
                if(patchGrid.rowMax < row) continue;

                patchGrid.drawRow(gl, row, y, xMin, xInc, colorscale, false, displayCurvature);
                if(nFirst == -1) nFirst = n;
            }
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    /** Draw any map edges on section */
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate)
    {
        if(!getIsVisible()) return;
        GL gl = glPanel3d.getGL();
        if(gl == null) return;
        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
        try
        {
            if(dirNo == StsCursor3d.ZDIR)
            {
                if(getDisplaySurfs())
                    displayPatchesNearXYZCursors(glPanel3d);
                //        displayPatchesNearZCursor(glPanel3d, dirCoordinate);
                return;
            }
            if(dirNo == StsCursor3d.YDIR)
            {
                if(!checkRowSortedPatchGrids()) return;
                gl.glDisable(GL.GL_LIGHTING);
                gl.glShadeModel(GL.GL_SMOOTH);
                StsColor drawColor = StsColor.BLACK; //getStsColor();
                drawColor.setGLColor(gl);
                gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
                glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
                int row = getNearestRowCoor(dirCoordinate);
                if(row == -1) return;
                float xMin = getXMin();
                float xInc = getXInc();
                for(StsPatchGrid patchGrid :rowSortedPatchGrids)
                {
                    if(patchGrid.rowMin > row) break;
                    if(patchGrid.rowMax < row) continue;
                    patchGrid.drawRow(gl, row, dirCoordinate, xMin, xInc, colorscale, true, displayCurvature);
                }
            }
            else if(dirNo == StsCursor3d.XDIR)
            {
                if(!checkColSortedPatchGrids()) return;
                gl.glDisable(GL.GL_LIGHTING);
                gl.glShadeModel(GL.GL_SMOOTH);
                StsColor drawColor = StsColor.BLACK; //getStsColor();
                drawColor.setGLColor(gl);
                gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
                glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
                int col = getNearestColCoor(dirCoordinate);
                if(col == -1) return;
                float yMin = getYMin();
                float yInc = getYInc();
                for(StsPatchGrid patchGrid :colSortedPatchGrids)
                {
                    if(patchGrid.colMin > col) break;
                    if(patchGrid.colMax < col) continue;
                    patchGrid.drawCol(gl, col, dirCoordinate, yMin, yInc, colorscale, true, displayCurvature);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawOnCursor3d", e);
        }
        finally
        {
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }

    }

    public void pickOnCursor3d(StsGLPanel3d glPanel3d)
    {
        StsCursor3d cursor3d = glPanel3d.window.getCursor3d();
        if(cursor3d == null) return;
        GL gl = glPanel3d.getGL();
        for(int dir = 0; dir < 2; dir++)
        {
            float dirCoordinate = cursor3d.getCurrentDirCoordinate(dir);
            drawOnCursor3d(glPanel3d, dir, dirCoordinate);
        }
    }

    public void display(StsGLPanel glPanel)
    {
        if(!getDisplaySurfs() || selectedPatchGrids == null) return;
        GL gl = glPanel.getGL();

        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();

        initializePatchDraw(gl);
        for(StsPatchGrid patchGrid : selectedPatchGrids)
            patchGrid.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);

        if(getDisplayVoxels())
        {
            displayVoxels(glPanel);
        }
    }

    public void displayVoxelsCursor(StsGLPanel glPanel3d, StsPoint[] points, boolean is3d)
    {
        //System.out.println("Display Voxels");
        GL gl = glPanel3d.getGL();
        if(gl == null) return;
        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
        StsGridPoint point1 = new StsGridPoint(points[0], this);
        StsGridPoint point2 = new StsGridPoint(points[3], this);
        int sameRow = StsGridPoint.getSameRow(point1, point2); // if not -1, this is the row these two points are on
        if(sameRow != -1)
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glColor4f(1.f, 1.f, 1.f, 1.f);

            float y;

            for(StsPatchGrid patchGrid : rowSortedPatchGrids)
            {
                int row1 = sameRow - 4;
                int row2 = sameRow + 4;
                y = yMin + (yInc * row1);
                for(int n = row1; n <= row2; n++)
                {
                    //if (n == 146)
                    //System.out.println("draw vox row "+n+" "+patchGrid.colMin+" "+patchGrid.colMax+" "+y);
                    patchGrid.drawRow(gl, n, y, xMin, xInc, colorscale, is3d, displayCurvature);
                    y += yInc;
                }
            }

            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
            return;
        }

    }

    public void displayVoxels(StsGLPanel glPanel3d)
    {
        //System.out.println("Display Voxels");

        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
            gl.glColor4f(1.f, 1.f, 1.f, 1.f);
            float xMin = getXMin();
            float xInc = getXInc();
            float yMin = getYMin();
            float yInc = getYInc();
            for(StsPatchGrid patchGrid : rowSortedPatchGrids)
                patchGrid.drawRowVox(gl, yMin, yInc, xMin, xInc, colorscale);
            gl.glEnable(GL.GL_LIGHTING);
            return;
        }

    }

    // surfaces near cursor only to keep clutter down
    public void displayPatchesNearZCursor(StsGLPanel glPanel3d, float z)
    {
        //System.out.println("Display Surfaces");
        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        initializePatchDraw(gl);
        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
        for(StsPatchGrid patchGrid : rowSortedPatchGrids)
        {
            if(patchGrid.isPatchGridNearZCursor(z))
                patchGrid.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);
        }
        return;
    }

    private void initializePatchDraw(GL gl)
    {
        gl.glEnable(GL.GL_LIGHTING);
        gl.glShadeModel(GL.GL_SMOOTH);
        gl.glEnable(GL.GL_NORMALIZE);
        gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE);
        // gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL, GL.GL_SEPARATE_SPECULAR_COLOR);
        //gl.glDisable(GL.GL_CULL_FACE);
        //gl.glDisable(GL.GL_POLYGON_STIPPLE);
        gl.glColor4f(1.f, 1.f, 1.f, 1.f);
    }

    public void displayPatchesNearXYZCursors(StsGLPanel glPanel3d)
    {
        StsCursor3d cursor3d = currentModel.win3d.getCursor3d();
        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        float x = cursor3d.getCurrentDirCoordinate(XDIR);
        float y = cursor3d.getCurrentDirCoordinate(YDIR);
        float z = cursor3d.getCurrentDirCoordinate(ZDIR);
        StsPoint cursorPoint = new StsPoint(x, y, z);
        boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
        if(currentCursorPoint != null && currentCursorPoint.equals(cursorPoint) && cursorPointPatch != null)
        {
            drawPatch(cursorPointPatch, displayCurvature, gl);
            return;
        }

        cursorPoint = null;
        currentCursorPoint = null;

        int volumeRow = getNearestRowCoor(y);
        if(volumeRow == -1) return;
        int volumeCol = getNearestColCoor(x);
        if(volumeCol == -1) return;
        int slice = getNearestSliceCoor(z);
        if(slice == -1) return;
        float dzPatch = largeFloat;
        cursorPointPatch = null;
        for(StsPatchGrid patchGrid : rowSortedPatchGrids)
        {
            if(patchGrid == null) continue;
            // if(patchGrid.values == null) continue;
            if(patchGrid.rowMin > volumeRow) break;
            if(patchGrid.rowMax >= volumeRow)
            {
                if(patchGrid.colMin <= volumeCol && patchGrid.colMax >= volumeCol)
                {
                    float dz = patchGrid.getZDistance(volumeRow, volumeCol, z);
                    if(dz < dzPatch)
                    {
                        cursorPointPatch = patchGrid;
                        dzPatch = dz;
                        currentCursorPoint = new StsPoint(x, y, z);
                    }
                }
            }
        }
        if(cursorPointPatch == null) return;

        drawPatch(cursorPointPatch, displayCurvature, gl);

        return;
    }

    private void drawPatch(StsPatchGrid patch, boolean displayCurvature, GL gl)
    {
        initializePatchDraw(gl);
        patch.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);
    }

//	public void setWizard(StsVolumeCurvatureWizard wizard) {
//		this.wizard = wizard;
//	}
//
//	public StsVolumeCurvatureWizard getWizard() {
//		return wizard;
//	}

    public int[] getPatchRangeForRow(int volumeRow)
    {
        int rowMin = -1;
        int rowMax = -1;
        int nPatchGrids = rowSortedPatchGrids.length;
        int n = 0;
        for(n = 0; n < nPatchGrids; n++)
        {
            StsPatchGrid patchGrid = rowSortedPatchGrids[n];
            if(patchGrid == null) continue;
            rowMax = n - 1;
            if(rowMin == -1)
                if(patchGrid.rowMin <= volumeRow && patchGrid.rowMax >= volumeRow) rowMin = n;
            else if(patchGrid.rowMin > volumeRow)
                break;
        }
        if(rowMin == -1)
            return new int[]{0, 0};
        else
            return new int[]{rowMin, rowMax};
    }

    public int[] getPatchRangeForCol(int col)
    {
        int colMin = -1;
        int colMax = -1;
        int nPatchGrids = colSortedPatchGrids.length;
        for(int n = 0; n < nPatchGrids; n++)
        {
            StsPatchGrid patchGrid = rowSortedPatchGrids[n];
            if(patchGrid == null) continue;
            if(colMin == -1)
            {
                if(patchGrid.colMin <= col && patchGrid.colMax >= col) colMin = n;
            }
            else if(patchGrid.colMin > col)
            {
                colMax = n - 1;
                break;
            }
        }
        if(colMin == -1 || colMax == -1)
            return new int[]{0, 0};
        else
            return new int[]{colMin, colMax};
    }

    public boolean getTraceCurvature(int volRow, int volCol, float[] buffer, int[] patchRange)
    {
        try
        {
            Arrays.fill(buffer, nullValue);
            int nPatchMin = patchRange[0];
            int nPatchMax = patchRange[1];
            boolean traceLoaded = false;
            for(int nPatch = nPatchMin; nPatch <= nPatchMax; nPatch++)
            {
                StsPatchGrid patchGrid = rowSortedPatchGrids[nPatch];
                if(patchGrid == null) continue;
                if(patchGrid.curvature == null) continue;
                int patchRow = volRow - patchGrid.rowMin;
                int patchCol = volCol - patchGrid.colMin;
                if(patchRow < 0 || patchCol < 0 || patchRow >= patchGrid.nRows || patchCol >= patchGrid.nCols) continue;
                float[][] pointsZ = patchGrid.getPointsZ();
                if(pointsZ == null) continue;
                float z = pointsZ[patchRow][patchCol];
                if(z == StsParameters.nullValue) continue;
                float val = patchGrid.curvature[patchRow][patchCol];
                if(val == nullValue) continue;
                int slice = getNearestSliceCoor(z);
                buffer[slice] = val;
                traceLoaded = true;
            }
            return traceLoaded;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getTraceCurvature", e);
            return false;
        }
    }

    public int getNearestSliceCoor(float z)
    {
        int slice = Math.round((z - zMin) / interpolatedZInc);
        if(slice < 0 || slice >= nInterpolatedSlices) return -1;
        return slice;
    }

    public int getPatchPointIndex(StsPatchPoint patchPoint)
    {
        return patchPoint.row*nCols + patchPoint.col;
    }

    public String toString() { return name; }

    public void treeObjectSelected()
    {
        getPatchVolumeClass().selected(this);
        currentModel.getGlPanel3d().checkAddView(StsView3d.class);
        currentModel.win3dDisplayAll();
    }


    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsFieldBean[] getDisplayFields()
    {
        //displayAttributeBean.setListItems(displayAttributes);
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return propertyFields;
    }


    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsColorscale getColorscale()
    {
        //setDataHistogram();
        return colorscale;
    }

    public void setColorscale(StsColorscale colorscale)
    {
        this.colorscale = colorscale;
        currentModel.win3dDisplayAll();
    }

    public void setDisplaySurfs(boolean displaySurfs)
    {
        if(this.displaySurfs == displaySurfs)
            return;
        this.displaySurfs = displaySurfs;
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplaySurfs()
    {
        return displaySurfs;
    }

    public void setDisplayVoxels(boolean displayVoxels)
    {
        if(this.displayVoxels == displayVoxels)
            return;
        this.displayVoxels = displayVoxels;
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayVoxels()
    {
        return displayVoxels;
    }

    public void setIsVisible(boolean vis)
    {
        super.setIsVisible(vis);
        currentModel.win3dDisplayAll();
    }

    public void setDataMin(float min)
    {
        dataMin = min;
        if(colorscale == null) return;
        colorscale.setRange(dataMin, dataMax);
    }

    public void setDataMax(float max)
    {
        dataMax = max;
        if(colorscale == null) return;
        colorscale.setRange(dataMin, dataMax);
    }

    public void addRemoveSelectedPatch(StsCursorPoint cursorPoint)
    {
        boolean displayChildPatches = getPatchVolumeClass().getDisplayChildPatches();
        float[] xyz = cursorPoint.point.v;
        int volumeRow = getNearestRowCoor(xyz[1]);
        int volumeCol = getNearestColCoor(xyz[0]);
        float z = xyz[2];
        StsPatchGrid selectedPatch = getNearestPatch(volumeRow, volumeCol, z);
        if(selectedPatch == null) return;
        
        float iline = getRowNumFromRow(volumeRow);
        float xline = getColNumFromCol(volumeCol);
        StsMessageFiles.logMessage("Picked patch parent id: " + selectedPatch.id + " child id: " + selectedPatch.idFinal + " at iline: " + iline + " xline: " + xline);
    /*
        if(cursorPoint.dirNo == StsCursor3d.YDIR)
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeRowCorrel(volumeRow, volumeCol));
        else //dirNo == XDIR
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeColCorrel(volumeRow, volumeCol));
    */
        int parentID = selectedPatch.id;
        // int nSheet = selectedPatch.nSheet;
        String addRemove;
        boolean removePatch = StsMath.arrayContains(selectedPatchGrids, selectedPatch);
        if(removePatch)
            addRemove = " removed ";
        else
            addRemove = " added ";
        int nGridsAdded = 0;
        if(!displayChildPatches)
        {
            if(removePatch)
                selectedPatchGrids = (StsPatchGrid[])StsMath.arrayDeleteElement(selectedPatchGrids, selectedPatch);
            else
                selectedPatchGrids = (StsPatchGrid[])StsMath.arrayAddElement(selectedPatchGrids, selectedPatch);
        }
        else
        {
            for(StsPatchGrid patchGrid : rowSortedPatchGrids)
            {
                if(patchGrid.id == parentID) // && patchGrid.nSheet == nSheet)
                {
                    nGridsAdded++;
                    if(removePatch)
                        selectedPatchGrids = (StsPatchGrid[])StsMath.arrayDeleteElement(selectedPatchGrids, patchGrid);
                    else
                        selectedPatchGrids = (StsPatchGrid[])StsMath.arrayAddElement(selectedPatchGrids, patchGrid);
                }
            }
        }
        StsMessageFiles.logMessage("Picked patch parent id: " + selectedPatch.id + addRemove + nGridsAdded + " child grids.");
    }

    private void clearSelectedPatches()
    {
        selectedPatchGrids = null;
    }

    public StsPatchGrid getNearestPatch(int volumeRow, int volumeCol, float z)
    {
        StsPatchGrid nearestPatch = null;
        float nearestPatchZ = largeFloat;
        //TODO make an iterator which returns patches which cross this volumeRow
        int[] patchRange = this.getPatchRangeForRow(volumeRow);
        int patchMin = patchRange[0];
        int patchMax = patchRange[1];
        for(int n = patchMin; n <= patchMax; n++)
        {
            StsPatchGrid patchGrid = rowSortedPatchGrids[n];
            float patchZ = patchGrid.getPointZ(volumeRow, volumeCol);
            if(patchZ == nullValue) continue;
            float dz = Math.abs(z - patchZ);
            if(dz < nearestPatchZ)
            {
                nearestPatch = patchGrid;
                nearestPatchZ = dz;
            }
        }
        return nearestPatch;
    }
/*
    public void addSelectedPatch(int patchID)
    {
        if(patchID < 0 || patchID >= rowSortedPatchGrids.length)
        {
            StsException.systemError(this, "addSelectedPatch", "patchID out of range: " + patchID);
            return;
        }
        StsPatchGrid selectedPatch = this.rowSortedPatchGrids[patchID];
        selectedPatchGrids = (StsPatchGrid[])StsMath.arrayAddElement(selectedPatchGrids, selectedPatch);
    }
*/
}