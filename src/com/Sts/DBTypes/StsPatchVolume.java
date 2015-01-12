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

    /** gridList contains patches which have been completed. At the end of a row, prevRowGrids contains disconnected grids which
     *  are added to gridList patches.  At the end of all rows, remaining grids are in rowGrids which are then added to gridList. */
    transient ArrayList<StsPatchGrid> gridList;

    /** rowGrids contains new patches and existing patches from previous row connected to the latest row;
     *  at the completion of the row, these become the previousRowGrids. At start of row, its initialized to empty.  If a new grid is created,
     *  it is added to rowGrids.  If an existing grid is connected to a point in the row, it is added to rowGrid and removed from prevRowGrid.
     *  At the end of the row, grids still connected are in rowGrid, and disconnected ones are in prevRowGrids. These disconnected grids are
     *  added to gridList.  prevRowGrids is then set to rowGrids and rowGrids initialized for the next row. */
    transient HashMap<Integer, StsPatchGrid> rowGrids = null;
    transient Iterator<StsPatchGrid> rowGridsIterator;

    /** prevRowGrids are the active patches in the previousRow; when making a connection to a point in the previous row, we look here for a patch. */
    transient HashMap<Integer, StsPatchGrid> prevRowGrids = null;
    transient Iterator<StsPatchGrid> prevRowGridsIterator;

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
    transient int row, col, volRow, volCol;
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

    transient float[] histogramValues;
    transient int nHistogramValues = 10000;

    /** if the row or col correl is < this value, that link and the associated edge is not drawn */
    static public float minLinkCorrel = 0.0f;

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

    static final boolean debug = false;
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

    /** This is central method for constructing the volume of patchGrids.
     *  For each row, we examine correlation with traces in same row & prev col and same col & prev row.
     *
     * @param pickPanel graphics panel with progress bar updated as each row completed
     */
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
        initializeLists();

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
        StsPatchGrid.staticInitialize();

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
        int nCroppedSlices = croppedBoundingBox.sliceMax - croppedSliceMin + 1;
        int nVolSlices = seismicVolume.nSlices;
        float[] traceValues = new float[nCroppedSlices];

        try
        {
            // row & col refer to the row and col in a croppedVolume over which picker is to run
            // volRow & volCol define the actual row and col in the volume (used only for reference)
            for (row = 0, volRow = croppedRowMin; row < nRows; row++, volRow++)
            {
                //statusArea.setProgress(row*40.f/nRows);
                TracePoints[] prevRowTracesPoints = rowTracePoints;
                rowTracePoints = new TracePoints[nCols];
                TracePoints prevRowTracePoints = null;
                incrementLists();
                // get row plane from seismicVolume at this volRow
                FloatBuffer rowFloatBuffer = seismicVolume.getRowPlaneFloatBuffer(volRow, croppedColMin);
                if (rowFloatBuffer == null) return;
                // if(croppedColMin > 0) rowFloatBuffer.position(croppedColMin * nVolSlices);
                for (col = 0, volCol = croppedColMin; col < nCols; col++, volCol++)
                {
                    // StsException.systemDebug(this, "constructPatchVolume", "col loop, col: " + col);
                    rowFloatBuffer.position(volCol * nSlices + croppedSliceMin);
                    rowFloatBuffer.get(traceValues);

                    TracePoints tracePoints = new TracePoints(row, col, traceValues, pickType, nSlices, croppedSliceMin, nVolSlices);
                    rowTracePoints[col] = tracePoints;
                    // prevColTracePoints are tracePoints in prev row & same col
                    TracePoints prevColTracePoints = null;
                    if (prevRowTracesPoints != null)
                        prevColTracePoints = prevRowTracesPoints[col];
                    addTracePatches(tracePoints, prevColTracePoints, prevRowTracePoints);
                    prevRowTracePoints = tracePoints;
                }
                processPrevRowGrids(row);
                if (progressPanel == null) continue;
                if (progressPanel.isCanceled())
                {
                    progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
                    return;
                }
                progressPanel.setValue(row + 1);
            }
            addRemainingGrids();
            finish();
            getPatchVolumeClass().setDisplayCurvature(true);
            progressPanel.finished();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructPatchVolume", e);
        }
    }

    /** creates correlated connections between this trace and traces at prev row & same col and prev col & same row
     *
     * @param trace new trace where correlated patchPoints are to be added.
     * @param prevColTrace prev trace in same col, prev row
     * @param prevRowTrace prev trace in same row, prev col
     */
    private void addTracePatches(TracePoints trace, TracePoints prevColTrace, TracePoints prevRowTrace)
    {
        if(prevRowTrace == null && prevColTrace == null) return;

        StsPatchPoint[] newPatchPoints = trace.tracePatchPoints;
        int nNewPatchPoints = newPatchPoints.length;
        if(nNewPatchPoints == 0) return;

        // For each patchPoint, check previous traces for points which are same type and just above or below this new point
        // and find which of these two possible prev trace points has the best correlation.
        // If this correlation is above the minCorrelation, add this new point to the prev point patch.
        // If the new point already has a patch (because it correlated with one of the other of the 4 otherTraces, then the addPatchPointToPatch
        // will merge the two patches.
        for(int i = 0; i < nIterations; i++)
        {
            initializeTraceIndices(trace, prevRowTrace, prevColTrace);
            float minStretchCorrelation = stretchCorrelations[i];
            for(int centerPointIndex = 0; centerPointIndex < nNewPatchPoints; centerPointIndex++)
            {
                StsPatchPoint centerPoint = newPatchPoints[centerPointIndex];
                if(centerPoint.getPatchGrid() != null) continue;
                CorrelationWindow window;
                window = new CorrelationWindow(trace, centerPoint, centerPointIndex);
                if(!window.initialize()) continue;
                if(prevColTrace != null) checkAddColConnection(window, prevColTrace, minStretchCorrelation);
                if(prevRowTrace != null) checkAddRowConnection(window, prevRowTrace, minStretchCorrelation);
            }
        }
    }

    private void checkAddRowConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
    {
        checkAddConnection(window, otherTrace, minStretchCorrelation, true);
    }

    private void checkAddColConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
    {
        checkAddConnection(window, otherTrace, minStretchCorrelation, false);
    }

    private void checkAddConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation, boolean isRow)
    {
        CorrelationWindow matchingWindow = window.findMatchingWindows(otherTrace, minStretchCorrelation);
        if(matchingWindow == null) return;
        StsPatchPoint centerPoint = window.centerPoint;
        StsPatchPoint otherCenterPoint = matchingWindow.centerPoint;
        otherTrace.centerPointIndex = matchingWindow.centerPointIndex;
        double distance = Math.abs(otherCenterPoint.slice - window.centerSlice);
        float correl = matchingWindow.stretchCorrelation;
        addPatchConnection(centerPoint, otherCenterPoint, correl, isRow);
    }
    /**
     * Given a newPatchPoint at newRow-newCol, which correlates with a prevPatchPoint at prevRow-prevCol which is possibly part of a patchGrid in the prevPatchGridsSet,
     * combine these two points in the same patch.  The prevPatchPoint may be on the previous col (same row), or previous row (same col).
     * If the previousPatchPoint is not part of an existing patchGrid (prevID == -1), then we will create a new patchGrid and add both points to it.
     * If the previousPatchPoint is part of a patchGrid we will add the newPatchPoint to this patchGrid, unless the newPatchPoint already belongs to another patchGrid
     * (this occurs when we first correlate with the previous column and find one patchGrid and also correlate with the previous row and find a different patchGrid).
     */
    public void addPatchConnection(StsPatchPoint newPatchPoint, StsPatchPoint otherPatchPoint, float correl, boolean isRow)
    {
        StsPatchGrid patchGrid = null;

        if(correl < minLinkCorrel) return;
        StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
        StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

        if(newPatchGrid == null)
        {
            if(otherPatchGrid == null) // prevPatchGrid doesn't exist, so create it and add otherPoint to it
            {
                patchGrid = StsPatchGrid.construct(this, newPatchPoint.pointType);
                patchGrid.addPatchPoint(otherPatchPoint);
            }
            else // otherPatchGrid does exist, so use it
            {
                patchGrid = otherPatchGrid;
                if(patchGrid == null) return;
            }
            // this point can't overlap this grid, so we don't need to call checkAddPatchPoint
            patchGrid.addPatchPoint(newPatchPoint);
            //patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
        }
        else // id != -1 means this point was just added to a patch from prev connection and the patchGrid would have been added to the rowGrids array
        {
            if(otherPatchGrid == null) // the otherPoint is not assigned to a patch; assign it to this one; don't add point to rowGrids unless it overlaps and addedGrid created
            {
                patchGrid = newPatchGrid;
                if(patchGrid == null) return; // patchGrid not found; systemDebugError was printed
                // otherPatchPoint doesn't have a patchGrid, but newPatchPoint does; try to add otherPatchPoint to newPatchGrid,
                // but if it overlaps, created an addedPatchGrid containing otherPatchPoint and a clone of newPatchPoint
                // return this addedPatchGrid or patchGrid (if no new grid added)
                patchGrid = patchGrid.checkAddPatchPoint(otherPatchPoint, newPatchPoint);
                //patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
                //checkAddPatchGridToRowGrids(patchGrid);
            }
            else if(otherPatchGrid.id == newPatchGrid.id) // otherPoint is already assigned to the same patch: addCorrelation
            {
                patchGrid = newPatchGrid;
                if(patchGrid == null) return; // patchGrid not found; systemDebugError was printed
                //checkAddPatchGridToRowGrids(patchGrid);
                //patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
            }
            // prevPoint and this point belong to different patches: merge newPatchGrid into prevPatchGrid and add connection
            // if we can't merge OK, then we create a new patch with newPoint and clone of connected otherPoint
            else
                patchGrid = checkMergePatchGrids(otherPatchPoint, newPatchPoint, correl);
        }

        if(patchGrid != null)
        {
            patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
            checkAddPatchGridToRowGrids(patchGrid);
        }
    }

    private StsPatchGrid checkMergePatchGrids(StsPatchPoint otherPatchPoint, StsPatchPoint newPatchPoint, float correl)
    {
        StsPatchGrid mergedGrid, removedGrid;

        StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
        if(otherPatchGrid == null) return null;
        StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
        if(newPatchGrid == null) return null;
        if(StsPatchGrid.mergePatchPointsOK(otherPatchGrid, newPatchGrid))
        {
            if(otherPatchGrid.id < newPatchGrid.id)
            {
                mergedGrid = otherPatchGrid;
                removedGrid = newPatchGrid;
            }
            else
            {
                mergedGrid = newPatchGrid;
                removedGrid = otherPatchGrid;
            }
            // merge shouldn't fail, so a return of false indicates a problem: bail out
            if(!mergedGrid.mergePatchPoints(removedGrid)) return null;

            //checkAddPatchGridToRowGrids(mergedGrid);
            //mergedGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
            removePatchGridFromLists(removedGrid);
            return mergedGrid;
        }
        else // can't merge: find larger grid and add a clone of the overlapped point from the smaller grid to it
        {
            StsPatchGrid changedGrid;
            StsPatchPoint clonedPoint;
            if(otherPatchGrid.nPatchPoints >= newPatchGrid.nPatchPoints)
            {
                changedGrid = otherPatchGrid;
                clonedPoint = newPatchPoint.clone();
                changedGrid.addPatchPoint(clonedPoint);
                changedGrid.addCorrelation(otherPatchPoint, clonedPoint, correl);
            }
            else
            {
                changedGrid = newPatchGrid;
                clonedPoint = otherPatchPoint.clone();
                changedGrid.addPatchPoint(clonedPoint);
                changedGrid.addCorrelation(clonedPoint, newPatchPoint, correl);
            }
            checkAddPatchGridToRowGrids(changedGrid);
            return null; // return null indicating this grid has been completely processed
        }
    }


    private void checkAddPatchGridToRowGrids(StsPatchGrid patchGrid)
    {
        int patchID = patchGrid.id;
        if(patchGrid.rowGridAdded) return;
        boolean debug = patchGrid.debug();
        StsPatchGrid value = rowGrids.put(patchID, patchGrid); // if return is null, no value exists at this key
        patchGrid.rowGridAdded = true;
        if(debug)
        {
            if(value == null)
                StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " added to rowGrids for row: " + row + " col: " + col);
            else
                StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " already exists for row: " + row + " col: " + col);
        }
    }

    private void removePatchGridFromLists(StsPatchGrid patchGrid)
    {
        StsPatchGrid value;
        int patchID = patchGrid.id;
        boolean debug = StsPatchGrid.debugPatchID != -1 && patchID == StsPatchGrid.debugPatchID;
        value = prevRowGrids.remove(patchID);
        if(debug)
        {
            if(value != null)
                StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from prevRowGrids for row: " + row);
            else
                StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in prevRowGrids for row: " + row);
        }
        value = rowGrids.remove(patchID);
        if(debug)
        {
            if(value != null)
                StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from rowGrids for row: " + row);
            else
                StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in rowGrids for row: " + row);
        }
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
        float stretchCorrelation;
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

        /**
         * On the other trace, starting from the point following the last correlated point,
         * find any and all picks of the center pick type on the other trace which match.
         * check that the above and below pick types match as well
         */
        CorrelationWindow newFindMatchingWindows(TracePoints otherTrace, float stretchCorrelation)
        {
            // sliceIndex for centerPoint of this window for which we want to find a matchingWindow on otherTrace
            int centerSlice = centerPoint.slice;
            // sliceIndex for last event found on otherTrace (or zero); we will search up and down from here for matches
            int otherTraceCenterPointIndex = otherTrace.centerPointIndex + 1;

            // search down otherTrace for event just below centerSlice of this window
            while (otherTraceCenterPointIndex < otherTrace.nTracePatchPoints)
            {
                StsPatchPoint otherTraceCenterPoint = otherTrace.tracePatchPoints[otherTraceCenterPointIndex];
                if (otherTraceCenterPoint.slice > centerSlice)
                    return findMatchingWindows(otherTrace, otherTraceCenterPointIndex, stretchCorrelation);
                otherTraceCenterPointIndex++;
            }
            return null;
        }

        /**
         *
         * @param otherTrace the trace on which we want to find matching windows
         * @param otherTraceCenterPointStartIndex index in otherTrace where we will start the search
         * @param stretchCorrelation minimum value used in matching windows @see matches()
         * @return best correlation window or null if none
         */
        CorrelationWindow findMatchingWindows(TracePoints otherTrace, int otherTraceCenterPointStartIndex, float stretchCorrelation)
        {
            // this is sliceIndex for centerPoint of this window for which we want to find a matchingWindow on otherTrace
            int centerSlice = centerPoint.slice;
            // otherTrace point closest to this centerPoint where we start our search for matches
            StsPatchPoint otherTraceCenterPoint = otherTrace.tracePatchPoints[otherTraceCenterPointStartIndex];
            // this is sliceIndex for last event found on otherTrace (or zero); we will search up and down from here for matches
            int otherTraceCenterPointLastIndex = otherTrace.centerPointIndex + 1;
            // search down otherTrace for event just below centerSlice of this window
            int nOtherTracePoints = otherTrace.nTracePatchPoints;

            CorrelationWindow matchingWindow = null; // this will be best matching window in otherTrace
            double bestCorrelation = 0.0; // this will be correlation between best otherTraceWindow and traceWindow
            CorrelationWindow otherWindow; // possible matchingWindow to check for best match

            int otherCenterSlice = otherTraceCenterPoint.slice;
            // otherCenterSlice is just below centerSlice; if the pointType matches, try a correlation
            if (otherTraceCenterPoint.pointType == centerPointType)
            {
                otherWindow = checkCreateMatchingWindow(otherTrace, otherTraceCenterPoint, otherCenterSlice, otherTraceCenterPointStartIndex, centerSlice, stretchCorrelation);
                if (otherWindow != null && otherWindow.stretchCorrelation > bestCorrelation)
                {
                    matchingWindow = otherWindow;
                    bestCorrelation = matchingWindow.stretchCorrelation;
                }

            }
            // if the point just below centerSlice doesn't have same patchType, check the point just above (as long as it is still below the last correlated event index
            else if (otherTraceCenterPointStartIndex > otherTraceCenterPointLastIndex)
            {
                StsPatchPoint nextOtherTraceCenterPoint = otherTrace.tracePatchPoints[otherTraceCenterPointStartIndex - 1];
                if (nextOtherTraceCenterPoint.pointType == centerPointType)
                {
                    otherTraceCenterPointStartIndex--;
                    otherTraceCenterPoint = nextOtherTraceCenterPoint;
                    // we have a pointType match just below, so try to match it
                    otherWindow = checkCreateMatchingWindow(otherTrace, otherTraceCenterPoint, otherCenterSlice, otherTraceCenterPointStartIndex, centerSlice, stretchCorrelation);
                    if (otherWindow != null && otherWindow.stretchCorrelation > bestCorrelation)
                    {
                        matchingWindow = otherWindow;
                        bestCorrelation = matchingWindow.stretchCorrelation;
                    }
                }
            }
            int nextOtherTraceCenterPointIndex;
            // now search up for a match, but don't go past an event which is already correlated
            for (nextOtherTraceCenterPointIndex = otherTraceCenterPointStartIndex - 1; nextOtherTraceCenterPointIndex > 0; nextOtherTraceCenterPointIndex--)
            {
                otherTraceCenterPoint = otherTrace.tracePatchPoints[nextOtherTraceCenterPointIndex];
                otherCenterSlice = otherTraceCenterPoint.slice;
                if (centerPointType == otherTraceCenterPoint.pointType)
                {
                    otherWindow = checkCreateMatchingWindow(otherTrace, otherTraceCenterPoint, otherCenterSlice, nextOtherTraceCenterPointIndex, centerSlice, stretchCorrelation);
                    if (otherWindow != null && otherWindow.stretchCorrelation > bestCorrelation)
                    {
                        matchingWindow = otherWindow;
                        bestCorrelation = matchingWindow.stretchCorrelation;
                    }
                    break;
                }
            }
            // now search down for a match, but don't go past an event which is already correlated
            for (nextOtherTraceCenterPointIndex = otherTraceCenterPointStartIndex + 1; nextOtherTraceCenterPointIndex < nOtherTracePoints; nextOtherTraceCenterPointIndex++)
            {
                otherTraceCenterPoint = otherTrace.tracePatchPoints[nextOtherTraceCenterPointIndex];
                if (otherTraceCenterPoint.getPatchGrid() != null) break;
                otherCenterSlice = otherTraceCenterPoint.slice;
                if (centerPointType == otherTraceCenterPoint.pointType)
                {
                    otherWindow = checkCreateMatchingWindow(otherTrace, otherTraceCenterPoint, otherCenterSlice, nextOtherTraceCenterPointIndex, centerSlice, stretchCorrelation);
                    if (otherWindow != null && otherWindow.stretchCorrelation > bestCorrelation)
                    {
                        matchingWindow = otherWindow;
                        bestCorrelation = matchingWindow.stretchCorrelation;
                    }
                    break;
                }
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
            float minusStretchFactor = 1.0f;
            float plusStretchFactor = 1.0f;
            if(windowType != WINDOW_BELOW)
            {
                float dzMinusOtherScalar = ((float)dSliceMinus) / otherWindow.dSliceMinus;
                minusStretchFactor = dzMinusOtherScalar;
                if(minusStretchFactor > 1.0f)
                    minusStretchFactor = 1 / minusStretchFactor;
            }
            else if(windowType != WINDOW_ABOVE)
            {
                float dzPlusOtherScalar = ((float)dSlicePlus) / otherWindow.dSlicePlus;
                plusStretchFactor = dzPlusOtherScalar;
                if(plusStretchFactor > 1.0f)
                    plusStretchFactor = 1 / plusStretchFactor;
            }
            float stretchFactor = Math.min(minusStretchFactor, plusStretchFactor);
            otherWindow.stretchCorrelation = stretchFactor;
            return stretchFactor >= stretchCorrelation;
        }

        boolean isZCenterOutsideWindow(int centerSlice)
        {
            return centerSlice < minSlice || centerSlice > maxSlice;
        }
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

    private void initializeLists()
    {
        gridList = new ArrayList<>(100);
        rowGrids = new HashMap<>();
        rowGridsIterator = rowGrids.values().iterator();
        prevRowGrids = new HashMap<>(); // not used on first row
    }

    private void incrementLists()
    {
        prevRowGrids = rowGrids;
        prevRowGridsIterator = rowGridsIterator;
        clearPrevRowGridsAddedFlags();
        rowGrids = new HashMap<>();
        rowGridsIterator = rowGrids.values().iterator();
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

    private void initializeTraceIndices(TracePoints trace, TracePoints prevRowTrace, TracePoints prevColTrace)
    {
        if(trace != null) trace.centerPointIndex = -1;
        if(prevRowTrace != null) prevRowTrace.centerPointIndex = -1;
        if(prevColTrace != null) prevColTrace.centerPointIndex = -1;
    }

    public StsPatchVolumeClass getPatchVolumeClass()
    {
        if(patchVolumeClass != null) return patchVolumeClass;
        patchVolumeClass = (StsPatchVolumeClass)getCreateStsClass();
        return patchVolumeClass;
    }

    private void finish()
    {
        // build patchGrid arrays. Overlapped grids will be split out into new grids
        // So construct an array of existing patchGrids and construct the various 2D arrays for each.
        // New grids generated will be added to the gridList
        int nGrids = gridList.size();
        rowSortedPatchGrids = new StsPatchGrid[nGrids];
        gridList.toArray(rowSortedPatchGrids);
        // debugCheckEmptyFraction(rowSortedPatchGrids);
        StsException.systemDebug(this, "finish", " Number of parent grids: " + nParentGrids + "Number of child grids: " + nGrids + " too small: " + nSmallGridsRemoved);
        // StsException.systemDebug(this, "finish", "max grid dimension: " + StsPatchGrid.maxGridSize);

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
        if(!isPersistent())
        {
            currentModel.getDatabase().blockTransactions();
            addToModel();
            currentModel.getDatabase().saveTransactions();
        }
        getPatchVolumeClass().setIsVisible(true);
        setIsVisible(true);
    }
/*
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
*/
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

    private void clearPrevRowGridsAddedFlags()
    {
        Iterator<StsPatchGrid> prevRowGridIterator = prevRowGrids.values().iterator();
        while(prevRowGridIterator.hasNext())
        {
            StsPatchGrid patchGrid = prevRowGridIterator.next();
            patchGrid.rowGridAdded = false;
        }
    }

    /** Called for the last row only; unconditionally add all rowGrids unless they are too small. */
    void addRemainingGrids()
    {
        StsPatchGrid[] patchGrids = rowGrids.values().toArray(new StsPatchGrid[0]);
        for(int  n = 0; n < patchGrids.length; n++)
        {
            StsPatchGrid patchGrid = patchGrids[n];
            if(!patchGrid.isTooSmall(nPatchPointsMin))
            {
                patchGrid.finish();
                gridList.add(patchGrid);
            }
        }
        if(StsPatchVolume.debug)
            StsException.systemDebug(this, "addRemainingGrids", "added " + patchGrids.length + " grids remaining.");
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
            if(patch != null) patch.clear();
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

        return Math.min(minusStretchFactor, plusStretchFactor);
    }

    private StsPatchGrid getPatchGrid(int id)
    {
        StsPatchGrid patchGrid = rowGrids.get(id);
        // if patchGrid exists in rowGrids, then it has already been added there and deleted from prevRowGrids
        if(patchGrid != null)
        {
            if(StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
                StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
                        " gotten from rowGrids at row: " + row + " col: " + col);
            return patchGrid;
        }
        // if patchGrid is not in rowGrids, then add it there and delete it from prevRowGrids
        else if(prevRowGrids != null)
        {
            patchGrid = prevRowGrids.get(id);
            if(patchGrid != null)
            {
                rowGrids.put(id, patchGrid);
                prevRowGrids.remove(id);
                if(StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
                    StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
                            " gotten and deleted from prevRowsGrids and added to rowGrids at row: " + row + " col: " + col);
                return patchGrid;
            }
        }
        StsException.systemError(this, "getPatchGrid",  "Couldn't get patchGrid for id " + id + " at row: " + " col: " + col);
        return null;
    }

    /**
     * This PatchGridSet is for the row before row just finished.
     * If a grid in this prev row is disconnected (doesn't have same patch in row just finished),
     * then either delete it if it is a small point, or add it to volume set.
     */
    void processPrevRowGrids(int row)
    {
        if(row == 0) return;
        StsPatchGrid[] prevRowPatchGrids = prevRowGrids.values().toArray(new StsPatchGrid[0]);
        int nDisconnectedGrids = 0;
        for(StsPatchGrid patchGrid : prevRowPatchGrids)
        {
            boolean disconnected = patchGrid.isDisconnected(row);
            if(!disconnected) continue;
            if(patchGrid.isTooSmall(nPatchPointsMin))
                nSmallGridsRemoved++;
            else
            {
                patchGrid.finish();
                gridList.add(patchGrid);
                nDisconnectedGrids++;
            }
            prevRowGrids.remove(patchGrid.id);
        }
        if(debug) StsException.systemDebug(this, "processPrevRowGrids", "prev row: " + (row - 1) + " added " + nDisconnectedGrids + " disconnected grids");
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
                    // displayPatchesNearXYZCursors(glPanel3d);
                    displayPatchesNearZCursor(glPanel3d, dirCoordinate);
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