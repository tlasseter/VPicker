package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SurfaceCurvature.*;
import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
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
	protected boolean filter = false;
	protected int boxFilterWidth = 1;
	transient public int nInterpolatedSlices;
	transient public float interpolatedZInc;
	transient public int interpolatedSliceMin;
	transient public int interpolatedSliceMax;
	transient StsPatchVolumeClass patchVolumeClass;
	transient public StsCroppedBoundingBox croppedBoundingBox;

	/**
	 * gridList contains patches which have been completed. At the end of a row, prevRowGrids contains disconnected grids which
	 * are added to gridList patches.  At the end of all rows, remaining grids are in rowGrids which are then added to gridList.
	 */
	transient ArrayList<StsPatchGrid> gridList;

	/**
	 * rowGrids contains new patches and existing patches from previous row connected to the latest row;
	 * at the completion of the row, these become the previousRowGrids. At start of row, its initialized to empty.  If a new grid is created,
	 * it is added to rowGrids.  If an existing grid is connected to a window in the row, it is added to rowGrid and removed from prevRowGrid.
	 * At the end of the row, grids still connected are in rowGrid, and disconnected ones are in prevRowGrids. These disconnected grids are
	 * added to gridList.  prevRowGrids is then set to rowGrids and rowGrids initialized for the next row.
	 */
	transient HashMap<Integer, StsPatchGrid> rowGrids = null;
	transient Iterator<StsPatchGrid> rowGridsIterator;

	/** prevRowGrids are the active patches in the previousRow; when making a connection to a window in the previous row, we look here for a patch. */
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
	/** pick window on adjoining trace cannot be more than this many wavelengths away from picked window */
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
	/** increment of stretch stretchCorrelation in interative pick */
	transient float autoCorInc;
	/** manual picking minimum acceptable cross-stretchCorrelation */
	transient float manualCorMin;
	/** minimum amplitude fraction of sample data max allowed for an event to be correlated */
	transient float minAmpFraction;
	/** sequence of stretchCorrelations: 1 if not iterative, max to min by inc if iterative */
	transient float[] stretchCorrelations;
	/** number of stretchCorrelations in sequence */
	transient int nIterations;
	/** min amplitudeRatio allowed for correlation */
	transient float minAmplitudeRatio;
	/** correlate using falseTypes (e.g., a false Max matches a Max, etc) */
	transient boolean useFalseTypes = false;
	/** double check connection by matching it back from selected prevWindow; accept match if backMatchWindow is null or the same */
	transient boolean checkBackMatch = true;
	/** row currently being computed: used for debugPatchGrid print out only */
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

	public static final int largeInt = 99999999;

	public static final String[] pickTypeNames = new String[]{"All", "Min+Max", "Maximum", "Minimum"}; //, "Zero-crossing+", "Zero-crossing-", "All"};
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
		return zMin + this.interpolatedZInc * slice;
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

	/** debugPatchGrid prints showing row operations */
	static public final boolean debug = true;
	/** turn on timer for values operation */
	static final boolean runTimer = false;
	/** millisecond timer */
	static StsTimer timer;
	/** debugPatchGrid for tracePoints linkList operations */
	static final boolean debugTracePointsLink = false;
	/** debugPatchGrid for CorrelationWindows */
	static final boolean debugCorrelationWindows = false;
	/** search for multiple windows */
	static final boolean searchForMultipleWindowMatches = true;
	/** print patch operations and draw only this patch */
	static boolean drawPatchBold = debug && StsPatchGrid.debugPatchGrid; // StsPatchGrid.debugPatchGrid;
	/** debug: allow point clone operations */
	static final boolean debugCloneOK = true;
	/** debug: add overlapping grids to group */
	static final boolean debugOverlapGridGroupOK = false;
	/** debug: connect closest points only */
	static final boolean debugConnectCloseOnly = true;
	static public String iterLabel = "";

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

	/**
	 * This is central method for constructing the volume of patchGrids.
	 * For each row, we examine stretchCorrelation with traces in same row & prev col and same col & prev row.
	 * @param pickPanel graphics panel with progress bar updated as each row completed
	 */
	public void constructPatchVolume(StsPatchPickPanel pickPanel)
	{
		StsProgressBar progressPanel = pickPanel.progressBar;
		windowSize = pickPanel.corWavelength;
		pickDifWavelengths = pickPanel.maxPickDif;
		minAmpFraction = pickPanel.minAmpFraction;
		float maxStretch = pickPanel.maxStretch;
		pickType = pickPanel.pickType;
		nPatchPointsMin = pickPanel.minPatchSize;
		useFalseTypes = pickPanel.useFalseTypes;
		checkBackMatch = pickPanel.checkBackMatch;
		isIterative = pickPanel.isIterative;
		autoCorMax = pickPanel.autoCorMax;
		autoCorMin = pickPanel.autoCorMin;
		autoCorInc = pickPanel.autoCorInc;
		manualCorMin = pickPanel.manualCorMin;
		minAmplitudeRatio = pickPanel.minAmpRatio;

		StsPatchGrid.initializeDebug(pickPanel);
		initializeLists();

		if (!isIterative)
		{
			nIterations = 1;
			stretchCorrelations = new float[]{manualCorMin};
		}
		else
		{
			nIterations = StsMath.ceiling(1 + (autoCorMax - autoCorMin) / autoCorInc);
			stretchCorrelations = new float[nIterations];
			float stretchCorrelation = autoCorMax;
			for (int n = 0; n < nIterations; n++, stretchCorrelation -= autoCorInc)
				stretchCorrelations[n] = stretchCorrelation;
		}

		rowSortedPatchGrids = new StsPatchGrid[0];
		colSortedPatchGrids = null;
		initialize();
		StsPatchGrid.staticInitialize();

		if (progressPanel != null)
			progressPanel.initialize(croppedBoundingBox.nRows);

		TracePoints[] rowTracePoints = null; //new TracePoints[nCols];
		//hack:  FIX!
		if (seismicVolume == null)
			seismicVolume = (StsSeismicVolume) currentModel.getCurrentObject(StsSeismicVolume.class);
		float absSeismicDataMax = Math.min(Math.abs(seismicVolume.dataMin), Math.abs(seismicVolume.dataMax));
		minDataAmplitude = minAmpFraction * absSeismicDataMax;
		initializeParameters();
		initializeToBoundingBox(croppedBoundingBox);
		initializeSliceInterpolation();

		int croppedRowMin = croppedBoundingBox.rowMin;
		int croppedRowMax = croppedBoundingBox.rowMax;
		int croppedColMin = croppedBoundingBox.colMin;
		int croppedColMax = croppedBoundingBox.colMax;
		int croppedSliceMin = croppedBoundingBox.sliceMin;
		int croppedSliceMax = croppedBoundingBox.sliceMax;
		int nCroppedSlices = croppedSliceMax - croppedSliceMin + 1;
		int nVolSlices = seismicVolume.nSlices;
		float[] traceValues = new float[nCroppedSlices];

		float[] nullTrace = new float[nCroppedSlices];
		//TODO fix seismic 3d io processing to save proper userNull
		Arrays.fill(nullTrace, seismicVolume.userNull);

		try
		{
			// row & col refer to the row and col in a croppedVolume over which picker is to run
			// volRow & volCol define the actual row and col in the volume (used only for reference)
			for (row = 0, volRow = croppedRowMin; volRow <= croppedRowMax; row++, volRow++)
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
				for (col = 0, volCol = croppedColMin; volCol <= croppedColMax; col++, volCol++)
				{
					// StsException.systemDebug(this, "constructPatchVolume", "col loop, col: " + col);
					rowFloatBuffer.position(volCol * nVolSlices + croppedSliceMin);
					rowFloatBuffer.get(traceValues);

					TracePoints tracePoints = null;
					if(!Arrays.equals(traceValues, nullTrace))
						tracePoints = TracePoints.constructor(this, row, col, traceValues);
					rowTracePoints[col] = tracePoints;
					if(tracePoints == null) continue;
						// prevColTracePoints are tracePoints in prev row & same col
					TracePoints prevColTracePoints = null;
					if (prevRowTracesPoints != null)
						prevColTracePoints = prevRowTracesPoints[col];

					// here we add the connected patchPoints
					tracePoints.connectWindows(prevColTracePoints, prevRowTracePoints);

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
			getPatchVolumeClass().setDisplayCurvature(false);
			progressPanel.finished();
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "constructPatchVolume", e);
		}
	}


	static void setIterLabel(int iter)
	{
	 	iterLabel = " Iter: " + Integer.toString(iter);
	}

	public boolean setFilter()
	{
		return filter;
	}

	public void setFilter(boolean filter)
	{
		this.filter = filter;
	}

	public int getBoxFilterWidth()
	{
		return boxFilterWidth;
	}

	public void setBoxFilterWidth(int boxFilterWidth)
	{
		this.boxFilterWidth = boxFilterWidth;
	}

	protected StsPatchGrid mergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint)
	{
		StsPatchGrid mergedGrid, removedGrid;

		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

		if (otherPatchGrid.id < newPatchGrid.id)
		{
			mergedGrid = otherPatchGrid;
			removedGrid = newPatchGrid;
		}
		else
		{
			mergedGrid = newPatchGrid;
			removedGrid = otherPatchGrid;
		}
		if (StsPatchGrid.debugPatchID != StsPatchGrid.NO_DEBUG && (mergedGrid.id == StsPatchGrid.debugPatchID || removedGrid.id == StsPatchGrid.debugPatchID))
			StsException.systemDebug(this, "mergePatchGrids", "MERGING GRID " + removedGrid.toGridString() + " TO GRID " + mergedGrid.toGridString());
		// merge shouldn't fail, so a return of false indicates a problem: bail out
		if (!mergedGrid.mergePatchPoints(removedGrid))
		{
			StsException.systemError(this, "mergePatchGrids", "Failed to merge removedGrid " + removedGrid.toGridString() + " to mergedGrid " + mergedGrid.toGridString());
			return null;
		}
		// removedGrid could be parent, child, or new; it is being merged into mergedGrid which could be any of these three as well
		// if removedGrid is a parent, then the first child needs to be made the parent
		// if removed grid is a child, it needs to be removed from the parent
		// if new, we don't have to do anything

		removedGrid.removeChildOrParent(mergedGrid);

		//checkAddPatchGridToRowGrids(mergedGrid);
		//mergedGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
		removePatchGridFromLists(removedGrid);
		if (StsPatchGrid.debugPatchGrid && removedGrid.id == StsPatchGrid.debugPatchID)
			StsPatchGrid.debugPatchID = mergedGrid.id;

		return mergedGrid;
	}
	protected StsPatchGrid checkAddOverlappingConnection(PatchPoint prevPatchPoint, PatchPoint patchPoint, Connection connection)
	{
		StsPatchGrid prevPatchGrid = prevPatchPoint.getPatchGrid();
		StsPatchGrid patchGrid = patchPoint.getPatchGrid();

		if(StsPatchVolume.debug && (StsPatchGrid.doDebugPoint(patchPoint) || StsPatchGrid.doDebugPoint(prevPatchPoint)))
			StsException.systemDebug(this, "mergeOverlappingPatchGrids", StsPatchVolume.iterLabel + " CONNECT point: " +
					patchPoint.toString() + " TO point: " + prevPatchPoint.toString());

		// if points for this connection mutually overlap the other grid, ignore this connection
		// alternatively, we could make a new patchGrid containing clones of these two connection points
		if (prevPatchGrid.patchPointOverlaps(patchPoint) && patchGrid.patchPointOverlaps(prevPatchPoint))
			return null;
		else if (prevPatchGrid.patchPointOverlaps(patchPoint))
			return patchGrid;
		else if(patchGrid.patchPointOverlaps(prevPatchPoint))
			return prevPatchGrid;
		else // connection itself doesn't overlap, so add connection to largest patch
			return StsPatchGrid.getLargestGrid(patchGrid, prevPatchGrid);
	}

	/** We wish to merge two patchGrids which are connected between these two points.  The two grids overlap.
	 *  If one of the points in the connection doesn't overlap the other grid, then add this point to the other grid.
	 *  Then check all the remaining points in the other grid, and add them to the grid if they don't overlap.
	 *  The remaining points in the other grid overlap the grid, so leave them there and add clones of any of the moved points.
	 *  If both points in the new connection overlap the other grid (mutually overlap), then add the non-overlapping points
	 *  from the small grid to the larger; for the remaining points, add clones of any of the moved points.
 	 * @param prevPatchPoint one of the two points in the new connection between different grids
	 * @param patchPoint the other point in the new connection
	 * @param connection
	 * @return
	 */
	protected StsPatchGrid mergeOverlappingPatchGrids(PatchPoint prevPatchPoint, PatchPoint patchPoint, Connection connection)
	{
		StsPatchGrid prevPatchGrid = prevPatchPoint.getPatchGrid();
		StsPatchGrid patchGrid = patchPoint.getPatchGrid();
		StsPatchGrid changedGrid = null;

		if(StsPatchVolume.debug && (StsPatchGrid.doDebugPoint(patchPoint) || StsPatchGrid.doDebugPoint(prevPatchPoint)))
			StsException.systemDebug(this, "mergeOverlappingPatchGrids", StsPatchVolume.iterLabel + " CONNECT point: " +
					patchPoint.toString() + " TO point: " + prevPatchPoint.toString());

		// if points for this connection mutually overlap the other grid, ignore this connection
		// alternatively, we could make a new patchGrid containing clones of these two connection points
		if (prevPatchGrid.patchPointOverlaps(patchPoint) && patchGrid.patchPointOverlaps(prevPatchPoint)) return null;
		// move overlapping points from smaller prevPatchGrid to larger patchGrid
		if(patchGrid.nPatchPoints >= prevPatchGrid.nPatchPoints)
		{
			patchGrid = patchGrid.moveNonOverlappingPointsFrom(prevPatchGrid, connection, true);
			// if(StsPatchVolume.debug) StsPatchGrid.debugCheckOverlappedGrids(prevPatchGrid, patchGrid);
			return patchGrid;
		}
		else  // move overlapping points from smaller patchGrid to larger prevPatchGrid
		{
			patchGrid = prevPatchGrid.moveNonOverlappingPointsFrom(patchGrid, connection, false);
			// if(StsPatchVolume.debug) StsPatchGrid.debugCheckOverlappedGrids(patchGrid, prevPatchGrid);
			return patchGrid;
		}
	}
 /*
	protected StsPatchGrid XXmergeOverlappingPatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, Connection connection)
	{
		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
		// Check if the otherPatchPoint doesn't overlap the newPatchGrid; if not:
		// move it to the newPatchGrid and leave a clone on the otherPatchGrid;
		// connection will move with the newPatchPoint
		if (!newPatchGrid.patchPointOverlaps(otherPatchPoint))
		{
			PatchPoint clonedPoint = otherPatchPoint.clone();
			newPatchGrid.addPatchPoint(otherPatchPoint);
			otherPatchGrid.setPoint(clonedPoint);
			if(StsPatchVolume.debugOverlapGridGroupOK) otherPatchGrid.combineChildGrids(otherPatchGrid);
		}
		// Check if the newPatchPoint doesn't overlap the otherPatchGrid; if not:
		// move it to the otherPatchGrid and leave a clone on the newPatchGrid;
		// remove this connection from the clonedPoint
		else if (!otherPatchGrid.patchPointOverlaps(newPatchPoint))
		{
			PatchPoint clonedPoint = newPatchPoint.clone();
			otherPatchGrid.addPatchPoint(newPatchPoint);
			newPatchGrid.setPoint(clonedPoint);
			clonedPoint.deleteConnection(connection);
			if (StsPatchVolume.debugOverlapGridGroupOK) otherPatchGrid.combineChildGrids(newPatchGrid);
		}
		else // each window overlaps the other grid, so we create a newGrid with clones of both
		{
			StsPatchGrid patchGrid = StsPatchGrid.construct(this, otherPatchGrid.patchType);
			PatchPoint newClonedPoint = newPatchPoint.cloneAndClear();
			patchGrid.addPatchPoint(newClonedPoint);
			PatchPoint otherClonedPoint = otherPatchPoint.cloneAndClear();
			patchGrid.addPatchPoint(otherClonedPoint);
			connection.resetConnection(newClonedPoint, otherClonedPoint);

			if(StsPatchVolume.debugOverlapGridGroupOK) patchGrid.combineChildGrids(newPatchGrid);
			return patchGrid;
		}
		StsPatchGrid patchGrid = null, smallPatchGrid = null;
		PatchPoint patchPoint, smallPatchPoint;

		if(StsPatchVolume.debug && (StsPatchGrid.doDebugPoint(newPatchPoint) || StsPatchGrid.doDebugPoint(otherPatchPoint)))
			StsException.systemDebug(this, "addPatchConnection", StsPatchVolume.iterLabel + " CONNECT point: " +
					newPatchPoint.toString() + " TO point: " + otherPatchPoint.toString());

		// move overlapping points from smaller grid to larger
		if(newPatchGrid.nPatchPoints >= otherPatchGrid.nPatchPoints)
		{
			patchGrid = newPatchGrid;
			smallPatchGrid = otherPatchGrid;
		}
		else
		{
			patchGrid = otherPatchGrid;
			smallPatchGrid = newPatchGrid;
		}
		patchGrid.moveNonOverlappingPointsFrom(smallPatchGrid, connection, false);
		if(StsPatchVolume.debug) StsPatchGrid.debugCheckOverlappedGrids(smallPatchGrid, patchGrid);
		return patchGrid;
	}

	protected StsPatchGrid XmergeOverlappingPatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, Connection connection)
	{
		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
		StsPatchGrid patchGrid = null, smallPatchGrid = null;
		PatchPoint patchPoint, smallPatchPoint;

		if(StsPatchVolume.debug && (StsPatchGrid.doDebugPoint(newPatchPoint) || StsPatchGrid.doDebugPoint(otherPatchPoint)))
			StsException.systemDebug(this, "addPatchConnection", StsPatchVolume.iterLabel + " CONNECT point: " +
					newPatchPoint.toString() + " TO point: " + otherPatchPoint.toString());

		// move overlapping points from smaller grid to larger
		boolean orderReversed = false;
		if(newPatchGrid.nPatchPoints >= otherPatchGrid.nPatchPoints)
		{
			patchGrid = newPatchGrid;
			smallPatchGrid = otherPatchGrid;
			patchPoint = newPatchPoint;
			smallPatchPoint = otherPatchPoint;
		}
		else
		{
			patchGrid = otherPatchGrid;
			smallPatchGrid = newPatchGrid;
			patchPoint = otherPatchPoint;
			smallPatchPoint = newPatchPoint;
			orderReversed = true;
		}
		patchGrid.moveNonOverlappingPointsFrom(smallPatchGrid, connection, false);
		if(StsPatchVolume.debug) StsPatchGrid.debugCheckOverlappedGrids(smallPatchGrid, patchGrid);
		// connection points may now be on same grid (one may have been moved or cloned)
		if(newPatchPoint.patchGrid == otherPatchPoint.patchGrid)
			return newPatchGrid;
		// see if the patchPoint on the largerGrid doesn't overlap a point on the smaller grid; if not, clone it and add to smaller grid
		if (!smallPatchGrid.patchPointOverlaps(patchPoint))
		{
			PatchPoint clonedPoint = patchPoint.clone();
			smallPatchGrid.addPatchPoint(patchPoint);
			patchGrid.setPoint(clonedPoint);
			connection.resetConnection(patchPoint, orderReversed);
			if(StsPatchVolume.debugOverlapGridGroupOK) patchGrid.combineChildGrids(otherPatchGrid);
			return smallPatchGrid;
		}
		else if (!patchGrid.patchPointOverlaps(smallPatchPoint))
		{
			PatchPoint clonedPoint = smallPatchPoint.clone();
			patchGrid.addPatchPoint(smallPatchPoint);
			smallPatchGrid.setPoint(clonedPoint);
			connection.resetConnection(smallPatchPoint, orderReversed);
			if (StsPatchVolume.debugOverlapGridGroupOK) patchGrid.combineChildGrids(newPatchGrid);
			return patchGrid;
		}
		else // each window overlaps the other grid, so we create a newGrid with clones of both
		{
			patchGrid = StsPatchGrid.construct(this, otherPatchGrid.patchType);
			PatchPoint newClonedPoint = newPatchPoint.cloneAndClear();
			patchGrid.addPatchPoint(newClonedPoint);
			PatchPoint otherClonedPoint = otherPatchPoint.cloneAndClear();
			patchGrid.addPatchPoint(otherClonedPoint);
			connection.resetConnection(newClonedPoint, false);

			if(StsPatchVolume.debugOverlapGridGroupOK) patchGrid.combineChildGrids(newPatchGrid);
			return patchGrid;
		}
	}

	StsPatchGrid cantMergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
		// can't merge: find larger grid and add a clone of the overlapped window from the smaller grid to it

		StsPatchGrid changedGrid;
		PatchPoint clonedPoint;
		if (otherPatchGrid.nPatchPoints >= newPatchGrid.nPatchPoints)
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
		return changedGrid; // return null indicating this grid has been completely processed
	}
*/

	protected void checkAddPatchGridToRowGrids(StsPatchGrid patchGrid)
	{
		if(patchGrid == null) return;
		int patchID = patchGrid.id;
		if (patchGrid.rowGridAdded) return;
		StsPatchGrid value = rowGrids.put(patchID, patchGrid); // if return is null, no value exists at this key
		patchGrid.rowGridAdded = true;
		if (debug && StsPatchGrid.debugPatchID == patchGrid.id)
		{
			if (value == null)
				StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " added to rowGrids for row: " + row + " col: " + col);
			else
				StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " already exists for row: " + row + " col: " + col);
		}
	}

	protected void removePatchGridFromLists(StsPatchGrid patchGrid)
	{
		StsPatchGrid value;
		int patchID = patchGrid.id;
		boolean debug = StsPatchGrid.debugPatchGrid && patchID == StsPatchGrid.debugPatchID;
		value = prevRowGrids.remove(patchID);
		if (debug)
		{
			if (value != null)
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from prevRowGrids for row: " + row);
			else
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in prevRowGrids for row: " + row);
		}
		value = rowGrids.remove(patchID);
		if (debug)
		{
			if (value != null)
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from rowGrids for row: " + row);
			else
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in rowGrids for row: " + row);
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
			if (colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				colorscale = new StsColorscale("Curvature", spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW), dataMin, dataMax);
				colorscale.setEditRange(dataMin, dataMax);
				colorscale.addActionListener(this);
			}
			colorscale.setRange(dataMin, dataMax);
		}
		catch (Exception e)
		{
			StsException.outputException("StsPatchcVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

	public void actionPerformed(ActionEvent e)
	{
		if (e.getSource() instanceof StsColorscale)
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

		if (windowSize == 0) // windowSize is 0.5 wavelengths, half-window is 0.25 wavelengths
		{
			nHalfSamples = 1;
			windowEndIsZeroCrossing = true;
			halfWindowSize = 0.25f;
			halfWindowPickDifFactor = pickDifWavelengths / 0.25f;
		}
		else
		{
			boolean isHalfWave = !StsMath.isEven(windowSize);
			if (isHalfWave)
			{
				// window size is odd, so window ends with zero-crossing; half-window size is windowSize/2.
				// we need to find (windowSize +1)/2 zero-crossings above and below window center (which is a max or min).
				nHalfSamples = (windowSize + 1) / 2;
				windowEndIsZeroCrossing = true;
				halfWindowPickDifFactor = pickDifWavelengths * 2 / windowSize;
			}
			else
			{
				// window size is even, so window ends with same window type as center; half-window size is windowSize/2.
				// we need to find windowSize/2 points above and below with same window type as window center (which is a max or min).
				nHalfSamples = windowSize / 2;
				halfWindowPickDifFactor = pickDifWavelengths / nHalfSamples;
				windowEndIsZeroCrossing = false;
			}
		}
	}

	public StsPatchVolumeClass getPatchVolumeClass()
	{
		if (patchVolumeClass != null) return patchVolumeClass;
		patchVolumeClass = (StsPatchVolumeClass) getCreateStsClass();
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
		for (StsPatchGrid grid : rowSortedPatchGrids)
			grid.resetIndex(nGrids++);

		// print out grid group summary
		int nParents = 0, nNew = 0, nChildren = 0;
		for (StsPatchGrid grid : rowSortedPatchGrids)
		{
			if (grid.isParent()) nParents++;
			else if (grid.isChild()) nChildren++;
			else nNew++;
		}
		StsException.systemDebug(this, "finish", "grid groups: " + nParents + " child grids: " + nChildren + " unattached grids: " + nNew);

		int num = getPatchVolumeClass().getSize();
		String newname = seismicName + ".patchVolume" + (num + 1);
		setName(newname);
		clearConstructionArrays();
		if (!isPersistent())
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
		if (rowSortedPatchGrids == null) return;
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			nPointsTotal += patchGrid.nPatchPoints;
	}

	private boolean checkColSortedPatchGrids()
	{
		if (colSortedPatchGrids != null) return true;
		if (rowSortedPatchGrids == null) return false;
		int nGrids = rowSortedPatchGrids.length;
		if (nGrids == 0) return false;
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
		while (high - low > 1)
		{
			probe = (high + low) / 2;
			int id = gridList.get(probe).originalID;
			if (id > target)
				high = probe;
			else
				low = probe;
		}
		if (low == -1 || gridList.get(low).originalID != target)
			return null;
		else
			return gridList.get(low);
	}

	private void clearPrevRowGridsAddedFlags()
	{
		Iterator<StsPatchGrid> prevRowGridIterator = prevRowGrids.values().iterator();
		while (prevRowGridIterator.hasNext())
		{
			StsPatchGrid patchGrid = prevRowGridIterator.next();
			patchGrid.rowGridAdded = false;
		}
	}

	/** Called for the last row only; unconditionally add all rowGrids unless they are too small. */
	void addRemainingGrids()
	{
		StsPatchGrid[] patchGrids = rowGrids.values().toArray(new StsPatchGrid[0]);
		for (int n = 0; n < patchGrids.length; n++)
		{
			StsPatchGrid patchGrid = patchGrids[n];
			if (!patchGrid.isTooSmall(nPatchPointsMin))
			{
				patchGrid.finish();
				gridList.add(patchGrid);
			}
		}
		StsException.systemDebug(this, "addRemainingGrids", "added " + patchGrids.length + " grids remaining.");
	}

	public void setCroppedBoundingBox(StsCroppedBoundingBox croppedBoundingBox)
	{
		this.croppedBoundingBox = croppedBoundingBox;
		croppedBoundingBox.setCroppedBoxRange();
	}

	/*
	* run the values calculation on the patches
	*/
	public void runCurvature(StsProgressPanel progressPanel, int filterSize, byte curveType, boolean runAllPatches)
	{
		this.filterSize = filterSize;

		if (runTimer)
		{
			timer = new StsTimer();
			timer.start();
		}
		StsPatchGrid[] runPatches = getRunCurvaturePatches(runAllPatches);
		int numPatches = runPatches.length;
		if (progressPanel != null)
			progressPanel.initialize(numPatches);
		int nValuePoints = 0;
		int nPoints = 0;
		float[] values = new float[nPointsTotal];
		double sum = 0;
		int progressUpdateInterval = Math.max(numPatches / 200, 1);
		int minNPoints = Math.min(filterSize * filterSize, StsQuadraticCurvature.minNPoints);
		for (int i = 0; i < numPatches; i++)
		{
			StsPatchGrid patch = runPatches[i];
			if (patch == null) continue;
			nPoints += patch.nPatchPoints;
			if (patch.computeCurvature(xInc, yInc, curveType, filterSize, minNPoints))
			{
				for (int row = patch.rowMin; row <= patch.rowMax; row++)
				{
					for (int col = patch.colMin; col <= patch.colMax; col++)
					{
						int ptCol = col - patch.colMin;
						int ptRow = row - patch.rowMin;
						float value = patch.values[ptRow][ptCol];
						if (value == StsPatchVolume.nullValue) continue;
						if (value == badCurvature || value == -badCurvature) continue;
						values[nValuePoints++] = value;
						sum += value;
					}
				}
			}

			if (progressPanel == null) continue;
			if (progressPanel.isCanceled())
			{
				progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
				clearPatches();
				return;
			}
			if (i % progressUpdateInterval == 0) progressPanel.setValue(i);
		}

		values = (float[]) StsMath.trimArray(values, nValuePoints);

		StsMessageFiles.infoMessage("Number values points: " + nValuePoints + " number of patch points " + nPointsTotal);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "Number values points: " + nValuePoints + " number of patch points " + nPointsTotal);


		double mean = sum / nValuePoints;
		int histogramDataInc = Math.max(1, nValuePoints / nHistogramValues);
		nHistogramValues = nValuePoints / histogramDataInc;
		histogramValues = new float[nHistogramValues + 1];
		double avgDev = 0;
		double sumSqr = 0;
		nHistogramValues = 0;
		for (int n = 0; n < nValuePoints; n++)
		{
			float value = values[n];
			double dif = value - mean;
			avgDev += Math.abs(dif);
			sumSqr += dif * dif;
			if (n % histogramDataInc == 0)
				histogramValues[nHistogramValues++] = value;
		}
		avgDev = avgDev / nValuePoints;
		double variance = (sumSqr) / nValuePoints;
		double stdDev = Math.sqrt(variance);
		StsMessageFiles.infoMessage("mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
		dataMin = (float) (mean - 2.0 * avgDev);
		dataMax = (float) (mean + 2.0 * avgDev);
		colorscale.setRange(dataMin, dataMax);
		StsMessageFiles.infoMessage("range set to +- 2.0*avg dev: " + dataMin + " to " + dataMax);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "range set to +- 2.0*std dev: " + dataMin + " to " + dataMax);

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

		if (runTimer) timer.stopPrint("Time to compute values for " + numPatches + " patches.");
	}

	private void clearPatches()
	{
		for (StsPatchGrid patch : rowSortedPatchGrids)
		{
			if (patch != null) patch.clear();
		}
	}

	private StsPatchGrid[] getRunCurvaturePatches(boolean runAllPatches)
	{
		if (runAllPatches || this.selectedPatchGrids == null) return rowSortedPatchGrids;
		else return selectedPatchGrids;
	}

	private ArrayList<TracePoints> getOtherTraces(TracePoints newTrace, TracePoints prevColTrace, TracePoints[] prevRowTraces)
	{
		ArrayList<TracePoints> otherTraces = new ArrayList<>();
		if (prevColTrace != null)
		{
			addOtherTrace(otherTraces, prevColTrace);
		}
		int col = newTrace.col;
		if (prevRowTraces != null)
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
		if (otherTrace.nTracePatchPoints == 0) return;
		otherTraces.add(otherTrace);
	}

	private StsPatchGrid getPatchGrid(int id)
	{
		StsPatchGrid patchGrid = rowGrids.get(id);
		// if patchGrid exists in rowGrids, then it has already been added there and deleted from prevRowGrids
		if (patchGrid != null)
		{
			if (StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
				StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
						" gotten from rowGrids at row: " + row + " col: " + col);
			return patchGrid;
		}
		// if patchGrid is not in rowGrids, then add it there and delete it from prevRowGrids
		else if (prevRowGrids != null)
		{
			patchGrid = prevRowGrids.get(id);
			if (patchGrid != null)
			{
				rowGrids.put(id, patchGrid);
				prevRowGrids.remove(id);
				if (StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
					StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
							" gotten and deleted from prevRowsGrids and added to rowGrids at row: " + row + " col: " + col);
				return patchGrid;
			}
		}
		StsException.systemError(this, "getPatchGrid", "Couldn't get patchGrid for id " + id + " at row: " + " col: " + col);
		return null;
	}

	/**
	 * This PatchGridSet is for the row before row just finished.
	 * If a grid in this prev row is disconnected (doesn't have same patch in row just finished),
	 * then either delete it if it is a small window, or add it to volume set.
	 */
	void processPrevRowGrids(int row)
	{
		if (row == 0) return;
		StsPatchGrid[] prevRowPatchGrids = prevRowGrids.values().toArray(new StsPatchGrid[0]);
		int nDisconnectedGrids = 0;
		for (StsPatchGrid patchGrid : prevRowPatchGrids)
		{
			boolean disconnected = patchGrid.isDisconnected(row);
			if (!disconnected) continue;
			if (patchGrid.isTooSmall(nPatchPointsMin))
				nSmallGridsRemoved++;
			else
			{
				patchGrid.finish();
				gridList.add(patchGrid);
				nDisconnectedGrids++;
			}
			prevRowGrids.remove(patchGrid.id);
		}
		StsException.systemDebug(this, "processPrevRowGrids", "prev row: " + (row - 1) + " added " + nDisconnectedGrids + " disconnected grids");
	}

	public StsColorscale getCurvatureColorscale()
	{
		return colorscale;
	}

	/* Draw any map edges on all 2d sections */
	public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped,
							   boolean xAxisReversed, boolean yAxisReversed)
	{
		if (!getIsVisible()) return;

		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		StsColor drawColor = StsColor.BLACK; //getStsColor();
		gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
		drawColor.setGLColor(gl);

		if (dirNo == StsCursor3d.XDIR) /* map edge is along a col	*/
		{
			if (!checkColSortedPatchGrids()) return;
			int col = getNearestColCoor(dirCoordinate);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glShadeModel(GL.GL_SMOOTH);

			float x = dirCoordinate;
			// int nFirst = -1;
			// int n = -1;

			for (StsPatchGrid patchGrid : colSortedPatchGrids)
			{
				if (patchGrid.colMin > col) break;
				// n++;
				if (patchGrid.colMax < col) continue;

				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawCol(gl, col, x, yMin, yInc, colorscale, false, displayCurvature);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				// if (nFirst == -1) nFirst = n;
			}
			gl.glEnable(GL.GL_LIGHTING);
		}
		else if (dirNo == StsCursor3d.YDIR)
		{
			if (!checkRowSortedPatchGrids()) return;
			int row = getNearestRowCoor(dirCoordinate);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			gl.glShadeModel(GL.GL_SMOOTH);

			float y = dirCoordinate;
			// int nFirst = -1;
			// int n = -1;
			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			{
				if (patchGrid.rowMin > row) break;
				// n++;
				if (patchGrid.rowMax < row) continue;
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawRow(gl, row, y, xMin, xInc, colorscale, false, displayCurvature, filter, boxFilterWidth);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				// if (nFirst == -1) nFirst = n;
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) break;
			}
			gl.glEnable(GL.GL_LIGHTING);
		}
	}

	/** Draw any map edges on section */
	public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate)
	{
		if (!getIsVisible()) return;
		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		try
		{
			if (dirNo == StsCursor3d.ZDIR)
			{
				if (getDisplaySurfs())
					// displayPatchesNearXYZCursors(glPanel3d);
					displayPatchesNearZCursor(glPanel3d, dirCoordinate);
				return;
			}
			if (dirNo == StsCursor3d.YDIR)
			{
				if (!checkRowSortedPatchGrids()) return;
				gl.glDisable(GL.GL_LIGHTING);
				gl.glShadeModel(GL.GL_SMOOTH);
				StsColor drawColor = StsColor.BLACK; //getStsColor();
				drawColor.setGLColor(gl);
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
				int row = getNearestRowCoor(dirCoordinate);
				if (row == -1) return;
				float xMin = getXMin();
				float xInc = getXInc();
				for (StsPatchGrid patchGrid : rowSortedPatchGrids)
				{
					if (patchGrid.rowMin > row) break;
					if (patchGrid.rowMax < row) continue;
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawRow(gl, row, dirCoordinate, xMin, xInc, colorscale, true, displayCurvature, setFilter(), getBoxFilterWidth());
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
					gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				}
			}
			else if (dirNo == StsCursor3d.XDIR)
			{
				if (!checkColSortedPatchGrids()) return;
				gl.glDisable(GL.GL_LIGHTING);
				gl.glShadeModel(GL.GL_SMOOTH);
				StsColor drawColor = StsColor.BLACK; //getStsColor();
				drawColor.setGLColor(gl);
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
				int col = getNearestColCoor(dirCoordinate);
				if (col == -1) return;
				float yMin = getYMin();
				float yInc = getYInc();
				for (StsPatchGrid patchGrid : colSortedPatchGrids)
				{
					if (patchGrid.colMin > col) break;
					if (patchGrid.colMax < col) continue;
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawCol(gl, col, dirCoordinate, yMin, yInc, colorscale, true, displayCurvature);
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
					gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				}
			}
		}
		catch (Exception e)
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
		if (cursor3d == null) return;
		GL gl = glPanel3d.getGL();
		for (int dir = 0; dir < 2; dir++)
		{
			float dirCoordinate = cursor3d.getCurrentDirCoordinate(dir);
			drawOnCursor3d(glPanel3d, dir, dirCoordinate);
		}
	}

	public void display(StsGLPanel glPanel)
	{
		if (!getDisplaySurfs() || selectedPatchGrids == null) return;
		GL gl = glPanel.getGL();

		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();

		initializePatchDraw(gl);
		boolean displayChildPatches = getPatchVolumeClass().getDisplayChildPatches();

		for (StsPatchGrid patchGrid : selectedPatchGrids)
		{
			if(patchGrid.parentGrid != null && displayChildPatches)
				patchGrid.parentGrid.drawPatchGrid(gl, displayChildPatches, displayCurvature, colorscale);
			else
				patchGrid.drawPatchGrid(gl, displayChildPatches, displayCurvature, colorscale);
		}
		if (getDisplayVoxels())
		{
			displayVoxels(glPanel);
		}
	}

	public void displayVoxelsCursor(StsGLPanel glPanel3d, StsPoint[] points, boolean is3d)
	{
		//System.out.println("Display Voxels");
		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		StsGridPoint point1 = new StsGridPoint(points[0], this);
		StsGridPoint point2 = new StsGridPoint(points[3], this);
		int sameRow = StsGridPoint.getSameRow(point1, point2); // if not -1, this is the row these two points are on
		if (sameRow != -1)
		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
			gl.glColor4f(1.f, 1.f, 1.f, 1.f);

			float y;

			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			{
				int row1 = sameRow - 4;
				int row2 = sameRow + 4;
				y = yMin + (yInc * row1);
				for (int n = row1; n <= row2; n++)
				{
					//if (n == 146)
					//System.out.println("draw vox row "+n+" "+patchGrid.colMin+" "+patchGrid.colMax+" "+y);
					patchGrid.drawRow(gl, n, y, xMin, xInc, colorscale, is3d, displayCurvature, setFilter(), getBoxFilterWidth());
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
		if (gl == null) return;

		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			gl.glColor4f(1.f, 1.f, 1.f, 1.f);
			float xMin = getXMin();
			float xInc = getXInc();
			float yMin = getYMin();
			float yInc = getYInc();
			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
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
		if (gl == null) return;

		initializePatchDraw(gl);
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
		{
			if (patchGrid.isPatchGridNearZCursor(z))
				patchGrid.draw(gl, displayCurvature, colorscale);
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
		if (gl == null) return;

		float x = cursor3d.getCurrentDirCoordinate(XDIR);
		float y = cursor3d.getCurrentDirCoordinate(YDIR);
		float z = cursor3d.getCurrentDirCoordinate(ZDIR);
		StsPoint cursorPoint = new StsPoint(x, y, z);
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		if (currentCursorPoint != null && currentCursorPoint.equals(cursorPoint) && cursorPointPatch != null)
		{
			drawPatch(cursorPointPatch, displayCurvature, gl);
			return;
		}

		cursorPoint = null;
		currentCursorPoint = null;

		int volumeRow = getNearestRowCoor(y);
		if (volumeRow == -1) return;
		int volumeCol = getNearestColCoor(x);
		if (volumeCol == -1) return;
		int slice = getNearestSliceCoor(z);
		if (slice == -1) return;
		float dzPatch = largeFloat;
		cursorPointPatch = null;
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
		{
			if (patchGrid == null) continue;
			// if(patchGrid.values == null) continue;
			if (patchGrid.rowMin > volumeRow) break;
			if (patchGrid.rowMax >= volumeRow)
			{
				if (patchGrid.colMin <= volumeCol && patchGrid.colMax >= volumeCol)
				{
					float dz = patchGrid.getZDistance(volumeRow, volumeCol, z);
					if (dz < dzPatch)
					{
						cursorPointPatch = patchGrid;
						dzPatch = dz;
						currentCursorPoint = new StsPoint(x, y, z);
					}
				}
			}
		}
		if (cursorPointPatch == null) return;

		drawPatch(cursorPointPatch, displayCurvature, gl);

		return;
	}

	private void drawPatch(StsPatchGrid patch, boolean displayCurvature, GL gl)
	{
		initializePatchDraw(gl);
		patch.draw(gl, displayCurvature, colorscale);
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
		for (n = 0; n < nPatchGrids; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			if (patchGrid == null) continue;
			rowMax = n - 1;
			if (rowMin == -1)
				if (patchGrid.rowMin <= volumeRow && patchGrid.rowMax >= volumeRow) rowMin = n;
				else if (patchGrid.rowMin > volumeRow)
					break;
		}
		if (rowMin == -1)
			return new int[]{0, 0};
		else
			return new int[]{rowMin, rowMax};
	}

	public int[] getPatchRangeForCol(int col)
	{
		int colMin = -1;
		int colMax = -1;
		int nPatchGrids = colSortedPatchGrids.length;
		for (int n = 0; n < nPatchGrids; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			if (patchGrid == null) continue;
			if (colMin == -1)
			{
				if (patchGrid.colMin <= col && patchGrid.colMax >= col) colMin = n;
			}
			else if (patchGrid.colMin > col)
			{
				colMax = n - 1;
				break;
			}
		}
		if (colMin == -1 || colMax == -1)
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
			for (int nPatch = nPatchMin; nPatch <= nPatchMax; nPatch++)
			{
				StsPatchGrid patchGrid = rowSortedPatchGrids[nPatch];
				if (patchGrid == null) continue;
				if (patchGrid.values == null) continue;
				int patchRow = volRow - patchGrid.rowMin;
				int patchCol = volCol - patchGrid.colMin;
				if (patchRow < 0 || patchCol < 0 || patchRow >= patchGrid.nRows || patchCol >= patchGrid.nCols)
					continue;
				float[][] pointsZ = patchGrid.getPointsZ();
				if (pointsZ == null) continue;
				float z = pointsZ[patchRow][patchCol];
				if (z == StsParameters.nullValue) continue;
				float val = patchGrid.values[patchRow][patchCol];
				if (val == nullValue) continue;
				int slice = getNearestSliceCoor(z);
				buffer[slice] = val;
				traceLoaded = true;
			}
			return traceLoaded;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "getTraceCurvature", e);
			return false;
		}
	}

	public int getNearestSliceCoor(float z)
	{
		int slice = Math.round((z - zMin) / interpolatedZInc);
		if (slice < 0 || slice >= nInterpolatedSlices) return -1;
		return slice;
	}

	public int getPatchPointIndex(PatchPoint patchPoint)
	{
		return patchPoint.getRow() * nCols + patchPoint.getCol();
	}

	public String toString()
	{
		return name;
	}

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
		if (objectPanel == null)
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
		if (this.displaySurfs == displaySurfs)
			return;
		this.displaySurfs = displaySurfs;
		if(!displaySurfs) clearSelectedPatches();
		currentModel.win3dDisplayAll();
	}

	public boolean getDisplaySurfs()
	{
		return displaySurfs;
	}

	public void setDisplayVoxels(boolean displayVoxels)
	{
		if (this.displayVoxels == displayVoxels)
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
		if (colorscale == null) return;
		colorscale.setRange(dataMin, dataMax);
	}

	public void setDataMax(float max)
	{
		dataMax = max;
		if (colorscale == null) return;
		colorscale.setRange(dataMin, dataMax);
	}

	public void addRemoveSelectedPatch(StsCursorPoint cursorPoint)
	{
		float[] xyz = cursorPoint.point.v;
		int volumeRow = getNearestRowCoor(xyz[1]);
		int volumeCol = getNearestColCoor(xyz[0]);
		float z = xyz[2];
		StsPatchGrid selectedPatch = getNearestPatch(volumeRow, volumeCol, z);
		if (selectedPatch == null) return;

		float iline = getRowNumFromRow(volumeRow);
		float xline = getColNumFromCol(volumeCol);
		StsMessageFiles.logMessage("Picked patch: " + selectedPatch.getPatchTypeString() + " " + selectedPatch.toFamilyString() + " at iline: " + iline + " xline: " + xline + " z: " + z);
	/*
		if(cursorPoint.dirNo == StsCursor3d.YDIR)
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeRowCorrel(volumeRow, volumeCol));
        else //dirNo == XDIR
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeColCorrel(volumeRow, volumeCol));
    */
		// int nSheet = selectedPatch.nSheet;

		StsPatchGrid selectedPatchParent = selectedPatch.getParentGrid();
		String addRemove;
		boolean removePatch = StsMath.arrayContains(selectedPatchGrids, selectedPatchParent);
		if (removePatch)
			addRemove = " removed ";
		else
			addRemove = " added ";

		String gridsString = selectedPatchParent.getGridDescription();
		StsMessageFiles.logMessage("Picked patch: " + selectedPatch.toFamilyString() + addRemove + " " + gridsString);

		if (removePatch)
			selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayDeleteElement(selectedPatchGrids, selectedPatchParent);
		else
		{
			selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayAddElement(selectedPatchGrids, selectedPatchParent);
		}
		currentModel.win3dDisplayAll();
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
		for (int n = patchMin; n <= patchMax; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			float patchZ = patchGrid.getVolumePointZ(volumeRow, volumeCol);
			if (patchZ == nullValue) continue;
			float dz = Math.abs(z - patchZ);
			if (dz < nearestPatchZ)
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

class TracePoints
{
	StsPatchVolume patchVolume;
	/** volume row for this trace */
	int row;
	/** volume col for this trace */
	int col;
	/** array of PatchPoints for this trace of various pointTypes (min, max, +zero-crossing, -zero-crossing; can be false or missing) */
	PatchPoint[] tracePatchPoints = new PatchPoint[0];
	/** length of tracePatchPoints array */
	int nTracePatchPoints;
	/** offset to first plus-zero-crossing in tracePatchPoints array */
	int zeroPlusOffset;
	/** a half-wave length window for each legitimate window (acceptable pointType) */
	CorrelationWindow[] windows;
	/** number of windows in this trace; one for which tracePatchPoint */
	int nWindows;
	/** double-linked list of connections between this trace and prevCol */
	ConnectionList colConnections;
	/** double-linked list of connections between this trace and prevCol */
	ConnectionList rowConnections;
	/** closest connections between this trace and the prevColTrace */
	CorrelationWindow[] colCloseConnectWindows;
	/** closest connections between this trace and the prevRowTrace */
	CorrelationWindow[] rowCloseConnectWindows;
	/**
	 * From the traceValues array, create an array of PatchPoint for each traceValue which qualifies as a pickType (all, min, max, +/- zero-crossing.
	 * This will be a sequential subset of the values array with nTracePatchPoints in this tracePatchPoints array
	 * @param patchVolume
	 * @param row volume row of this trace
	 * @param col volume col of this trace
	 * @param traceValues original seismic values for this trace
	 */
	private TracePoints(StsPatchVolume patchVolume, int row, int col, float[] traceValues) throws StsException
	{
		this.patchVolume = patchVolume;
		this.row = row;
		this.col = col;
		// tracePoints - uniform cubic interpolation of traceValues
		float[] tracePoints = StsTraceUtilities.computeCubicInterpolatedPoints(traceValues, patchVolume.nInterpolationIntervals);
		float z = patchVolume.croppedBoundingBox.zMin;
		if (tracePoints == null) return;
		int nTracePoints = tracePoints.length;
		tracePatchPoints = new PatchPoint[nTracePoints];
		int nTracePatchPoint = 0;
		byte[] tracePointTypes = StsTraceUtilities.getPointTypes(tracePoints);

		// create tracePoints from values and pointTypes for legitimate events (zero+, max, zero-, min)
		// count the number missing so we have the final array length with missing points added
		int nTotalMissing = 0;
		int nMissing;
		PatchPoint prevPoint, nextPoint = null;

		try
		{
			for (int n = 0; n < nTracePoints; n++, z += patchVolume.interpolatedZInc)
			{
				byte tracePointType = tracePointTypes[n];
				if (StsTraceUtilities.isMaxMinZeroOrFalseMaxMin(tracePointType))
				{

					prevPoint = nextPoint;
					nextPoint = new PatchPoint(this, n, z, tracePoints[n], tracePointType, nTracePatchPoint);
					tracePatchPoints[nTracePatchPoint] = nextPoint;
					nTracePatchPoint++;
					nMissing = getNMissingPoints(prevPoint, nextPoint);
					nTotalMissing += nMissing;
				}
			}
			// insert any missing points so that we have a series of zero+, max, zero-, min points with any missing filled as ZPM, MXM, MNM, ZMM

			nTracePatchPoints = nTracePatchPoint;
			if (nTotalMissing > 0)
			{
				int nTotalPoints = nTracePatchPoints + nTotalMissing;
				PatchPoint[] addedPoints = new PatchPoint[nTotalPoints];
				nextPoint = tracePatchPoints[0];
				nTracePatchPoint = 0;
				for (int n = 1; n < nTracePatchPoints; n++)
				{
					prevPoint = nextPoint;
					nextPoint = tracePatchPoints[n];

					addedPoints[nTracePatchPoint] = prevPoint.resetIndex(nTracePatchPoint);
					nTracePatchPoint++;
					nMissing = getNMissingPoints(prevPoint, nextPoint);
					if( nMissing > 0)
					{
						byte pointTypeStart = StsTraceUtilities.pointTypesAfter[prevPoint.pointType];
						byte pointTypeEnd = StsTraceUtilities.pointTypesBefore[nextPoint.pointType];
						byte missingType = pointTypeStart;
						for (int i = 0; i < nMissing; i++, missingType++)
						{
							if (missingType > 4) missingType -= 4;
							addedPoints[nTracePatchPoint] = new PatchPoint(prevPoint, missingType, nTracePatchPoint);
							nTracePatchPoint++;
						}
					}
				}
				addedPoints[nTracePatchPoint] = nextPoint.resetIndex(nTracePatchPoint);
				nTracePatchPoint++;
				//if (nTracePatchPoint != nTotalPoints)
				//	StsException.systemError(this, "new TracePoints", " nTotalPoints " + nTotalPoints + " not equal to nTracePatchPoints " + nTracePatchPoint);
				tracePatchPoints = addedPoints;
				nTracePatchPoints = nTracePatchPoint;
			}
			else
				tracePatchPoints = (PatchPoint[]) StsMath.trimArray(tracePatchPoints, nTracePatchPoints);

			if(tracePatchPoints.length < 2)  throw new StsException("TracePoints.constructor()", "Less than 2 tracePatchPoints.");

			zeroPlusOffset = StsTraceUtilities.zeroPlusOffset[tracePatchPoints[0].pointType];
			constructTraceWindows();
		}
		catch (StsException stse)
		{
			throw new StsException("TracePoints.constructor()", stse.getMessage());
		}
		catch (Exception e)
		{
			StsException.outputWarningException(TracePoints.class, "constructor", e);
			throw new StsException("TracePoints.constructor()", e.getMessage());
		}
	}

	public static TracePoints constructor(StsPatchVolume patchVolume, int row, int col, float[] traceValues)
	{
		try
		{
			return new TracePoints(patchVolume, row, col, traceValues);
		}
		catch(StsException e)
		{
			return null;
		}
	}

	static public boolean pointsNotInSequence(PatchPoint prevPoint, PatchPoint nextPoint)
	{
		if(prevPoint == null || nextPoint == null) return false;
		byte prevType = StsTraceUtilities.coercedPointTypes[prevPoint.pointType];
		byte nextType = StsTraceUtilities.coercedPointTypes[nextPoint.pointType];
		return StsTraceUtilities.pointTypesAfter[prevType] != nextType;
	}

	static int getNMissingPoints(PatchPoint prevPoint, PatchPoint nextPoint)
	{
		if(prevPoint == null || nextPoint == null) return 0;
		byte prevType = StsTraceUtilities.coercedPointTypes[prevPoint.pointType];
		byte nextType = StsTraceUtilities.coercedPointTypes[nextPoint.pointType];
		if(StsTraceUtilities.pointTypesAfter[prevType] == nextType) return 0;
		byte pointTypeStart = StsTraceUtilities.pointTypesAfter[prevPoint.pointType];
		byte pointTypeEnd = StsTraceUtilities.pointTypesBefore[nextPoint.pointType];
		return StsTraceUtilities.getNumPointTypesBetweenInclusive(pointTypeStart, pointTypeEnd);
	}

	/**
	 * creates correlated connections between this trace and traces at prev row & same col and prev col & same row
	 * @param prevColTrace prev trace in same col, prev row
	 * @param prevRowTrace prev trace in same row, prev col
	 */
	void connectWindows(TracePoints prevColTrace, TracePoints prevRowTrace)
	{
		if (prevRowTrace == null && prevColTrace == null) return;

		if (windows == null || windows.length < 1) return;

		// The ConnectionLists for this trace is initialized with top and bot inactive connections to prev row and col traces.
		// These inactive connections are used to limit the search process.
		// a connection is from the first window in this patchPointsList to the first window in the corresponding row or col trace
		if(prevColTrace != null)
		{
			colConnections = new ConnectionList(this, prevColTrace);
			prevColTrace.colConnections = colConnections;
		}
		if(prevRowTrace != null)
		{
			rowConnections = new ConnectionList(this, prevRowTrace);
			prevRowTrace.rowConnections = rowConnections;
		}
		// create initial guesses of connections from the windows on this trace to windows on prevCol and prevRow traces
		// guesses are the closest windows vertically on the other traces to the window on this trace
		if(prevColTrace != null)
			createClosestConnectWindow(this, prevColTrace, false);
		if(prevRowTrace != null)
			createClosestConnectWindow(this, prevRowTrace, true);
		// Iterate from maxCorrelation down to minCorrelation in nIterations steps.  At each iteration,
		// make connections from prevCol and prevRow traces to this trace */
		for (int iter = 0; iter < patchVolume.nIterations; iter++)
		{
			StsPatchVolume.setIterLabel(iter);
			connectWindows(patchVolume, prevColTrace, prevRowTrace, iter);
		}
	}

	/** This prevents crossing connections.
	 *  Create lists of closest points of same type from trace to otherTrace and back.
	 *  If connections are identical, then retain in the trace->otherTrace list; if not null out.
	 * @param trace connections will be from this trace to otherTrace
	 * @param prevTrace connections list back will be compared with trace list
	 * @param isRow indicates trace and otherTrace are on the same row (other trace is prevCol)
	 */
	static void createClosestConnectWindow(TracePoints trace, TracePoints prevTrace, boolean isRow)
	{
		CorrelationWindow window, prevWindow, backWindow;
		int n;
		CorrelationWindow[] windows = trace.windows;
		int nWindows = trace.nWindows;
		CorrelationWindow[] connectWindows = trace.createClosestConnectWindows(prevTrace);
		CorrelationWindow[] otherConnectWindows = prevTrace.createClosestConnectWindows(trace);
		try
		{
			CorrelationWindow lastWindow = null;
			for (n = 0; n < nWindows; n++)
			{
				window = windows[n];
				prevWindow = connectWindows[n];

				if(prevWindow == null)
					continue;
				int otherIndex = prevWindow.windowIndex;
				backWindow = otherConnectWindows[otherIndex];
				if (backWindow != window)
				{
					connectWindows[n] = null;
					otherConnectWindows[otherIndex] = null;
				}
				else if(lastWindow != null && prevWindow.getCenterPoint().slice <= lastWindow.getCenterPoint().slice)
					connectWindows[n] = null;
				else
					lastWindow = prevWindow;
			}
			// check for crossing connections; remove connection if it crosses
			if (isRow)
			{
				trace.rowCloseConnectWindows = connectWindows;
				prevTrace.rowCloseConnectWindows = otherConnectWindows;
			}
			else
			{
				trace.colCloseConnectWindows = connectWindows;
				prevTrace.colCloseConnectWindows = otherConnectWindows;
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(TracePoints.class, "createClosestConnectWindow", e);
		}
	}

	static void checkCrossings(CorrelationWindow[] connectWindows, int nextIndex)
	{
		if(nextIndex == 0) return;

		int prevSlice = -1, slice, nextSlice;
		if(nextIndex > 1)
			prevSlice = connectWindows[nextIndex - 2].getCenterPoint().slice;
		slice = connectWindows[nextIndex - 1].getCenterPoint().slice;
		nextSlice = connectWindows[nextIndex].getCenterPoint().slice;
		if(prevSlice >=slice || nextSlice <= slice)
			connectWindows[nextIndex-1] = null;
	}

	CorrelationWindow[] createClosestConnectWindows(TracePoints prevTrace)
	{
		CorrelationWindow[] closeConnectWindows = new CorrelationWindow[nWindows];

		// assign closest prevWindow to this window regardless of pointType
		CorrelationWindow[] prevWindows = prevTrace.windows;
		int nPrevWindows = prevTrace.nWindows;
		CorrelationWindow prevWindowAbove = prevWindows[0];
		CorrelationWindow prevWindowBelow = prevWindows[1];
		int otherNextIndex = 2;
		for (int i = 0; i < nWindows; i++)
		{
			CorrelationWindow window = windows[i];
			if (window.isBelowOrEqual(prevWindowAbove) && window.isAboveOrEqual(prevWindowBelow))
			{
				CorrelationWindow prevWindow = window.getClosestWindow(prevWindowAbove, prevWindowBelow);
				closeConnectWindows[i] = prevWindow;
			}
			else if (window.isAboveOrEqual(prevWindowAbove)) continue;
			else // window is below prevWindowBelow, so move prevWindows down
			{
				while (window.isBelowOrEqual(prevWindowBelow) && otherNextIndex < nPrevWindows)
				{
					prevWindowAbove = prevWindowBelow;
					prevWindowBelow = prevWindows[otherNextIndex++];
				}
				if (window.isBelowOrEqual(prevWindowAbove) && window.isAboveOrEqual(prevWindowBelow))
				{
					CorrelationWindow prevWindow = window.getClosestWindow(prevWindowAbove, prevWindowBelow);
					closeConnectWindows[i] = prevWindow;
				}
			}
		}
		// now adjust guesses to nearest prevWindow of this pointType
		for (int i = 0; i < nWindows; i++)
		{
			CorrelationWindow window = windows[i];
			CorrelationWindow prevWindow = closeConnectWindows[i];
			if(prevWindow == null) continue;
			int pointTypeOffset = getPointTypeDif(window, prevWindow);
			if(pointTypeOffset == 0) continue; // we have the correct type: so don't adjust
			if (pointTypeOffset >= 2) pointTypeOffset -= 4;
			int prevWindowIndex = prevWindow.windowIndex + pointTypeOffset;
			if (prevWindowIndex < 0)
				prevWindowIndex += 4;
			else if(prevWindowIndex >= nPrevWindows)
				prevWindowIndex -= 4;
			prevWindow = prevWindows[prevWindowIndex];
			closeConnectWindows[i] = prevWindowIsInside(window, prevWindow);
		}
		// now compute correlation factors
		for (int i = 0; i < nWindows; i++)
		{
			CorrelationWindow prevWindow = closeConnectWindows[i];
			if(prevWindow == null) continue;
			windows[i].computeCorrelation(prevWindow, patchVolume.autoCorInc);
		}
		return closeConnectWindows;
	}

	private CorrelationWindow prevWindowIsInside(CorrelationWindow window, CorrelationWindow prevWindow)
	{
		int prevWindowCenterSlice = prevWindow.getCenterPoint().slice;
		int prevPointSlice = getPrevWindow(window).getCenterPoint().slice;
		if(prevWindowCenterSlice < prevPointSlice) return null;
		int nextPointSlice = getNextWindow(window).getCenterPoint().slice;
		if(prevWindowCenterSlice > nextPointSlice) return null;
		return prevWindow;
	}

	private void connectWindows(StsPatchVolume patchVolume, TracePoints prevColTrace, TracePoints prevRowTrace, int iter)
	{
		CorrelationWindow matchingWindow, backMatchingWindow;
		Connection colConnection, rowConnection;
		CorrelationWindow otherClosestWindow, closestWindow;
		CorrelationWindow connectedWindowAbove, connectedWindowBelow;
		CorrelationWindow window;

		try
		{
			reinitializeTraceIndices(prevRowTrace, prevColTrace);
			for (int n = 0; n < nWindows; n++)
			{
				window = windows[n];

				if(StsPatchVolume.debug && StsPatchGrid.debugPoint && (StsPatchGrid.doDebugPoint(window.getCenterPoint())))
					StsException.systemDebug(this, "connectWindows", StsPatchVolume.iterLabel + " MATCH THIS WINDOW: " + window.toString());

				colConnection = null;
				if (prevColTrace != null && !window.hasColConnection())
				{
					otherClosestWindow = colCloseConnectWindows[n];
					connectedWindowAbove = colConnections.connectionAbove.prevWindow;
					connectedWindowBelow = colConnections.connectionBelow.prevWindow;
					matchingWindow = TracePoints.connectWindows(patchVolume, window, otherClosestWindow, prevColTrace, connectedWindowAbove, connectedWindowBelow, iter);
					if (matchingWindow != null && patchVolume.checkBackMatch)
					{
						closestWindow = prevColTrace.colCloseConnectWindows[matchingWindow.windowIndex];
						connectedWindowAbove = prevColTrace.colConnections.connectionAbove.window;
						connectedWindowBelow = prevColTrace.colConnections.connectionBelow.window;
						backMatchingWindow = TracePoints.connectWindows(patchVolume, matchingWindow, closestWindow, this, connectedWindowAbove, connectedWindowBelow, iter);
						if (backMatchingWindow != null && backMatchingWindow != window && backMatchingWindow.stretchCorrelation >= matchingWindow.stretchCorrelation)
							matchingWindow = null;
					}
					if (matchingWindow != null)
						colConnection = new Connection(matchingWindow, window);
				}
				rowConnection = null;
				if (prevRowTrace != null && !window.hasRowConnection())
				{
					otherClosestWindow = rowCloseConnectWindows[n];
					connectedWindowAbove = rowConnections.connectionAbove.prevWindow;
					connectedWindowBelow = rowConnections.connectionBelow.prevWindow;
					matchingWindow = TracePoints.connectWindows(patchVolume, window, otherClosestWindow, prevRowTrace, connectedWindowAbove, connectedWindowBelow, iter);
					if (matchingWindow != null && patchVolume.checkBackMatch)
					{
						closestWindow = prevRowTrace.rowCloseConnectWindows[matchingWindow.windowIndex];
						connectedWindowAbove = prevRowTrace.rowConnections.connectionAbove.window;
						connectedWindowBelow = prevRowTrace.rowConnections.connectionBelow.window;
						backMatchingWindow = TracePoints.connectWindows(patchVolume, matchingWindow, closestWindow, this, connectedWindowAbove, connectedWindowBelow, iter);
						if (backMatchingWindow != null && backMatchingWindow != window && backMatchingWindow.stretchCorrelation >= matchingWindow.stretchCorrelation)
							matchingWindow = null;
					}
					if (matchingWindow != null)
						rowConnection = new Connection(matchingWindow, window);
				}
				// Check if we have a cycle skip.
				// If the rowConnection.prevPoint overlaps the currentGrid at this point which was connected via the colConnection,
				// then we have a cycle skip. If we do, allow either the rowConnection or the colConnection, which ever has
				// the highest correlation, but not both. If we use the rowConnection and the colConnection already exists, we
				// need to delete the colConnection from the connected point..

				// if colConnection is not null, this is a new colConnection
				// Set currentColConnection to either this new colConnection or the existing colConnection
				// The rowConnection must be new as

				Connection currentRowConnection = rowConnection;
				if(currentRowConnection != null)
					currentRowConnection = window.getRowConnection();
				Connection currentColConnection = colConnection;
				if(currentColConnection != null)
					currentColConnection = window.getColConnection();
				if(currentRowConnection != null && currentColConnection != null)
				{
					StsPatchGrid colPatchGrid = currentColConnection.getPrevPatchGrid();
					StsPatchGrid rowPatchGrid = currentRowConnection.getPrevPatchGrid();
					if(colPatchGrid != null && rowPatchGrid != null && colPatchGrid != rowPatchGrid)
					{
						if(rowConnectionIsBetter(currentRowConnection, currentColConnection))
						{
							if(rowConnection == null) // this must be an existing connections
								window.deleteRowConnection();
							else // rowConnection != null: a new connection, so remove it
								rowConnection = null;

						}
						else // eliminate the col connection
						{
							if (colConnection == null)
								window.deleteColConnection();
							else
								colConnection = null;
						}
					}
				}
				processNewConnections(window, colConnection, rowConnection);
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "connectWindows(vol,trace,trace,iter)", e);
		}
	}

	static boolean rowConnectionIsBetter(Connection rowConnection, Connection colConnection)
	{
		CorrelationWindow rowWindow = rowConnection.window;
		CorrelationWindow colWindow = colConnection.window;

		if(rowWindow.amplitudeRatio > 2*colWindow.amplitudeRatio) return true;
		if(rowWindow.amplitudeRatio < 0.5*colWindow.amplitudeRatio) return false;

		return rowWindow.stretchCorrelation > colWindow.stretchCorrelation;
	}

	/** from the current connectionAbove for this connectionList, find index offset to this window and apply it to the otherTrace window
	 *  to get the middle matchingWindow.  Check connections to otherTrace windows above and below.
	 *
	 * @param patchVolume
	 * @param newWindow newWindow we want to match on this otherTrace
	 * @param centerOtherWindow center candidate window for matching
	 * @param otherTrace other trace to which we want to make a connection from this trace
	 * @param otherConnectWindowAbove otherTrace windows on connections above
	 * @param otherConnectWindowBelow otherTrace windows on connections below
	 * @param iter iteration we are on
	 * @return best connection between newWindow and a matchingWindow on this trace; return null of non exists or don't qualify
	 */
	static private CorrelationWindow connectWindows(StsPatchVolume patchVolume, CorrelationWindow newWindow,
													CorrelationWindow centerOtherWindow, TracePoints otherTrace,
													CorrelationWindow otherConnectWindowAbove, CorrelationWindow otherConnectWindowBelow,
													int iter)
	{
		// matchingWindow we wish to find
		CorrelationWindow matchingWindow = null;
		// index of candidate centerOtherWindow
		int centerOtherWindowIndex;
		// indexes for the two other matchingWindow candidates above and below the center (offsets of -4 and +4)
		int aboveWindowIndex, belowWindowIndex;
		// candidate matching windows above and below
		CorrelationWindow aboveOtherWindow, belowOtherWindow;
		// having found a matchingWindow, try matching it back to the newTrace to see if we find a different connection with a better stretchCorrelation
		// if we do find it, we ignore this match completely and let the search find it directly (rather than backMatching to find it)
		CorrelationWindow backMatchingWindow;

		try
		{
			if(centerOtherWindow == null) return null;

			float correlation = patchVolume.stretchCorrelations[iter];
			float minAmplitudeRatio = patchVolume.minAmplitudeRatio;
			// set a penalty except for the last iteration
			float correlPenalty = 0.0f;
			if(iter < patchVolume.nIterations-1)
				correlPenalty = patchVolume.autoCorInc;

			// centerOtherWindow must be between bounding connections above and below and cannot cross
			// nextWindow is already between them, so move centerOtherWindow up or down to be between as well
			// new window selected must be of same type so move is +/- 4 index

			if(centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
			{
				if(StsPatchVolume.debugConnectCloseOnly) return null;
				//CorrelationWindow prevWindow = centerOtherWindow;
				while(centerOtherWindow != null && centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
				{
					int index = centerOtherWindow.windowIndex + 4;
					if(index >= otherTrace.nWindows)
						return null;
					//prevWindow = centerOtherWindow;
					centerOtherWindow = otherTrace.windows[index];
					if(centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
						return null;
				}
				//centerOtherWindow = newWindow.getClosestWindow(prevWindow, centerOtherWindow);
			}
			else if(centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
			{
				if(StsPatchVolume.debugConnectCloseOnly) return null;
				//CorrelationWindow prevWindow = centerOtherWindow;
				while(centerOtherWindow != null && centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
				{
					int index = centerOtherWindow.windowIndex - 4;
					if(index < 0)
						return null;
					//prevWindow = centerOtherWindow;
					centerOtherWindow = otherTrace.windows[index];
					if(centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
						return null;
				}
				//centerOtherWindow = newWindow.getClosestWindow(prevWindow, centerOtherWindow);
			}
			// centerOtherWindow is between above and below connection points on otherTrace, so compute stretchCorrelation with window on this trace
			// this centerOtherWindow has already been determined to be the closest if their are two bracketing windows (@see
			if (newWindow.correlationOK(correlation, correlPenalty, minAmplitudeRatio))
			{
				matchingWindow = centerOtherWindow;
				correlation = centerOtherWindow.stretchCorrelation;
			}

			if(StsPatchVolume.debugConnectCloseOnly)
				return matchingWindow;

			// try to make a match with prevWindow above centerOtherWindow
			centerOtherWindowIndex = centerOtherWindow.windowIndex;
			aboveWindowIndex = centerOtherWindowIndex - 4;
			if (aboveWindowIndex > otherConnectWindowAbove.windowIndex) // index must be below connectionAbove
			{
				aboveOtherWindow = otherTrace.windows[aboveWindowIndex];
				if (newWindow.correlationOK(correlation, correlPenalty, minAmplitudeRatio))
				{
					matchingWindow = aboveOtherWindow;
					correlation = aboveOtherWindow.stretchCorrelation;
				}
			}

			belowWindowIndex = centerOtherWindowIndex + 4;
			if (belowWindowIndex < otherConnectWindowBelow.windowIndex)
			{
				belowOtherWindow = otherTrace.windows[belowWindowIndex];
				if (newWindow.correlationOK(correlation, correlPenalty, minAmplitudeRatio))
				{
					matchingWindow = belowOtherWindow;
					correlation = belowOtherWindow.stretchCorrelation;
				}
			}
			return matchingWindow;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(TracePoints.class, "connectWindows", e);
			return null;
		}
	}

	static CorrelationWindow getCenterOtherWindow(CorrelationWindow newWindow, TracePoints otherTrace, ConnectionList connectionList)
	{
		// center candidate window for matching
		CorrelationWindow centerOtherWindow;
		// bounding connections above and below which this connection cannot cross
		Connection connectionAbove, connectionBelow;
		// trace this new window is on
		TracePoints newTrace;
		// points on connections above and below
		PatchPoint newPointAbove, otherPointAbove, otherPointBelow;
		// offset from connectionAbove.point.traceIndex to the newWindow.centerPoint.traceIndex
		int newWindowIndexOffset;
		// index of candidate centerOtherWindow
		int centerOtherWindowIndex;
		// offset from the pointType at the parallel offset on the otherTrace from connectionAbove to the pointType we want on otherTrace
		int pointTypeOffset;
		// total offset including windowIndexOffset and pointTypeOffset
		int offset;

		connectionAbove = connectionList.connectionAbove;
		newTrace = newWindow.getCenterPoint().tracePoints;
		newPointAbove = connectionAbove.getWindowCenterPoint(newTrace);
		// newWindowIndexOffset is offset from connectionAbove to this newWindow.centerPoint
		newWindowIndexOffset = newWindow.windowIndex - newPointAbove.traceIndex;
		otherPointAbove = connectionAbove.getWindowCenterPoint(otherTrace);
		// compute offset to point on otherTrace which is parallel to connectionAbove and displaced down to newWindow
		// so the point is located newWindowIndexOffset below the otherPointAbove.traceIndex
		centerOtherWindowIndex = otherPointAbove.traceIndex + newWindowIndexOffset;
		if(centerOtherWindowIndex >= otherTrace.nWindows)
			centerOtherWindowIndex -= 4;
		centerOtherWindow = otherTrace.windows[centerOtherWindowIndex];
		// depending on the type here, we want to further offset from the type at this point to the desired type
		// which is equal to the numerical difference in type byte values
		// offset is always positive and is 0,1,2, or 3; if 0 or 1, use this offset; otherwise for 2, or 3, offset -2 or -1 respectively (offset-4)
		pointTypeOffset = getPointTypeDif(newWindow, centerOtherWindow);
		if (pointTypeOffset >= 2) pointTypeOffset -= 4;
		// the total offset is the shift from the connetionAbove (newWindowIndexOffset) plus the pointTypeOffset
		offset = newWindowIndexOffset + pointTypeOffset;
		if (offset < 1) offset += 4;
		centerOtherWindowIndex = otherPointAbove.traceIndex + offset;
		connectionBelow = connectionList.connectionBelow;
		otherPointBelow = connectionBelow.getWindowCenterPoint(otherTrace);

		if(centerOtherWindowIndex < otherPointBelow.traceIndex)
			return otherTrace.windows[centerOtherWindowIndex];

		while(centerOtherWindowIndex > otherPointAbove.traceIndex)
		{
			centerOtherWindowIndex -= 4;
			if(centerOtherWindowIndex < otherPointBelow.traceIndex)
				return otherTrace.windows[centerOtherWindowIndex];
		}
		return null;
	}

	static int getPointTypeDif(CorrelationWindow newWindow, CorrelationWindow prevWindow)
	{
		int typeDif = newWindow.getCenterPoint().pointType - prevWindow.getCenterPoint().pointType;
		if(typeDif < 0) typeDif += 4;
		return typeDif;
	}
	/**
	 * Given a newPatchPoint at newRow-newCol, which correlates with a prevPatchPoint at prevRow-prevCol which is possibly part of a patchGrid in the prevPatchGridsSet,
	 * combine these two points in the same patch.  The prevPatchPoint may be on the previous col (same row), or previous row (same col).
	 * If the previousPatchPoint is not part of an existing patchGrid (prevID == -1), then we will create a new patchGrid and add both points to it.
	 * If the previousPatchPoint is part of a patchGrid we will add the newPatchPoint to this patchGrid, unless the newPatchPoint already belongs to another patchGrid
	 * (this occurs when we first correlate with the previous column and find one patchGrid and also correlate with the previous row and find a different patchGrid).
	 */
	public Connection addPatchConnection(Connection connection, ConnectionList connectionList, boolean isRow)
	{
		StsPatchGrid patchGrid = null;
		CorrelationWindow window = connection.window;
		CorrelationWindow prevWindow = connection.prevWindow;
		float correlation = prevWindow.stretchCorrelation;
		// if (stretchCorrelation < patchVolume.minLinkCorrel) return null;

		if(connectionList.connectionsCross(connection)) return null;

		PatchPoint newPatchPoint = window.getCenterPoint();
		PatchPoint otherPatchPoint = prevWindow.getCenterPoint();
		double distance = Math.abs(otherPatchPoint.slice - newPatchPoint.slice);

		if(distance > 20)
			StsException.systemDebug(this, "addPatchConnection", "DISTANCE LARGE FOR Connection: " + connection.toString());

		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

		if(StsPatchVolume.debug && (StsPatchGrid.doDebugPoint(newPatchPoint) || StsPatchGrid.doDebugPoint(otherPatchPoint)))
			StsException.systemDebug(this, "addPatchConnection", StsPatchVolume.iterLabel + " CONNECT point: " +
					newPatchPoint.toString() + " TO point: " + otherPatchPoint.toString());

		// normally we can insert a new connectedPoint in the trace patchPointsList and split the connected interval;
		// but if we have cloned this new window and it is already connected, don't add/split the trace again
		//boolean splitIntervalOK = true;
		if (newPatchGrid == null)
		{
			if (otherPatchGrid == null) // prevPatchGrid doesn't exist, so create it and add otherPoint to it
			{
				patchGrid = StsPatchGrid.construct(patchVolume, newPatchPoint.getPointType(patchVolume.useFalseTypes));
				patchGrid.addPatchPoint(otherPatchPoint);
			}
			else // otherPatchGrid does exist, so use it
			{
				// if this newPatchPoint overlaps the otherPatchGrid, we can't add it;
				// So create a new patch and add a clone of the otherPatchPoint
				if (StsPatchVolume.debugCloneOK && otherPatchGrid.patchPointOverlaps(newPatchPoint)) // return null;
				{
					patchGrid = StsPatchGrid.construct(patchVolume, newPatchPoint.pointType);
					//PatchPoint clonedOtherPatchPoint = otherPatchPoint.cloneAndClear();
					//patchGrid.addPatchPoint(clonedOtherPatchPoint);
					//connection.setPrevPoint(clonedOtherPatchPoint);
					otherPatchGrid.combineChildGrids(patchGrid);
					//splitIntervalOK = false;
				}
				else // no overlap, so we will only need to add the newPatchPoint to it (below else)
					patchGrid = otherPatchGrid;
			}
			patchGrid.addPatchPoint(newPatchPoint);
		}
		else // newPatchGrid != null which means this window was just added to a patch from prevColTrace and the patchGrid would have been added to the rowGrids array
		{
			if (otherPatchGrid == null) // the otherPoint is not assigned to a patch; assign it to this one; don't add window to rowGrids unless it overlaps and addedGrid created
			{
				patchGrid = newPatchGrid;
				// otherPatchPoint doesn't have a patchGrid, but newPatchPoint does; try to add otherPatchPoint to newPatchGrid,
				// but if it overlaps, created a new patchGrid containing otherPatchPoint and a clone of newPatchPoint
				if (StsPatchVolume.debugCloneOK && patchGrid.patchPointOverlaps(otherPatchPoint))
				{
					patchGrid = StsPatchGrid.construct(patchVolume, patchGrid.patchType);
					//PatchPoint clonedNextPatchPoint = newPatchPoint.cloneAndClear();
					//patchGrid.addPatchPoint(clonedNextPatchPoint);
					//connection.setNextPoint(clonedNextPatchPoint);
					newPatchGrid.combineChildGrids(patchGrid);
				}
				patchGrid.addPatchPoint(otherPatchPoint);

				// patchGrid = patchGrid.checkAddPatchPoint(otherPatchPoint, newPatchPoint);
				//patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
				//checkAddPatchGridToRowGrids(patchGrid);
			}
			else if (otherPatchGrid.id == newPatchGrid.id) // otherPoint is already assigned to the same patch: addCorrelation
			{
				patchGrid = newPatchGrid;
				//if (patchGrid == null) return null; // patchGrid not found; systemDebugError was printed
				//checkAddPatchGridToRowGrids(patchGrid);
				//patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
			}
			// prevPoint and this window belong to different patches: merge newPatchGrid into prevPatchGrid and add connection (below)
			// if the grids don't overlap.
			// If they do overlap, then we merge as many points from the smaller grid to the larger (@see mergeOverlappingPatchGrids);
			// the connection is added in mergeOverlappingPatchGrids (two connections may have been generated,
			// @see moveNonOverlappingPointsTo for details); connection is fully-handled, so we don't redundantly handle connection below.
			else
			{
				if (StsPatchGrid.mergePatchPointsOK(otherPatchGrid, newPatchGrid))
				{
					patchGrid = patchVolume.mergePatchGrids(otherPatchPoint, newPatchPoint);
					// if patchGrid==null error occurred: systemError written in mergePatchGrids routine
				}
				else
				{
					patchGrid = patchVolume.checkAddOverlappingConnection(otherPatchPoint, newPatchPoint, connection);
					//	patchGrid = patchVolume.mergeOverlappingPatchGrids(otherPatchPoint, newPatchPoint, connection);
					patchVolume.checkAddPatchGridToRowGrids(patchGrid);
				}
			}
		}
		if (patchGrid == null) return null;

		// patchGrid is the grid this connection is to go on to;
		// if both points are on the patchGrid, then we can add the connection to the connection.window.centerPoint (nextPoint);
		// if points are on different grids, we have two cases from processing above (mergeOverlappingPatchGrids takes care of itself):
		//   1) nextPoint is on the patchGrid and prevPoint is on the otherGrid:
		//      replace the connection.prevPoint with  a cloneAndClear copy of the prevPoint
		//   2) prevPoint is on the patchGrid and nextPoint is on the otherGrid:
		//      replace the nextPoint.connection with a cloned connection which has the isMoved flag set to true;
		//      the connection.nextPoint will be replaced by a clone of the nextPoint added to the patchGrid
		// checkResetClonedPoints(connection, otherPatchPoint, newPatchPoint, isRow);
		processConnection(connection, patchGrid);
		// connection.addConnectionToPoint();
		// ??? if merging overlapped grids, we have two adjust grids; do we need to call checkAddPatchGridToRowGrids for the other as well?
		patchVolume.checkAddPatchGridToRowGrids(patchGrid);
		//if (splitIntervalOK)
		return connection;
		//else
		//	return null;
	}

	/** patchGrid is the grid this connection is to go on to;
	 * if both points are on the patchGrid, then we can add the connection to the connection.window.centerPoint (nextPoint);
	 * if points are on different grids, we have two cases from processing above:
	 *   1) nextPoint is on the patchGrid and prevPoint is on the otherGrid:
	 *      replace the connection.prevPoint with  a cloneAndClear copy of the prevPoint added to the patchGrid
	 *   2) prevPoint is on the patchGrid and nextPoint is on the otherGrid:
	 *      replace the nextPoint.connection with a cloned connection which has the isMoved flag set to true;
	 *      the connection.nextPoint will be replaced by a clearAndClone of the nextPoint added to the patchGrid
	 * @param connection connection to be processed
	 * @param patchGrid grid this connection is to go on to
	 */
	void processConnection(Connection connection, StsPatchGrid patchGrid)
	{
		StsPatchGrid nextGrid = connection.getPatchGrid();
		StsPatchGrid prevGrid = connection.getPrevPatchGrid();

		if(nextGrid == prevGrid)
			connection.addConnectionToPoint();
		else if(nextGrid == patchGrid)
		{
			PatchPoint prevPoint = connection.getPrevPoint();
			PatchPoint clonedPrevPoint = prevPoint.cloneAndClear();
			patchGrid.addPatchPoint(clonedPrevPoint);
			connection.setPrevPoint(clonedPrevPoint);
			connection.addConnectionToPoint();
		}
		else if(prevGrid == patchGrid)
		{
			Connection clonedConnection = connection.clone();
			clonedConnection.isMoved = true;
			PatchPoint nextPoint = connection.getNextPoint();
			nextPoint.setConnection(clonedConnection);
			PatchPoint clonedNextPoint = nextPoint.cloneAndClear();
			patchGrid.addPatchPoint(clonedNextPoint);
			connection.setNextPoint(clonedNextPoint);
			connection.addConnectionToPoint();
		}
		else
		{
			StsException.systemError(this, "processConnection", "GRID ERROR: Grid " + patchGrid.id +
					" IS NOT nextGrid " + nextGrid.id + " OR " + prevGrid.id);
		}
	}
	private CorrelationWindow findOtherNearestWindow(CorrelationWindow window, ConnectionList connectionList)
	{
		// first guess on window closest is at the same window array index
		// Use this as a start and search up and down for nearest
		int slice = window.getCenterPoint().slice;
		int otherIndex = getBoundedIndex(window.windowIndex);
		int otherIndexAbove = connectionList.getOtherPointIndexAbove();
		int otherIndexBelow = connectionList.getOtherPointIndexBelow();
		otherIndex = StsMath.limitBetweenExclusive(otherIndex, otherIndexAbove, otherIndexBelow);
		CorrelationWindow prevWindow = windows[otherIndex];

		int otherSlice = prevWindow.getCenterPoint().slice;
		byte pointType = StsTraceUtilities.coercedPointTypes[window.getCenterPoint().pointType];

		// find first nearest by searching forward and back
		CorrelationWindow windowBelow, windowAbove;
		Iterator<CorrelationWindow> forwardIterator = new WindowPointTypeForwardIterator(pointType, prevWindow, connectionList);
		windowBelow = forwardIterator.next();
		Iterator<CorrelationWindow> backwardIterator = new WindowPointTypeForwardIterator(pointType, prevWindow, connectionList);
		windowAbove = forwardIterator.next();

		if(windowAbove == null || windowAbove.getCenterPoint().slice < otherSlice)
		{
			if(windowBelow.getCenterPoint().slice <= otherSlice)
				return windowBelow;
			else
			   windowAbove = windowBelow;
		}
		else if(windowBelow == null)
		{
			if(windowAbove.getCenterPoint().slice >= otherSlice)
				return windowAbove;
			else
				windowAbove = windowBelow;
		}
		if (otherSlice < slice)
		{
			forwardIterator = new WindowPointTypeForwardIterator(pointType, prevWindow, connectionList);
			if(!forwardIterator.hasNext()) return null;
			CorrelationWindow nextWindow = forwardIterator.next();
			int nextDistance = nextWindow.getCenterPoint().slice - slice;
			while(forwardIterator.hasNext())
			{
				prevWindow = nextWindow;
				int prevDistance = nextDistance;
				nextWindow = forwardIterator.next();

				nextDistance = nextWindow.getCenterPoint().slice - slice;
				if(prevDistance <= 0 && nextDistance > 0)
				{
					if(-prevDistance < nextDistance)
						return prevWindow;
					else
						return nextWindow;
				}
			}
			return nextWindow;
		}
		else // otherSlice > slice
		{
			Iterator<CorrelationWindow> backwardsIterator = new WindowPointTypeBackwardIterator(pointType, prevWindow, connectionList);
			CorrelationWindow nextWindow = backwardsIterator.next();
			if(nextWindow == null) return null;
			int nextDistance = nextWindow.getCenterPoint().slice - slice;
			while(backwardsIterator.hasNext())
			{
				prevWindow = nextWindow;
				int prevDistance = nextDistance;
				nextWindow = backwardsIterator.next();

				nextDistance = nextWindow.getCenterPoint().slice - slice;
				if(prevDistance >= 0 && nextDistance < 0)
				{
					if(prevDistance > nextDistance)
						return prevWindow;
					else
						return nextWindow;
				}
			}
			return nextWindow;
		}
	}

	CorrelationWindow getWindowOfTypeBelow(CorrelationWindow window, ConnectionList connectionList)
	{
		Iterator<CorrelationWindow> forwardIterator = new WindowPointTypeForwardIterator(window.getCenterPoint().pointType, window, connectionList);
		return forwardIterator.next();
	}

	class WindowPointTypeForwardIterator implements Iterator<CorrelationWindow>
	{
		byte pointType;
		CorrelationWindow window;
		int otherIndexBelow;

		WindowPointTypeForwardIterator(byte pointType, CorrelationWindow window, ConnectionList connectionList)
		{
			this.pointType = pointType;
			this.window = window;
			otherIndexBelow = connectionList.getOtherPointIndexBelow();
			initialize(connectionList);
		}

		void initialize(ConnectionList connectionList)
		{
			if(window.getCenterPoint().pointType == pointType) return;
			window = getWindowBelowOfType(window, pointType, otherIndexBelow);
		}

		public boolean hasNext()
		{
			return window != null;
		}

		public CorrelationWindow next()
		{
			CorrelationWindow currentWindow = window;
			window = getWindowBelowOfType(window, pointType, otherIndexBelow);
			return currentWindow;
		}

		public void remove() {}
	}

	CorrelationWindow getWindowOfTypeAbove(CorrelationWindow window, ConnectionList connectionList)
	{
		Iterator<CorrelationWindow> backwardIterator = new WindowPointTypeBackwardIterator(window.getCenterPoint().pointType, window, connectionList);
		return backwardIterator.next();
	}

	CorrelationWindow getNextWindow(CorrelationWindow window)
	{
		int index = window.windowIndex + 1;
		if(index >= nWindows) return window;
		return windows[index];
	}

	CorrelationWindow getPrevWindow(CorrelationWindow window)
	{
		int index = window.windowIndex - 1;
		if(index < 0) return window;
		return windows[index];
	}

	class WindowPointTypeBackwardIterator implements Iterator<CorrelationWindow>
	{
		byte pointType;
		CorrelationWindow window;
		int otherIndexBelow;

		WindowPointTypeBackwardIterator(byte pointType, CorrelationWindow window, ConnectionList connectionList)
		{
			this.pointType = pointType;
			this.window = window;
			otherIndexBelow = connectionList.getOtherPointIndexBelow();
			if(window.getCenterPoint().pointType == pointType) return;
			this.window = getWindowAboveOfType(window, pointType, otherIndexBelow);
		}

		public boolean hasNext()
		{
			return window != null;
		}

		public CorrelationWindow next()
		{
			CorrelationWindow currentWindow = window;
			window = getWindowAboveOfType(window, pointType, otherIndexBelow);
			return currentWindow;
		}

		public void remove() {}
	}

	CorrelationWindow getWindowBelowOfType(CorrelationWindow window, byte pointType, int otherIndexBelow)
	{
		if(window == null) return null;
		int otherIndex = window.windowIndex;
		for(int i = otherIndex+1; i <= otherIndexBelow-1; i++)
			if(windows[i].getCenterPoint().pointType == pointType) return windows[i];
		return null;
	}

	CorrelationWindow getWindowAboveOfType(CorrelationWindow window, byte pointType, int otherIndexAbove)
	{
		int otherIndex = window.windowIndex;
		for(int i = otherIndex-1; i >= otherIndexAbove+1; i--)
			if(windows[i].getCenterPoint().pointType == pointType) return windows[i];
		return null;
	}

	int getBoundedIndex(int windowIndex)
	{
		return StsMath.minMax(windowIndex, 0, nWindows);
	}

	/** For this new window, we may have a new row and/or col connection or no connection.
	 *  If we have any new connection, then split our bounded connection interval at the new window
	 *  and move the interval down to this new interval. If the window was cloned, we need to use the
	 *  original window (window.clonedPoint) for these operations as it has the trace links for the
	 *  split and move operations.  If no connections, just move the interval down.
	 * @param window
	 * @param colConnection connection from the prevColPoint (same col, prev row) to this new window or its clone
	 * @param rowConnection connection from the prevRowPoint (same row, prev col) to this new window or its clone
	 */
	private void processNewConnections(CorrelationWindow window, Connection colConnection, Connection rowConnection)
	{
		try
		{
			// if we have a new row and/or col connection, split the connection interval
			// by inserting this connection into it
			if (colConnection != null)
			{
				colConnection = addPatchConnection(colConnection, colConnections, false);
				if(colConnection != null) colConnections.insert(colConnection);

			}
			if (rowConnection != null)
			{
				rowConnection = addPatchConnection(rowConnection, rowConnections, true);
				if(rowConnection != null) rowConnections.insert(rowConnection);
			}
			PatchPoint windowCenterPoint = window.getCenterPoint();
			if(windowCenterPoint.getColConnection() != null) colConnections.movePatchInterval(windowCenterPoint.getColConnection());
			if(windowCenterPoint.getRowConnection() != null) rowConnections.movePatchInterval(windowCenterPoint.getRowConnection());
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "processNewConnections", e);
		}
	}
	/*
		private Connection checkAddRowConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
		{
			if (otherTrace == null) return null;
			if (window.centerPoint.rowConnection != null) return null;
			return window.checkAddConnection(otherTrace, minStretchCorrelation, true);
		}

		private Connection checkAddColConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
		{
			if (otherTrace == null) return null;
			if (window.centerPoint.colConnection != null) return null;
			return window.checkAddConnection(otherTrace, minStretchCorrelation, false);
		}

		Connection addConnection(boolean isRow, CorrelationWindow prevWindow, CorrelationWindow window, float stretchCorrelation)
		{
			return new Connection(prevWindow, window, stretchCorrelation);

		}
    */
	private void constructTraceWindows() throws StsException
	{
		windows = new CorrelationWindow[nTracePatchPoints];
		nWindows = 0;
		PatchPoint prevPoint;
		PatchPoint point = null;
		PatchPoint nextPoint = tracePatchPoints[0];
		CorrelationWindow window;

		try
		{
			for (int n = 1; n < nTracePatchPoints - 1; n++)
			{
				prevPoint = point;
				point = nextPoint;
				nextPoint = tracePatchPoints[n];

				window = checkCreateWindow(prevPoint, point, nextPoint, nWindows);
				if (window == null) continue;
				if(nWindows != window.windowIndex)
					StsException.systemDebug(this, "constructTraceWindows", "Index out of sequence.");
				windows[nWindows] = window;
				nWindows++;
			}
			if (nWindows < nTracePatchPoints)
				windows = (CorrelationWindow[]) StsMath.trimArray(windows, nWindows);
			if(nWindows == 0)
				throw new StsException("TracePoints.constructTraceWindows()", "nWindows == 0");
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "constructTraceWindows", e);
			throw new StsException("TracePoints.constructTraceWindows()", e.getMessage());
		}
	}

	private CorrelationWindow checkCreateWindow(PatchPoint prevPoint, PatchPoint point, PatchPoint nextPoint, int windowIndex)
	{
		if(!arePointTypesOK(prevPoint, point, nextPoint)) return null;
		return new CorrelationWindow(prevPoint, point, nextPoint, windowIndex);
	}

	private boolean arePointTypesOK(PatchPoint prevPoint, PatchPoint point, PatchPoint nextPoint)
	{
		if(prevPoint != null)
		{
			if(nextPoint != null)
				return StsTraceUtilities.arePointTypesOK(prevPoint.pointType, point.pointType, nextPoint.pointType);
			else
				return StsTraceUtilities.arePointTypesAboveOK(prevPoint.pointType, point.pointType);
		}
		else // prevPoint == null
		{
			if(nextPoint != null)
				return StsTraceUtilities.arePointTypesBelowOK(point.pointType, nextPoint.pointType);
			else
				return false;
		}
	}

	private void reinitializeTraceIndices(TracePoints prevRowTrace, TracePoints prevColTrace)
	{
		if (this != null) reinitializeTraceIndexing();
		if (prevRowTrace != null) prevRowTrace.reinitializeTraceIndexing();
		if (prevColTrace != null) prevColTrace.reinitializeTraceIndexing();
	}

	void reinitializeTraceIndexing()
	{
		if(rowConnections != null) rowConnections.reinitializeTraceIndexing();
		if(colConnections != null) colConnections.reinitializeTraceIndexing();
	}

	/**
	 * currentPoint is the last currentPoint on this trace in the previous search operation, so is a good starting window for this search
	 * @param slice slice for which we want to find the nearest tracePoint
	 * @return the nearestTracePoint
	 */
	/*
		private PatchPoint nearestPatchPoint(int slice)
		{
			int distance;
			// if currentNearestPoint is above slice, search down for nearest
			if (currentPoint.slice < slice)
			{
				distance = slice - currentPoint.slice;
				PatchPoint point = currentPoint;

				for (int index = currentPoint.traceIndex + 1; index < nTracePatchPoints; index++)
				{
					PatchPoint lastPoint = point;
					int lastDistance = distance;
					point = tracePatchPoints[index];
					// if window is now below slice, then we have bracketed window: set and return currentPoint
					if (point.slice >= slice)
					{
						distance = point.slice - slice;
						if (distance < lastDistance)
							currentPoint = point;
						else
							currentPoint = lastPoint;
						return currentPoint;
					}
					else
						distance = slice - point.slice;
				}
				// didn't bracket, so slice is still below last window; return last window
				currentPoint = point;
			}
			// if currentNearestPoint is below slice, search up for nearest
			else if (currentPoint.slice > slice)
			{
				distance = currentPoint.slice - slice;
				PatchPoint point = currentPoint;

				for (int index = currentPoint.traceIndex - 1; index >= 0; index--)
				{
					PatchPoint lastPoint = point;
					int lastDistance = distance;
					point = tracePatchPoints[index];

					// if window is now above slice, then we have bracketed window: set and return currentPoint
					if (point.slice <= slice)
					{
						distance = slice - point.slice;
						if (distance < lastDistance)
							currentPoint = point;
						else
							currentPoint = lastPoint;
						return currentPoint;
					}
				}
				currentPoint = point;
			}
			return currentPoint;
		}

		private PatchPoint getOtherConnectedPatchPointAbove(boolean isRow)
		{
			return patchPointsList.getConnectionAbove(isRow).otherPoint;
		}

		private PatchPoint getOtherConnectedPatchPointBelow(boolean isRow)
		{
			return patchPointsList.getConnectionBelow(isRow).otherPoint;
		}
	*/
}

class PatchPoint implements Comparable<PatchPoint>, Cloneable
{
	/** trace this point is on */
	TracePoints tracePoints;
	float value;
	float z = StsParameters.nullValue;
	byte pointType;
	StsPatchGrid patchGrid;
	int slice;

	PatchPoint next, prev;
	private Connection rowConnection;
	private Connection colConnection;
	/** index of this window in the trace containing it */
	int traceIndex;
	/** correl factor between this window and next in row. Note that rowConnection is from this window back. */
	// float rowCorrel;
	/** correl factor between this window and next in col.  Note that colConnection is from this window back. */
	// float colCorrel;
	/** cloned window for debugging.  Point this window was cloned from if cloned. */
	PatchPoint clonedPoint;

	/** first window above which has a connected patch */
//		PatchPoint connectionAbove = null;

	/** first window below which has a connected patch */
//		PatchPoint connectionBelow = null;

	PatchPoint()
	{
	}

	/** constructor for first and last links in doubly-linked list of PatchPoints */
	// PatchPoint(int traceIndex) { this.traceIndex = traceIndex; }

	/** constructor for first and last links in doubly-linked ConnectionList */
	PatchPoint(TracePoints trace, int slice)
	{
		this.tracePoints = trace;
		this.slice = slice;
		traceIndex = slice;
	}

	PatchPoint(TracePoints tracePoints, int slice, float z, float value, byte pointType, int traceIndex)
	{
		this.tracePoints = tracePoints;
		this.slice = slice;
		this.z = z;
		this.value = value;
		this.pointType = pointType;
		this.traceIndex = traceIndex;
	}

	PatchPoint(PatchPoint patchPoint, byte pointType, int traceIndex)
	{
		tracePoints = patchPoint.tracePoints;
		slice = patchPoint.slice;
		z = patchPoint.z;
		value = patchPoint.value;
		this.pointType = pointType;
		this.traceIndex = traceIndex;
	}


	PatchPoint(TracePoints tracePoints, int slice, int traceIndex)
	{
		this.tracePoints = tracePoints;
		this.slice = slice;
		this.traceIndex = traceIndex;
	}

	public boolean hasConnection()
	{
		return getRowConnection() != null || getColConnection() != null;
	}

	public boolean hasConnection(boolean isRow)
	{
		if (isRow)
			return getRowConnection() != null;
		else
			return getColConnection() != null;
	}

	public boolean hasRowConnection()
	{
		return getRowConnection() != null;
	}

	public boolean hasColConnection()
	{
		return getColConnection() != null;
	}

	public Connection getConnection(boolean isRow)
	{
		if (isRow) return getRowConnection();
		else return getColConnection();
	}

	PatchPoint getPrevConnectedColPoint()
	{
		if(getRowConnection() == null) return null;
		return getRowConnection().prevWindow.getCenterPoint();
	}

	PatchPoint getPrevConnectedRowPoint(Connection newConnection)
	{
		if(getColConnection() == null) return null;
		return getColConnection().prevWindow.getCenterPoint();
	}

	public int compareTo(PatchPoint otherPoint)
	{
		if (slice > otherPoint.slice) return 1;
		else if (slice < otherPoint.slice) return -1;
		else return 0;
	}

	public boolean isSameRowCol(PatchPoint point)
	{
		if(point != null)
			return getRow() == point.getRow() && getCol() == point.getCol();
		StsException.systemError(this, "isSameRowCol", "POINT IS NULL!");
		return false;
	}

	protected PatchPoint clone()
	{
		try
		{
			PatchPoint clonedPoint = (PatchPoint) super.clone();
			clonedPoint.clonedPoint = this;
			return clonedPoint;
		}
		catch (Exception e)
		{
			StsException.systemError(this, "clone");
			return null;
		}
	}

	/**
	 * A new point needs to be connected to a grid which already has a point at this location.
	 * So clone the window and clear any connection data.  This point will be added to the otherGrid
	 * defined by the otherPoint or to a new grid if there is no grid associated with the otherPoint.
	 * @return the cloned window
	 */
	protected PatchPoint cloneAndClear()
	{
		try
		{
			PatchPoint clonedPoint = clone();
			clonedPoint.clearConnectionData();
			return clonedPoint;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "cloneAndClear", e);
			return null;
		}
	}

	boolean isCloneFromGrid(StsPatchGrid patchGrid)
	{
		if(clonedPoint == null) return false;
		return clonedPoint.patchGrid == patchGrid;
	}

	void clearConnectionData()
	{
		deleteRowConnection();
		deleteColConnection();
		patchGrid = null;
	}

	void setConnection(Connection connection)
	{
		if(StsPatchVolume.debug)
		{
			PatchPoint point = connection.window.getCenterPoint();
			PatchPoint prevPoint = connection.prevWindow.getCenterPoint();
			if(StsPatchGrid.doDebugPoint(point) || StsPatchGrid.doDebugPoint(prevPoint))
				StsException.systemDebug(this, "setConnection", StsPatchVolume.iterLabel + " ADJUST CONNECT point: " +
					point.toString() + " TO point: " + prevPoint.toString());
		}

		if(connection.isRow)
			setRowConnection(connection);
		else
			setColConnection(connection);
	}

	void deleteConnection(boolean isRow)
	{
		if(isRow)
			setRowConnection(null);
		else
			setColConnection(null);
	}

	void deleteConnection(Connection connection)
	{
		if(StsPatchVolume.debug)
		{
			PatchPoint point = connection.window.getCenterPoint();
			PatchPoint prevPoint = connection.prevWindow.getCenterPoint();
			if(StsPatchGrid.doDebugPoint(point) || StsPatchGrid.doDebugPoint(prevPoint))
				StsException.systemDebug(this, "deleteConnection", StsPatchVolume.iterLabel + " DELETE CONNECTION: " +
						connection.toString());
		}

		if(getRowConnection() == connection)
			setRowConnection(null);
		else if(getColConnection() == connection)
			setColConnection(null);
		else
			StsException.systemError(this, "deleteConnection", "CONNECTION DOESN'T EXIST for point: " + toString() +
					" CONNECTION: " + connection.toString());
	}

	void setRowConnectionMoved()
	{
		rowConnection = rowConnection.clone();
		rowConnection.isMoved = true;
	}

	void setColConnectionMoved()
	{
		colConnection = colConnection.clone();
		colConnection.isMoved = true;
	}

	Integer hashCode(int nVolumeCols)
	{
		return new Integer(getRow() * nVolumeCols + getCol());
	}

	int getSlice()
	{
		return slice;
	}

	int getID()
	{
		if (patchGrid == null) return -1;
		else return patchGrid.id;
	}

	int getIndex(int nVolumeCols)
	{
		return getCol() + getRow() * nVolumeCols;
	}

	public String toString()
	{
		if (clonedPoint != null)
			return pointToString() + " cloned from " + clonedPoint.patchToString();
		else
			return pointToString();
	}

	private String patchToString()
	{
		int id = -1;
		if (patchGrid != null) id = patchGrid.id;
		return "id " + id + " ";
	}

	private String pointToString()
	{
		String connectString = (getRowConnection() != null ? " RC" : "") + (getColConnection() != null ? " CC" : "") + " ";
		return patchToString() + "r " + getRow() + " c " + getCol() + " s " + slice + " v " + value +
				" i " + traceIndex + " z " + z + " t " + connectString + StsTraceUtilities.typeStrings[pointType];
	}

	static public String staticToString(PatchPoint point)
	{
		if(point == null) return "NULL";
		else return point.toString();
	}

	String nullOrToString(String string, PatchPoint patchPoint)
	{
		if (patchPoint == null) return " " + string + " null";
		else return " " + string + " " + patchPoint.toString();
	}

	StsPatchGrid getPatchGrid()
	{
		return patchGrid;
	}

	public void setPatchGrid(StsPatchGrid patchGrid)
	{
		this.patchGrid = patchGrid;
	}

	public byte getPointType(boolean useFalseTypes)
	{
		if (!useFalseTypes) return pointType;
		return StsTraceUtilities.coercedPointTypes[pointType];
	}

	public PatchPoint resetIndex(int index)
	{
		traceIndex = index;
		return this;
	}

	final protected int getRow() { return tracePoints.row; }
	final protected int getCol()
	{
		return tracePoints.col;
	}

	/** connection from this tracePoint to the tracePoint on the adjacent trace at this row, col-1 (same row) */
	public Connection getRowConnection()
	{
		return rowConnection;
	}

	public void setRowConnection(Connection rowConnection)
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(this))
				StsException.systemDebug(this, "setRowConnection", StsPatchVolume.iterLabel + " ROW CONNECTION SET FOR Point: " +
						toString() + " rowConnection: " + rowConnection.toString());
		}
		this.rowConnection = rowConnection;
	}

	/** connection from this tracePoint to the tracePoint on the adjacent trace at this row-1, col (same col) */
	public Connection getColConnection()
	{
		return colConnection;
	}

	public void setColConnection(Connection colConnection)
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(this))
				StsException.systemDebug(this, "setColConnection", StsPatchVolume.iterLabel + " COL CONNECTION SET FOR Point: " +
						toString() + " colConnection: " + colConnection.toString());
		}
		this.colConnection = colConnection;
	}

	public void deleteRowConnection()
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(this))
				StsException.systemDebug(this, "deleteRowConnection", StsPatchVolume.iterLabel +
						" ROW CONNECTION DELETED FOR Point: " + toString());
		}
		rowConnection = null;
	}

	public void deleteColConnection()
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(this))
				StsException.systemDebug(this, "deleteColConnection", StsPatchVolume.iterLabel +
						" COL CONNECTION DELETED FOR Point: " + toString());
		}
		colConnection = null;
	}

}

/** A CorrelationWindow has a pointCenter and is bounded by trace points above and below of the appropriate point types.
 *  The window is essentially a half-wave.  Windows have the type of the center point.  A MAX window for example has
 *  a ZP (zero-plus) point above an a ZM (zero-minus) point below.  The skewness of the half-wave is considered in matching
 *  it to other windows so we retain for the window the minus and plus half-widths defined by the slice value difference.
 *  When matched with another windows, the stretchCorrelation is computed as the average of how much each side of the haf-wave
 *  has to be stretched to match the other.  A window with half-widths of -4 and +2 when matched with one with half-widths of
 *  -3 and +4 would have stretch ratios of .75 and 0.5 with an average of .667.  The instantaneous amplitude has been computed
 *  for each window and is used in computing the amplitude ratio (always <= 1.0) between the two windows.  The correlation is
 *  acceptable if this ratio is above a minimum.  The amplitudeRatio is saved in the CorrelationWindow.
 */
class CorrelationWindow implements Cloneable
{
	private PatchPoint centerPoint;
	/** windowType: BELOW if no window above center, ABOVE if no window below or CENTER */
	byte windowType;
	/** slice difference from center window to top window */
	int dSliceMinus;
	/** slice difference from bot window to center window */
	int dSlicePlus;
	/** stretchCorrelation between this window and the connected window */
	float stretchCorrelation;
	/** amplitudeRatio between this window and the connected window */
	float amplitudeRatio;
	/** index of this window in the windows array */
	int windowIndex;
	/** instantaneous amplitude for this window; average of min and max values of window (1 or 2 points) */
	float amplitude;

	public static final byte CENTER = 0;
	public static final byte ABOVE = -1;
	public static final byte BELOW = 1;

	CorrelationWindow(PatchPoint prevPoint, PatchPoint centerPoint, PatchPoint nextPoint, int windowIndex)
	{
		this.setCenterPoint(centerPoint);
		this.windowIndex = windowIndex;
		float dValueMinus = 0.0f, dValuePlus = 0.0f;
		if(prevPoint != null)
		{
			dSliceMinus = centerPoint.slice - prevPoint.slice;
			dValueMinus = Math.abs(centerPoint.value - prevPoint.value);
			if (nextPoint != null)
			{
				dSlicePlus = nextPoint.slice - centerPoint.slice;
				dValuePlus = Math.abs(centerPoint.value - nextPoint.value);
				windowType = CENTER;
			}
			else // nextPoint == null && prevPoint != null
			{
				dSlicePlus = dSliceMinus;
				dValuePlus = dValueMinus;
				windowType = ABOVE;
			}
		}
		else if(nextPoint != null)
		{
			dSlicePlus = nextPoint.slice - centerPoint.slice;
			dValuePlus = Math.abs(centerPoint.value - nextPoint.value);
			dSliceMinus = dSlicePlus;
			dValueMinus = dValuePlus;
			windowType = BELOW;
		}
		else // error both prev and next points are null: shouldn't happen: print message, don't throw exception
		{
			StsException.systemError(this, "constructor", "prev and next points are both null!");
		}
		amplitude = (dValueMinus + dValuePlus)/2;
	}

	public CorrelationWindow clone()
	{
		try
		{
			CorrelationWindow window = (CorrelationWindow) super.clone();
			// window.centerPoint = centerPoint.clone();
			return window;
		}
		catch (Exception e)
		{
			StsException.systemError(this, "clone");
			return null;
		}
	}

	boolean hasRowConnection() { return getCenterPoint().getRowConnection() != null; }
	boolean hasColConnection() { return getCenterPoint().getColConnection() != null; }

	Connection getRowConnection()  { return getCenterPoint().getRowConnection(); }
	Connection getColConnection()  { return getCenterPoint().getColConnection(); }

	boolean isAboveOrEqual(CorrelationWindow prevWindow)
	{
		return getCenterPoint().slice <= prevWindow.getCenterPoint().slice;
	}

	boolean isBelowOrEqual(CorrelationWindow prevWindow)
	{
		return getCenterPoint().slice >= prevWindow.getCenterPoint().slice;
	}

	/** check for closest of two windows where one must be above or equal to and the other must be below or equal to this window.
	 *  two windows can't be at same slice nor the order switched (above is below, below is above) as debug sanity checks
	 *
	 * @param windowAbove window above or equal in slice value to this window
	 * @param windowBelow window below or equal in slice value to this window
	 * @return
	 */
	CorrelationWindow getClosestWindow(CorrelationWindow windowAbove, CorrelationWindow windowBelow)
	{
		if(windowAbove == null)
			return windowBelow;
		else if(windowBelow == null)
			return windowAbove;
		int difAbove = getCenterPoint().slice - windowAbove.getCenterPoint().slice;
		int difBelow = windowBelow.getCenterPoint().slice - getCenterPoint().slice;
		if(StsPatchVolume.debug && (difAbove < 0 || difBelow < 0))
		{
			StsException.systemDebug(this, "getClosestWindow", "window not between windows above and below");
			return null;
		}
		if(StsPatchVolume.debug && (difAbove == 0 && difBelow == 0))
		{
//			StsException.systemDebug(this, "getClosestWindow", "other windows are the same, so can't be between");
			return null;
		}
		if(difAbove <= difBelow) return windowAbove;
		else					 return windowBelow;
	}

	void computeCorrelation(CorrelationWindow prevWindow, float correlPenalty)
	{
		amplitudeRatio = amplitude/prevWindow.amplitude;
		if(amplitudeRatio > 1.0f) amplitudeRatio = 1.0f/amplitudeRatio;

		// check stretchCorrelation stretch
		if (windowType == BELOW)
			stretchCorrelation = computePlusStretchFactor(prevWindow);
		else if (windowType == ABOVE)
			stretchCorrelation = computeMinusStretchFactor(prevWindow);
		else
			stretchCorrelation = (computePlusStretchFactor(prevWindow) + computeMinusStretchFactor(prevWindow)) / 2;
	}
	/** Correlation between two windows is OK if the amplitudeRatio is above a minimum (minAmplitudeRatio*correlation where
	 *  correlation is the current iteration correlation value) and
	 *  the stretchCorrelation is greater than the current correlation including penalties (false min or max).
	 *  store the stretchCorrelation value in the prevWindow
	 * @param correlation
	 * @param correlPenalty
	 * @param minAmplitudeRatio
	 * @return true if amplitudeRatio and stretchCorrelation are >= correlation
	 */
	boolean correlationOK(float correlation, float correlPenalty, float minAmplitudeRatio)
	{
		if(amplitudeRatio < minAmplitudeRatio*correlation)
			return false;

		float totalPenalty = 0.0f;
		if(correlPenalty > 0.0f)
		{
			if(windowType != CENTER) totalPenalty = correlPenalty;
			if(StsTraceUtilities.isPointTypeFalse(getCenterPoint().pointType) || StsTraceUtilities.isPointTypeFalse(this.getCenterPoint().pointType))
				totalPenalty += correlPenalty;
		}
		stretchCorrelation -= totalPenalty;
		return stretchCorrelation >= correlation;
	}

		/*
			private Connection checkAddConnection(TracePoints otherTrace, float minStretchCorrelation, boolean isRow)
			{
				CorrelationWindow otherMatchingWindow = findOtherMatchingWindow(otherTrace, isRow, minStretchCorrelation);
				if (otherMatchingWindow == null) return null;
				return addPatchConnection(otherMatchingWindow, isRow);
			}

			private ConnectionList getConnectionList(boolean isRow)
			{
				if(isRow) return rowConnections;
				else	  return colConnections;
			}

			boolean pointTypesMatch(CorrelationWindow prevWindow)
			{
				byte otherCenterType = prevWindow.centerPointType;
				if (centerPointType == otherCenterType) return true;
				if (!useFalseTypes) return false;

				byte centerType = StsTraceUtilities.coercedPointTypes[centerPointType];
				otherCenterType = StsTraceUtilities.coercedPointTypes[otherCenterType];
				return centerType == otherCenterType;
			}
        */
	/**
	 * check the above and below types to see that they match.
	 * We are assuming the centers have already been checked for matches
	 * @param window the window
	 * @param prevWindow prevWindow we are comparing it to
	 * @return true if centerTypes, and above and below types match
	 */
		/*
			boolean windowTypesMatch(CorrelationWindow window, CorrelationWindow prevWindow)
			{
				byte above = window.pointAbove.getPointType();
				byte below = window.pointBelow.getPointType();
				byte otherAbove = prevWindow.pointAbove.getPointType();
				byte otherBelow = prevWindow.pointBelow.getPointType();
				return above == otherAbove && below == otherBelow;
			}
        */
	/** if we have two windows with the exactly identical centerPoint, they must be equivalent if not equal windows. */
	boolean sameAs(CorrelationWindow prevWindow)
	{
		return prevWindow.getCenterPoint() == getCenterPoint();
	}

	float computeMinusStretchFactor(CorrelationWindow prevWindow)
	{
		float minusStretchFactor = ((float) dSliceMinus) / prevWindow.dSliceMinus;
		if (minusStretchFactor > 1.0f)
			minusStretchFactor = 1 / minusStretchFactor;
		return minusStretchFactor;
	}

	float computePlusStretchFactor(CorrelationWindow prevWindow)
	{
		float plusStretchFactor = ((float) dSlicePlus) / prevWindow.dSlicePlus;
		if (plusStretchFactor > 1.0f)
			plusStretchFactor = 1 / plusStretchFactor;
		return plusStretchFactor;
	}
	/*
		boolean isCenterSliceOutsideWindow(int centerSlice)
		{
			return centerSlice < minSlice || centerSlice > maxSlice;
		}
	*/
	public String toString()
	{
		return " index " + windowIndex + " centerPoint: " + PatchPoint.staticToString(centerPoint); // + " stretchCorrelation: " + stretchCorrelation + " amplitudeRatio: " + amplitudeRatio;
	}

	/** window in center of this window */
	public PatchPoint getCenterPoint()
	{
		if(centerPoint == null)
			StsException.systemError(this, "getCenterPoint", "CENTER POINT IS NULL!");
		return centerPoint;
	}

	public void setCenterPoint(PatchPoint centerPoint)
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(centerPoint))
				StsException.systemDebug(this, "setCenterPoint", StsPatchVolume.iterLabel + " CENTER POINT SET: " +
						toString() + " TO point: " + centerPoint.toString());
		}
		if(centerPoint == null)
			StsException.systemError(this, "setCenterPoint", "NEW CERNTER POINT IS NULL!");
		this.centerPoint = centerPoint;
	}

	void deleteRowConnection()
	{
	 	centerPoint.setRowConnection(null);
	}

	void deleteColConnection()
	{
		centerPoint.setColConnection(null);
	}
	/*
	private float computeStretchCorrelation(CorrelationWindow prevWindow)
	{
		if (prevWindow == null) return 0.0f;

		TracePoints traceOther = prevWindow.getTracePoints();
		int centerOther = prevWindow.centerSlice;
		int minOther = prevWindow.minSlice;
		int maxOther = prevWindow.maxSlice;

		// translate and stretch/shrink pointsOther z values so they line up with pointsNew

		int dSliceMinusOther = centerOther - minOther;
		int dSliceMinusNew = centerSlice - minSlice;
		float dSliceMinusOtherScalar = (float) dSliceMinusNew / dSliceMinusOther;
		// if(dzMinusOtherScalar < minStretchLimit || dzMinusOtherScalar > maxStretchLimit) return 0.0f;

		float minusStretchFactor = dSliceMinusOtherScalar;
		if (minusStretchFactor > 1.0f)
			minusStretchFactor = 1 / minusStretchFactor;

		int dSlicePlusOther = maxOther - centerOther;
		int dSlicePlusNew = maxSlice - centerSlice;
		float dSlicePlusOtherScalar = (float) dSlicePlusNew / dSlicePlusOther;
		// if(dzPlusOtherScalar < minStretchLimit || dzPlusOtherScalar > maxStretchLimit) return 0.0f;

		float plusStretchFactor = dSlicePlusOtherScalar;
		if (plusStretchFactor > 1.0f)
			plusStretchFactor = 1 / plusStretchFactor;

		return Math.min(minusStretchFactor, plusStretchFactor);
	}
	*/
}
/**
 * tracePoints has a series of row connections and col connections which are maintained during construction of this trace.
 * They will be removed after construction.  When checking on a new connection between this trace and the otherTrace
 * which is either row or col aligned with this trace, we bracket the search by connections above and below to prevent
 * crossing them.  On completion, if a new connection is created, it is added to to either trace.rowConnections or trace.colConnections.
 * These connection lists are currently double-linked, but could perhaps be only single-linked.
 */
class Connection implements Cloneable
{
	Connection next, prev;
	/** connected window on this trace */
	CorrelationWindow window;
	/** connected window on other trace */
	CorrelationWindow prevWindow;
	/** connection is in same row */
	boolean isRow;
	/** average of slice values for two connected points; since connections can't cross or be identical, it will provide correct order */
	float sliceAvg;
	/** stretch correlation between these two points */
	//float stretchCorrelation;
	/** amplitude ratio between these two Points  */
	//float amplitudeRatio;
	/** isMoved indicates this connection exists but has been moved to other patchGrid; a clone is left here indicating that
	 *  the connection exists and doesn't need to be regenerated on a subsequent iteration; because it isMoved, however, an
	 *  actual connection doesn't exist; the correlation between the two connected points is zero. */
	boolean isMoved = false;
  	/** construct connection from window to prevWindow. Clone the window since it may be shared by two connections
	 *  and one may affect the other. The window and cloned window share the same centerPoint which has the row and col connections.
	 * @param prevWindow
	 * @param window
	 */
	Connection(CorrelationWindow prevWindow, CorrelationWindow window)
	{
		this.window = window.clone();
		this.prevWindow = prevWindow;
		sliceAvg = (this.window.getCenterPoint().slice + prevWindow.getCenterPoint().slice)/2.0f;
		//this.stretchCorrelation = window.stretchCorrelation;
		//this.amplitudeRatio = window.amplitudeRatio;
		isRow = prevWindow.getCenterPoint().getRow() == this.window.getCenterPoint().getRow();
	}

	PatchPoint getWindowCenterPoint(TracePoints tracePoints)
	{
		if(window.getCenterPoint().tracePoints == tracePoints) return window.getCenterPoint();
		else if(prevWindow.getCenterPoint().tracePoints == tracePoints) return prevWindow.getCenterPoint();
		StsException.systemDebug(this, "getWindowCenterPoint.  tracePoints ", tracePoints.toString() + " not found in " + toString());
		return null;
	}

	public PatchPoint getNextPoint()
	{
		return window.getCenterPoint();
	}

	public PatchPoint getPrevPoint()
	{
		return prevWindow.getCenterPoint();
	}

	public void setNextPoint(PatchPoint patchPoint)
	{
		window.setCenterPoint(patchPoint);
	}

	public void setPrevPoint(PatchPoint patchPoint)
	{
		prevWindow.setCenterPoint(patchPoint);
	}

	StsPatchGrid getPatchGrid()
	{
		return window.getCenterPoint().patchGrid;
	}

	StsPatchGrid getPrevPatchGrid()
	{
		return prevWindow.getCenterPoint().patchGrid;
	}

	public void movePoint(PatchPoint newPatchPoint)
	{
		PatchPoint oldPoint = window.getCenterPoint();
		oldPoint.deleteConnection(this);
		newPatchPoint.setConnection(this);
	}

	public final void addConnectionToPoint()
	{
		PatchPoint patchPoint = window.getCenterPoint();
		if (isRow)
		{
			patchPoint.setRowConnection(this);
			// otherPatchPoint.rowCorrel = connection.stretchCorrelation;
		}
		else
		{
			patchPoint.setColConnection(this);
			// otherPatchPoint.colCorrel = connection.stretchCorrelation;
		}
	}

	public void resetConnections(PatchPoint fromPoint, PatchPoint toPoint)
	{
		resetConnection(fromPoint);
		resetConnection(toPoint);
	}
	/** We have an existing connection which we need to reset as the connected points have changed (cloned or moved to another grid).
	 *  The connection is from window.centerPoint to prevWindow.centerPoint and the connection is specified only as a member
	 *  (rowConnection or colConnection of window.centerPoint.
	 *
	 * @param patchPoint the new point to be set as connection start at window.centerPoint.
	 */
	public void resetConnection(PatchPoint patchPoint)
	{
		if(StsPatchVolume.debug)
		{
			if(StsPatchGrid.doDebugPoint(patchPoint))
				StsException.systemDebug(this, "resetConnection", StsPatchVolume.iterLabel + " RESET CONNECTION POINT: " +
						toString() + " TO point: " + patchPoint.toString());
		}

		if(patchPoint.getID() == -1)
		{
			StsException.systemError(this, "resetConnection", "POINT HAS NO PATCH ID: " + patchPoint.toString());
			return;
		}

		if(window.getCenterPoint().isSameRowCol(patchPoint))
		{
			window.setCenterPoint(patchPoint);
			patchPoint.setConnection(this);
		}
		else if(prevWindow.getCenterPoint().isSameRowCol(patchPoint))
		{
			prevWindow.setCenterPoint(patchPoint);
		}
		else // shouldn't occur
			StsException.systemError(this, "resetConnection", "POINT NOT FOUND FOR CONNECTION: " + toString() + " POINT: " + patchPoint.toString());
	}

	public void restorePoints(PatchPoint nextPoint, PatchPoint prevPoint)
	{
		if(window.getCenterPoint() == nextPoint && prevWindow.getCenterPoint() == prevPoint) return;
		StsException.systemError(this, "restorePoints", "HAD TO RESTORE POINTS: CHECK!");
		window.setCenterPoint(nextPoint);
		prevWindow.setCenterPoint(prevPoint);
	}

	public void moveOtherPoint(PatchPoint otherPatchPoint)
	{
		prevWindow.setCenterPoint(otherPatchPoint);
	}
	/** Given a newPoint replacing an oldPoint, replace the oldPoint with new on this connection.
	 *  We don't know which end of the connection it is on, so check both and clonedPoint of both.
	 * @param newPoint
	 * @param oldPoint
	 */
	public void checkResetClonedPoint(PatchPoint newPoint, PatchPoint oldPoint)
	{
		if(window.getCenterPoint() == oldPoint)
		{
			window.setCenterPoint(newPoint);
			newPoint.setConnection(this);
			oldPoint.deleteConnection(this);
		}
		else if(window.getCenterPoint().clonedPoint == oldPoint)
		{
			window.setCenterPoint(newPoint);
			newPoint.setConnection(this);
			oldPoint.deleteConnection(this);
		}
		else if(prevWindow.getCenterPoint() == oldPoint)
			prevWindow.setCenterPoint(newPoint);
		else if(prevWindow.getCenterPoint().clonedPoint == oldPoint)
			prevWindow.setCenterPoint(newPoint);
		else
			StsException.systemError(this, "checkResetClonedPoint", "POINT NOT FOUND ON Connection: " + toString() +
					"OLD POINT: " + oldPoint.toString());
	}

	public void checkResetClonedPoints(PatchPoint otherPatchPoint, PatchPoint newPatchPoint)
	{
		if (prevWindow.getCenterPoint() != otherPatchPoint) // otherPatchPoint may have been cloned
			prevWindow.setCenterPoint(otherPatchPoint);
		else if (window.getCenterPoint() != newPatchPoint) // newPatchPoint may have been cloned
		{
			window.getCenterPoint().deleteConnection(this);
			newPatchPoint.setConnection(this);
		}
		else
			StsException.systemError(this, "checkResetClonedPoints", "POINT NOT FOUND ON Connection: " + toString() +
					" otherPoint: " + otherPatchPoint.toString() + " newPatchPoint: " + newPatchPoint.toString());
	}

	public Connection clone()
	{
		try
		{
			Connection clonedConnection = (Connection)super.clone();
			clonedConnection.window = (CorrelationWindow)window.clone();
			clonedConnection.prevWindow = (CorrelationWindow)prevWindow.clone();
			return clonedConnection;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "clone", e);
			return null;
		}
	}

	public String toString()
	{
		String rowColString = (isRow) ? "ROW" : "COL";
		return rowColString + " CONNECTION window: " + window.toString() + "     TO OTHER window: " + prevWindow.toString() +
				" sliceAvg: " + sliceAvg; //  + " correl: " + stretchCorrelation + " amplitudeRatio: " + amplitudeRatio;
	}

	static public String staticToString(Connection connection)
	{
		if(connection == null)
			return "CONNECTION IS NULL";
		return connection.toString();
	}
}
/**
 * doubly linked list of Connection[s].  There are two lists for each trace: rowConnections to the prev col on same row,
 * and colConnections to the prev row on the same col.  Connections in the list are in order, and must not cross or be identical.
 * Order is determined by the avg of the two slice values of the connected points (@see Connection).
 */
class ConnectionList
{
	/** first window in link list (connected to first actual window in list) */
	final Connection first;
	/** last window in link list (connected to last actual window in list) */
	final Connection last;
	/** last connected window in linked list just above current window */
	Connection connectionAbove;
	/** connected window just below connectionAbove in linked list */
	Connection connectionBelow;
	/** last connected window that was set; a convenient starting window for any search */
	// Connection currentConnection;

	/** Insert inactive row and col connections at the top and bottom of the connectionLists.
	 * New connections are added to these lists in order of the sliceAvg of the connection.
 	 * @param trace connection is from trace back to the otherTrace
	 * @param prevTrace trace connected which is either prevRow or prevCol
	 */
	ConnectionList(TracePoints trace, TracePoints prevTrace)
	{

		CorrelationWindow firstWindow = trace.windows[0].clone();
		firstWindow.windowIndex -= 1;
		firstWindow.getCenterPoint().slice -= 1;

		CorrelationWindow lastWindow = trace.windows[trace.nWindows-1].clone();
		lastWindow.windowIndex = trace.nWindows;
		lastWindow.getCenterPoint().slice += 1;

		CorrelationWindow firstPrevWindow = prevTrace.windows[0].clone();
		firstPrevWindow.windowIndex = -1;
		firstPrevWindow.getCenterPoint().slice -= 1;

		CorrelationWindow lastPrevWindow = prevTrace.windows[prevTrace.nWindows-1].clone();
		lastPrevWindow.windowIndex = prevTrace.nWindows;
		lastPrevWindow.getCenterPoint().slice += 1;

		first = new Connection(firstPrevWindow, firstWindow);
		last = new Connection(lastPrevWindow, lastWindow);
		first.next = last;
		last.prev = first;
		// currentConnection = first;
		connectionAbove = first;
		connectionBelow = last;
		connectionAbove.next = last;
		connectionBelow.prev = first;
	}

	void reinitializeTraceIndexing()
	{
		connectionAbove = first;
		connectionBelow = first.next;
		// currentConnection = first;
	}

	int getOtherPointIndexAbove()
	{
		return connectionAbove.prevWindow.windowIndex;
	}

	int getOtherPointIndexBelow()
	{
		return connectionBelow.prevWindow.windowIndex;
	}

	/** we have moved down to a new existing correlated patchPoint; set the interval to the one between this patchPoint and the window below */
	void movePatchInterval(Connection connectionAbove)
	{
		this.connectionAbove = connectionAbove;
		if(connectionAbove.next != null)
			connectionBelow = connectionAbove.next;
	}

	/**
	 * we are inserting this connectedPoint in an interval between connectionAbove and connectionBelow, either of which could be null
	 * meaning that it could be an openInterval with above and/or below undefined.  The interval (open or closed) is
	 * split into two subintervals and the current interval is set to the lower subinterval.
	 * @param connection between pointAbove and pointBelow where interval is to be split into two subIntervals.
	 */
	void insert(Connection connection)
	{
		if ( StsPatchVolume.debug && (connection == connectionAbove || connection == connectionBelow))
		{
			StsException.systemDebug(this, "insert", " connection " + connection.toString() + " same as " + connectionAbove.toString() + " or " + connectionBelow.toString());
			return;
		}
		connectionAbove.next = connection;
		connection.prev = connectionAbove;
		connection.next = connectionBelow;
		connectionBelow.prev = connection;

		// connectionAbove = connection;
		// connectionBelow = connection.next;
	}

	boolean connectionsCross(Connection connection)
	{
		if (!connectionsCross(connection, connectionAbove) && !connectionsCross(connection, connectionBelow))
			return false;

		if(StsPatchVolume.debug && StsPatchGrid.debugPoint && (StsPatchGrid.doDebugPoint(connection.window.getCenterPoint()) || StsPatchGrid.doDebugPoint(connection.prevWindow.getCenterPoint())))
			StsException.systemDebug(this, "connectionCrosses", StsPatchVolume.iterLabel + connection.toString());

		return true;
	}

	boolean connectionsCross(Connection c1, Connection c2)
	{
		int crosses = StsMath.signProduct(c1.window.getCenterPoint().slice - c2.window.getCenterPoint().slice, c1.window.getCenterPoint().slice - c2.window.getCenterPoint().slice);
		return crosses < 0;
	}
}