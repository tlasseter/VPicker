package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

public class StsHorpickWizard extends StsWizard
{
    public StsSelectHorpick selectHorpick;
    public StsDefineHorpick defineHorpick;
    public StsRunHorpick runHorpick;
    public StsFinishHorpick finishHorpick;

    private StsCursorPoint seedPoint = null;
    private StsPickPatch pickPatch = null;
    private int pickPatchIndex = -1;  // index of current pickPatch
    private int currentPatchNumber = -1;  // fixed number assigned to pickPatch
//    private boolean selectSeed = false;
    private StsHorpick selectedHorpick = null;
    /** type of pick to make: see pickTypes (convenience copy from StsPickPatch) */
    public byte pickType = PICK_MAX;
    /** window size for correlation around pick point (convenience copy from StsPickPatch) */
    public int windowSize = defaultWindowSize;
    /** minimum acceptable cross-correlation for pick point (convenience copy from StsHorpick) */
    public float minCorrel;
    /** avg distance between picks from trace to trace on test area */
    public double avgTestPickDif;
    /** max distance between picks from trace to trace on test area */
    public float maxTestPickDif;
    /** number of pickDifs; we will divide the sum (stored in avgTestPickDif) by this number to get avg */
    public int nTestPickDifs;
    /** min wavelength in samples for test area */
    public float minWaveLength;
    /** max wavelength in samples for test area */
    public float maxWaveLength = 0.0f;
    /** colorscale used for corCoefs; horpick gets colorscale from here to create colorDisplayList */
//    public StsColorscalePanel corCoefsColorscalePanel = null;

    /** criteria for merging patches to surface (convenience copy from StsPickPatch) */
//    public byte stopCriteria = STOP_CRITERIA_RERUN;
    /** current propertyName being isVisible (or to be isVisible) on picked surface */
    private String displayPropertyName = StsHorpick.displayPropertyPatchColor;
    /** list of properties which can be isVisible on picked surface */
    private String[] displayPropertyNames = null;

    private StsWizardStep[] mySteps =
    {
        selectHorpick = new StsSelectHorpick(this),
        defineHorpick = new StsDefineHorpick(this),
        runHorpick = new StsRunHorpick(this),
        finishHorpick = new StsFinishHorpick(this)
    };

    /** horizon grid being picked */
    StsSurface surface = null;
	/** list of surfaces that have been picked in this wizard session */
	ArrayList processedSurfaces = new ArrayList(4);
    /** volume being auto-picked */
    StsSeismicVolume seismicVolume;
    /** patch grid z values */
//    float[][] pointsZ = null;
    /** patch grid null types */
//    byte[][] pointsNull = null;
    /** best correlation for this trace and neighboring traces on surface */
//    float[][] surfaceCorCoefs = null;
    /** either patch or surface corCoefs depending on process */
    float[][] corCoefs = null;
    /** row of point we want to interpolate */
    int iCenter;
    /** col of point we want to interpolate */
    int jCenter;
    /** surface grid z values */
    float[][] pointsZ = null;
    /** surface grid null types */
    byte[][] pointsNull = null;
    /** index of patch at each row-col */
    byte[] patchIndexes = null;
    /** print mainDebug information */
    boolean debug = false;
    /** size of bounding box for interpolation */
    int rowMin, rowMax, colMin, colMax;
    /** number of surface grid rows */
    int nRows;
    /** number of surface grid cols */
    int nCols;

    /** number of possible picks on a new trace */
    int nPossiblePicks = 0;
    /** location of possiblePicks defined by float sampleNumber */
    float[] possiblePicks = null;
    /** Array of traces (allocated as needed) */
    Trace[][] traces = null;
    /** When memory-limited, don't save traces in trace array (useTraces = false) */
    boolean useTraces = true;
    /** windowSize plus margins */
    int maxWindowLength = 0;
    /** max offset between adjacent traces for same event: controlled by pickPatch, but this is the latest value which is assigned to any new pickPatch */
    float maxPickDif = defaultMaxPickDif;
    /** total number of samples in a traces (the vertical size of the seismic volume) */
    int nSamples = 0;
    /** number of firstSample in trace window */
    int nFirstSample = 0;
    /** number of lastSample in trace window */
    int nLastSample = 0;
    /** state flag for memory allocation checks */
    byte memoryCheckFlag = MEMORY_CHECK_INITIAL;
    /** number of traces which have been created */
    int nAllocatedTraces = 0;
    /** number of traces created for the current patch; assigned to patch when pick is completed. */
    int nPatchTraces = 0;
    /** maxMemory available */
    long maxMemory = 0;
    /** interval in trace number between memory checks */
    long traceCheckInterval = Long.MAX_VALUE;
    /** applied to all pickPatches if true */
    boolean selectAllPatches = true;
    /** indicates this is a test patch */
    boolean isTestPatch = false;
    /** in range 0-254, this is the data zero crossing */
//    float unsignedByteAverage;

    Thread pickSurfaceThread = null;

    boolean stopThread = false;

//    static boolean saveCenteredValues = false;
    /** adjustedValues for prev trace are saved here instead of in each trace */
    static float[] staticPrevAdjustedValues = null;
    /** adjustedValues for next trace are saved here instead of in each trace */
    static float[] staticNextAdjustedValues = null;
    // mainDebug statistics
//    int[] posStat = new int[100];
//    int[] corStat = new int[100];

    int traceInstanceByteSize;

    int iterationDisplayInterval = 0;

    public static final byte PICK_MAX = StsRunHorpickPanel.PICK_MAX;
    public static final byte PICK_MIN = StsRunHorpickPanel.PICK_MIN;
    public static final byte PICK_ZERO_PLUS = StsRunHorpickPanel.PICK_ZERO_PLUS;
    public static final byte PICK_ZERO_MINUS = StsRunHorpickPanel.PICK_ZERO_MINUS;

    /** before picking a trace, check whether a pick already exists on surface at this i-j
     *  and responding according to the criterion selected. */
    static final byte STOP_CRITERIA_STOP = 0; // don't pick if surface has one
    static final byte STOP_CRITERIA_REPLACE = 1; // replace one on surface
    static final byte STOP_CRITERIA_IF_SAME_Z = 2; // replace unless at same Z within a tolerance
    static final byte STOP_CRITERIA_RERUN = 3; // points added on previous run considered as added on this run

    public static final byte MEMORY_CHECK_INITIAL = 0;
    public static final byte MEMORY_CHECK_NO = -1;
    public static final byte MEMORY_CHECK_YES = 1;

    static int [] N = new int [] { 1, 0 };
    static int [] S = new int [] { -1, 0 };
    static int [] E = new int [] { 0, 1 };
    static int [] W = new int [] { 0, -1 };
    static int [] NE = new int [] { 1, 1 };
    static int [] SE = new int [] { -1, 1 };
    static int [] SW = new int [] { -1, -1 };
    static int [] NW = new int [] { 1, -1 };

    static int[][] forwardBottom = new int[][] { N, W, NW, NE };
    static int[][] reverseBottom = new int[][] { S, E, SE, SW };
    static int[][] forwardRight  = new int[][] { W, S, SW, NW };
    static int[][] reverseRight  = new int[][] { E, N, NE, SE };
    static int[][] forwardTop    = new int[][] { S, E, SE, SW };
    static int[][] reverseTop    = new int[][] { N, W, NW, NE };
    static int[][] forwardLeft   = new int[][] { E, N, NE, SE };
    static int[][] reverseLeft   = new int[][] { W, S, SW, NW };
    static int[][] forwardBL1    = new int[][] { NE };
    static int[][] forwardBL2    = new int[][] { N, E};
    static int[][] reverseBL     = new int[][] { S, W, SE, SW, NW };
    static int[][] forwardBR1    = new int[][] { W, NW };
    static int[][] forwardBR2    = new int[][] { N };
    static int[][] reverseBR     = new int[][] { E, S, NE, SE, SW };
    static int[][] forwardTR1    = new int[][] { S, SW };
    static int[][] forwardTR2    = new int[][] { W };
    static int[][] reverseTR     = new int[][] { N, E, SE, NE, NW };
    static int[][] forwardTL1    = new int[][] { E, SE };
    static int[][] forwardTL2    = new int[][] { S };
    static int[][] reverseTL     = new int[][] { N, W, NW, SW, NE };

    static final float nullValue = StsParameters.nullValue;
    static final float largeFloat = StsParameters.largeFloat;

    static final byte SURF_NULL = StsSurface.SURF_PNT;
    static final byte SURF_GAP_NOT_FILLED = StsSurface.SURF_GAP_NOT_FILLED;
    static final byte SURF_GAP_SET = StsSurface.SURF_GAP_SET;

    static final float minFraction = .8f;
    static final float maxFraction = .9f;

    /** this is the most that one trace might be offset against another */
    static final int maxLag = 2;
    static final int maxPossibleWindow = 100;
    static final int defaultWindowSize = 21;
    static final float defaultMaxPickDif = 1.0f;

    /** scratch values used in trace window calculations */
    static float[] floatValues;

    public StsHorpickWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Select/Define Horpick");
        dialog.getContentPane().setSize(550, 575);
        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        if(pickSurfaceThread != null && pickSurfaceThread.isAlive())
        {
            pickSurfaceThread.stop();
        }
        deleteTraces();
//        if(line2d != null) line2d.deleteSliceBlocks();
        success = true;
        if(surface != null)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
            model.instanceChange(surface, "surfacePick");
        }
		if(processedSurfaces.size() > 0)
		{
			StsSurface processedSurface;
			for(int n = 0; n < processedSurfaces.size(); n++)
			{
				processedSurface = (StsSurface) processedSurfaces.get(n);
//                processedSurface.interpolateNullZPointsRadiusWeighted();
				processedSurface.interpolateDistanceTransform();
				processedSurface.constructGrid();
                processedSurface.setDefaultSurfaceTexture(model.win3d.getGlPanel3d());
			}
		}

        if(selectedHorpick != null)
        {
            this.setHorpickVisible(false);

            if(!selectedHorpick.isComplete())
                selectedHorpick.delete();
            else
            {
                selectedHorpick.cleanup();
                model.instanceChange(selectedHorpick, "horpick");
            }
        }
//        restoreMemoryFractions();
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectHorpick)
        {
            if(selectedHorpick != null) selectedHorpick.cleanup();
//            selectedHorpick = selectHorpick.getSelectedHorpick();
            if(selectedHorpick == null)
            {
                disableFinish();
                new StsMessage(frame, StsMessage.INFO, "No horizon selected.");
                return;
            }
            else
            {
                selectedHorpick.setWizard(this);
                setHorpickVisible(true);
                corCoefs = selectedHorpick.getCorCoefs();
 //               setCorCoefsColorscaleRange();
//                enableFinish();
                model.setCurrentObject(selectedHorpick);
                StsGLPanel3d glPanel3d = model.win3d.getGlPanel3d();
                StsView currentView = glPanel3d.getView();
                if(!(currentView instanceof StsView3d) || !(currentView instanceof StsViewCursor))
                                 glPanel3d.checkAddView(StsView3d.class);
//                model.win3d.getViewSelectToolbar().setButtonEnabled(model.glPanel3d.getCurrentView(), true);
                surface = selectedHorpick.getSurface();
                surface.initializeSurface();
                pointsNull = surface.getPointsNull();
                surface.setSurfaceTexture(selectedHorpick.getPatchColorSurfaceTexture(), glPanel3d);
//                selectedHorpick.displayState = StsHorpick.BUILD;
                seismicVolume = selectedHorpick.getSeismicVolume();
//                line2d.setMemoryFractions(minFraction, maxFraction);
//                enableFinish();
                model.win3dDisplay();
                gotoStep(runHorpick);
            }
        }
        else if(currentStep == defineHorpick)
        {
            selectedHorpick = defineHorpick.getDefinedHorpick();
            if(selectedHorpick == null)
            {
                disableFinish();
                new StsMessage(frame, StsMessage.INFO, "No horizon defined.");
                return;
            }
            model.setCurrentObject(selectedHorpick);
            setHorpickVisible(true);
//            setCorCoefsColorscaleRange();
            selectedHorpick.setWizard(this);
            model.win3dDisplay();
            gotoStep(runHorpick);
        }
        else if(currentStep == runHorpick)
        {
            gotoStep(finishHorpick);
//            finishHorpick.panel.colorscaleEnabled(false);
            enableFinish();
        }
    }

    private void setHorpickVisible(boolean isVisible)
    {
        StsHorpickClass horpickClass = (StsHorpickClass)model.getStsClass(StsHorpick.class);
        horpickClass.setIsVisible(isVisible);
    }

	public void setSelectedHorpick(StsHorpick horpick) { selectedHorpick = horpick; }
	public StsHorpick getSelectedHorpick()
	{
		return selectedHorpick;
	}


/*
    public void finish()
    {
        next();
        super.finish();
    }
*/
    public float getMaxPickDif()
    {
        if(pickPatch == null) return defaultMaxPickDif;
        return pickPatch.getMaxPickDif();
    }
    public void setMaxPickDif(float maxPickDif)
    {
        this.maxPickDif = maxPickDif;
        if(pickPatch == null) return;
        pickPatch.setMaxPickDif(maxPickDif);
    }

    public byte getPickType() { return pickType; }
    public void setPickType(byte pickType)
    {
        this.pickType = pickType;
        if(pickPatch == null) return;
        pickPatch.setPickType(pickType);
    }

    public StsHorpick getHorpick() { return selectedHorpick; }
    public StsSeismicVolume[] getSeismicVolumes()
    {
        return (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
    }

    public void createNewHorpick()
    {
        gotoStep(defineHorpick);
    }

/*
    public void setSelectSeed(boolean isSelected)
    {
        selectSeed = isSelected;
    }
*/
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
		if(selectedHorpick == null) return true;
//        if(selectedHorpick == null || !selectSeed) return true;
//        if(selectedHorpick == null || !runHorpick.panel.selectButton.isSelected()) return true;

        if( mouse.getCurrentButton() != StsMouse.LEFT ) return true;
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        // set temporary cursorPoint on mouse PRESS or DRAG
        if(buttonState == StsMouse.PRESSED || buttonState == StsMouse.DRAGGED)
        {
            StsView currentView = glPanel3d.getView();
            if(currentView instanceof StsView3d)
            {
                StsCursor3d cursor3d = glPanel3d.getCursor3d();
                seedPoint = cursor3d.getCursorPoint(glPanel3d, mouse);
				if (seedPoint != null)
					cursor3d.setCurrentDirNo(seedPoint.dirNo);
				else
					new StsMessage(model.win3d, StsMessage.WARNING, "Cannot add a pick (no pick)");
            }
            else if(currentView instanceof StsViewCursor)
                seedPoint = ((StsViewCursor)currentView).getCursorPoint(mouse);
            else
                return true;
        }
        // permanently add this point when mouse.RELEASED
        else if(StsMouse.isButtonStateReleasedOrClicked(buttonState))
        {
            surface = selectedHorpick.getSurface();
            if(surface == null)
            {
                StsClass volumeClass = this.model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume");
                seismicVolume  = (StsSeismicVolume)volumeClass.getCurrentObject();
                if(seismicVolume == null)
                {
                    StsException.systemError("StsHorpick constructor failed. no seismic volume found.");
                    return false;
                }
                selectedHorpick.setSeismicVolume(seismicVolume);
//                line2d.setMemoryFractions(minFraction, maxFraction);
                surface = seismicVolume.constructSurface(selectedHorpick.getName(), selectedHorpick.getStsColor(), StsSurface.IMPORTED);
                model.refreshObjectPanel(surface); // adds surface to object panel
//                surface.initializeGrid(StsSurface.NULL_GAP_NOT_FILLED, true);
                selectedHorpick.setSurface(surface);
            }

			if (seedPoint == null) return false; // jbw no pick
			StsGridPoint gridPoint = new StsGridPoint(seedPoint.point, surface);
			byte surfacePointNull = surface.getPointNull(gridPoint.row, gridPoint.col);
			if(surface.isPointNotNull(gridPoint.row, gridPoint.col))
			{
				int patchIndex = selectedHorpick.getPatchIndex(gridPoint.row, gridPoint.col);
				new StsMessage(model.win3d, StsMessage.WARNING, "Cannot add a pick as another pick patch already exists here.\n" +
							   "Delete pick number " + patchIndex + " first and then add this pick.");
				return true;
			}
			byte pickType = runHorpick.panel.getPickType();
			pickPatch = StsPickPatch.constructor(selectedHorpick, pickType, this);
			gridPoint.setNullType(StsGridPoint.SURF_PNT);
			pickPatch.setGridPoint(gridPoint);
			pickPatchIndex = pickPatch.getIndex();
			selectedHorpick.addPatch(pickPatch);
            surface.setSurfaceTexture(selectedHorpick.getPatchColorSurfaceTexture(), glPanel);
            runHorpick.panel.seedSelected(pickPatch);
        }
        model.win3dDisplay();
        return true;
    }

    public void runPickSurface()
    {
        Runnable runPickSurface = new Runnable()
        {
            public void run()
            {
                pickSurfaceProcess();
            }
        };

        pickSurfaceThread = new Thread(runPickSurface);
        pickSurfaceThread.start();
//        runPickSurface.run();
    }

    private void pickSurfaceProcess()
    {
        try
        {
    //        memoryCheckFlag = MEMORY_CHECK_INITIAL;
	        StsObjectRefList patches = selectedHorpick.getPatches();
	        int nPatches = patches.getSize();
			if(nPatches == 0)
			{
				new StsMessage(model.win3d, StsMessage.WARNING, "No seeds have been picked: can't run.  Pick seeds first, please.");
				return;
			}
	        if(seismicVolume == null)
			{
				StsException.systemError("StsHorpickWizard failed. No line2d has been defined.");
				return;
			}
            seismicVolume.clearCache();
//            line2d.setupSliceBlocks();

            rowMin = surface.getRowMin();
            rowMax = surface.getRowMax();
            colMin = surface.getColMin();
            colMax = surface.getColMax();

            nRows = rowMax - rowMin + 1;
            nCols = colMax - colMin + 1;

            checkInitializeArrays();
            checkInitializeSurface();
            isTestPatch = false;
//            setCorCoefsColorscaleRange();

            if(isIterative())
            {
                iterativePickSurface();
//                deleteTraces();
            }
            else
            {
                minCorrel = selectedHorpick.getManualCorMin();
//                corCoefsColorscalePanel.getColorscale().setRange(minCorrel, 1.0f);

                if(selectAllPatches) // indicates all patches to be run
                {
                    boolean displayOK = iterationDisplayInterval == 0;
                    for(int n = 0; n < nPatches; n++)
                    {
//                        setStatusMessage("Processing Patch " + n + " of " + nPatches);
                        pickPatch = (StsPickPatch)patches.getElement(n);
                        pickSurface();
                        if(displayOK) displayPickSurface();
                    }
                    if(!displayOK) displayPickSurface();
                }
                else // run just the current selected patch
                {
                    pickSurface();
                    displayPickSurface();
                }
                surface.computeZRange();
//                setStatusMessage("Constructing Grid...");
//                selectedHorpick.saveSurface();
//                model.win3dDisplay();
            }
			processedSurfaces.add(surface);
        }
        catch(Exception e)
        {
            StsException.outputException("StsHorpickWizard.iterativeProcess() failed.", e, StsException.WARNING);
        }
        finally
        {
            setStatusMessage("AutoTrack picking completed.");
            StsMessageFiles.infoMessage("AutoTrack picking completed.");
            runHorpick.panel.pickCompleted();
            printMemorySummary();
//            if(line2d != null) line2d.deleteSliceBlocks();
            deleteTraces();
        }
    }

    private void printMemorySummary()
    {
        if(seismicVolume != null) seismicVolume.printMemorySummary();
        int memoryUsed = nAllocatedTraces*traceInstanceByteSize;
        StsMessageFiles.infoMessage("Number Traces in memory: " + nAllocatedTraces + " memory used: " + memoryUsed);
        selectedHorpick.printMemorySummary();
    }

    public boolean isIterative()
    {
        if(selectedHorpick == null) return false;
        return selectedHorpick.getIsIterative();
    }
/*
    private void setCorCoefsColorscaleRange()
    {
        float minCor = getCurrentMinCorrel();
        corCoefsColorscalePanel.getColorscale().setRange(minCor, 1.0f);
        corCoefsColorscalePanel.reinitialize();
    }
*/
    public float getMinCorrel()
    {
        if(selectedHorpick == null) return 0.0f;
        return selectedHorpick.getMinCorrel();
    }

    public float getMinCorFilter()
    {
        if(selectedHorpick == null) return 0.0f;
        return selectedHorpick.getMinCorrelFilter();
    }

    public void iterativePickSurface()
    {
        if(selectedHorpick == null) return;

        if(selectAllPatches)
        {
            StsObjectRefList patches = selectedHorpick.getPatches();
            int nPatches = patches.getSize();
            int nIterations = 1 + (int)((selectedHorpick.autoCorMax - selectedHorpick.autoCorMin)/selectedHorpick.autoCorInc);
            float cor = selectedHorpick.autoCorMax;
            for(int i = 0; i < nIterations; i++)
            {
                // if interval is 0, display every patch on every iteration
                boolean displayOK = iterationDisplayInterval == 0;
                minCorrel = cor;
                for(int n = 0; n < nPatches; n++)
                {
//                    setStatusMessage("Patch " + pickPatchIndex + " of " + nPatches + " cor: " + cor);
                    pickPatch = (StsPickPatch)patches.getElement(n);
                    pickSurface();
                    if(displayOK) displayPickSurface();
                }
                cor -= selectedHorpick.autoCorInc;
                // if interval > 0, display only at the end of the 0th, each ith interval, and the last
                if(!displayOK)
                {
                    displayOK = nIterations == 0 || i%iterationDisplayInterval == 0 || i == nIterations-1;
                    if(displayOK) displayPickSurface();
                }
            }
        }
        else // run just the current selected patch
        {
            float cor = selectedHorpick.autoCorMax;
            boolean displayOK = iterationDisplayInterval == 0;
            while(cor >= selectedHorpick.autoCorMin)
            {
                minCorrel = cor;
                pickSurface();
                if(displayOK) displayPickSurface();
                cor -= selectedHorpick.autoCorInc;
            }
            if(!displayOK) displayPickSurface();
        }
        surface.computeZRange();
    }

    private boolean pickSurface()
    {
        pickPatchIndex = pickPatch.getIndex();
        currentPatchNumber = pickPatch.getIndex();
//        selectedHorpick.setMinCorrel(cor);
        initializePickPatchParameters();
        if(!autoPickSurface(selectedHorpick, pickPatch, debug))
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "Auto-Picker failed for patch " + currentPatchNumber);
            return false;
        }
//        setStatusMessage("Patch " + pickPatchIndex + " of " + nPatches + ": constructing Grid.");
        return true;
    }

    private void displayPickSurface()
    {
        surface.initializeSurface();
//        surface.interpolateDistanceTransform();
        surface.setIsVisible(true);
        model.win3dDisplay();
    }

    private void deleteTraces()
    {
        traces = null;
        nAllocatedTraces = 0;
    }

    public void deleteHorpick()
    {
        selectedHorpick = null;
    }

    /** Set the min and max fractions of maxMemory which
     *  is used by the line2d in allocating data planes.
     *  MaxFraction is the fraction of totalMemory at which
     *  planes will be deallocated until the minFraction is reached.
     */
/*
    private void setMemoryFractions()
    {
        if(line2d != null) line2d.setMemoryFractions(minFraction, maxFraction);
    }

    private void restoreMemoryFractions()
    {
        if(line2d != null) line2d.restoreMemoryFractions();
    }
*/
//    public void setRerun(boolean rerun) { this.rerun = rerun; }

    public boolean computeTestCorrelRange()
    {
        if(pickPatch == null) return false;

        surface.setIsVisible(false);
//        if(pickPatch.correlTestMin != 0.0f) return true;
        isTestPatch = true;

//        byte pickType = runHorpick.panel.getPickType();
//        pickPatch.setPickType(pickType);
        StsGridPoint gridPoint = pickPatch.getSeedPoint();
        this.iCenter = gridPoint.row;
        this.jCenter = gridPoint.col;

        nRows = surface.getNRows();
        rowMin = gridPoint.row - 20;
        rowMax = rowMin + 40;
        if(rowMin < 0)
        {
            rowMin = 0;
            rowMax = Math.min(40, nRows-1);
        }
        else if(rowMax >= nRows)
        {
            rowMax = nRows - 1;
            rowMin = Math.max(rowMax-40, 0);
        }
        nCols = surface.getNCols();
        colMin = gridPoint.col - 20;
        colMax = colMin + 40;
        if(colMin < 0)
        {
            colMin = 0;
            colMax = Math.min(40, nCols-1);
        }
        else if(colMax >= nCols)
        {
            colMax = nCols - 1;
            colMin = Math.max(colMax-40, 0);
        }

        seismicVolume.clearCache();
//        line2d.setupSliceBlocks();

//        checkInitializeArrays();
        checkInitializeSurface();
        corCoefs = new float[nRows][nCols];
        patchIndexes = new byte[nRows*nCols];

        pickPatchIndex = pickPatch.getIndex();
        currentPatchNumber = pickPatch.getIndex();
        initializeTestPickPatchParameters();
        if(!autoPickSurface(selectedHorpick, pickPatch, false)) return false;

        float corMin = largeFloat;
        float corMax = -largeFloat;
        boolean hasCor = false;
        for(int row = rowMin; row <= rowMax; row++)
        {
            for(int col = colMin; col <= colMax; col++)
            {
                float cor = corCoefs[row][col];
                if(cor != 0.0f)
                {
                    corMin = Math.min(corMin, cor);
                    corMax = Math.max(corMax, cor);
                    hasCor = true;
                }
            }
        }

        pickPatch.correlTestMin = corMin;
        pickPatch.correlTestMax = corMax;
        setTestPickPatchParameters();

        // reset seed point so not used in interpolation
        pointsNull[iCenter][jCenter] = SURF_GAP_NOT_FILLED;

        deleteTraces();
		surface.setIsVisible(true);

        return true;
    }

    private void initializeTestPickPatchParameters()
    {
        windowSize = defaultWindowSize;
        maxPickDif = 5*defaultMaxPickDif;
    }

    private void initializePickPatchParameters()
    {

        windowSize = pickPatch.getWindowSize();
        maxPickDif = pickPatch.getMaxPickDif();
    }

    private void setTestPickPatchParameters()
    {
        avgTestPickDif /= nTestPickDifs;
        pickPatch.avgTestPickDif = (float)avgTestPickDif;
        pickPatch.maxTestPickDif = maxTestPickDif;
//        maxPickDif = (float)StsMath.ceiling(maxTestPickDif);
        setMaxPickDif(maxTestPickDif);
        pickPatch.minWaveLength = minWaveLength;
        pickPatch.maxWaveLength = maxWaveLength;
        if(maxWaveLength != 0.0f) setWindowSize(StsMath.ceiling(2*maxWaveLength));
    }

    private boolean autoPickSurface(StsHorpick horpick, StsPickPatch pickPatch, boolean debug)
    {
        StsTimer timer = null;
        if(debug)
        {
            timer = new StsTimer();
            timer.start();
        }
/*
        if(rerun == true)
        {
            surface.deleteDisplayLists();
        }
*/
/*
        pickType = pickPatch.getPickType();

        if(isTestPatch) maxPickDif = 5*maxPickDif;

        windowSize = pickPatch.getWindowSize();
*/
        windowSize = 2*(windowSize/2) + 1;
        maxWindowLength = windowSize + 2*StsMath.ceiling(maxPickDif) + 2*maxLag;

        traceInstanceByteSize = 36 + maxWindowLength;

//        nLags = 2*maxLag + 1;

        floatValues = new float[maxWindowLength];

//        pointsZ = surface.getPointsZ();
//        pointsNull = surface.getPointsNull();

//        patchIndexes = selectedHorpick.getPatchIndexes();
//        corCoefs = selectedHorpick.getCorCoefs();

//        this.line2d = horpick.getSeismicVolume();
        nSamples = seismicVolume.getNSlices();

        avgTestPickDif = 0;
        maxTestPickDif = 0.0f;
        nTestPickDifs = 0;
        minWaveLength = largeFloat;
        maxWaveLength = 0.0f;

        if(useTraces && traces == null) traces = new Trace[nRows][nCols];

        StsGridPoint gridPoint = pickPatch.getSeedPoint();

        this.iCenter = gridPoint.row;
        this.jCenter = gridPoint.col;

        float pickZ = pickPatch.getPickZ();
        float nPickF = seismicVolume.getSliceCoor(pickZ);

//        unsignedByteAverage = line2d.getUnsignedByteAverage();

        // adjust pick depending on pickType
        nPickF = findNearestSeedPick(iCenter, jCenter, nPickF);
        if(nPickF == nullValue) return false;
        pickZ = seismicVolume.getZCoor(nPickF);
//        setPickZ(pickZ);
        pickPatch.setPickZ(pickZ);
        patchIndexes[iCenter*nCols+jCenter] = (byte)pickPatchIndex;

        surface.interpolateNeighbors(iCenter, jCenter, pickZ);

        nPatchTraces = 0;
        pickSurface(nPickF);
        if(debug) timer.stopPrint("pick surface: ");
//        if(nTestPickDifs > 0) avgTestPickDif /= nTestPickDifs;
        pickPatch.setNPatchTraces(nPatchTraces);
//        if(isTestPatch) maxPickDif = maxPickDif/5;

        return true;
    }

    public void pickPreferences()
    {
        StsTextAreaDialog dialog = new StsTextAreaDialog(model.win3d, "Display Update Interval:", new Integer(iterationDisplayInterval).toString(), 1, 40);
        dialog.setVisible(true);
        iterationDisplayInterval = new Integer(dialog.getText()).intValue();
    }

    public void deletePatch()
    {
        if(selectAllPatches) deleteAllPatches();
        else              deletePatch(pickPatch);
        model.win3dDisplay();
    }

    private void deleteAllPatches()
    {
        selectedHorpick.deleteAllPatches();
        runHorpick.panel.deleteAllPicks();
        deleteArrays();
        deleteSurface();
    }

    private void deletePatch(StsPickPatch pickPatch)
    {
        selectedHorpick.deletePatch(pickPatch);
        runHorpick.panel.deletePick(pickPatch);
        surface.initializeSurface();
    }

    public void deletePicks()
    {
        int nPatches = selectedHorpick.getPatches().getSize();
        if(nPatches == 0) return;

        if(selectAllPatches || nPatches == 1)
        {
            deleteArrays();
            deleteSurface();
        }
        else
        {
            selectedHorpick.deletePicks(pickPatch);
            surface.initializeSurface();
        }
        model.win3dDisplay();
    }

    private void checkInitializeArrays()
    {
        selectedHorpick.checkInitializeArrays();
        patchIndexes = selectedHorpick.getPatchIndexes();
        corCoefs = selectedHorpick.getCorCoefs();
    }

    private void checkInitializeSurface()
    {
        if(surface == null) return;
        surface.checkInitializeGrid(StsSurface.SURF_GAP_NOT_FILLED, true);
        pointsZ = surface.getPointsZ();
        pointsNull = surface.getPointsNull();
    }

    private void deleteSurface()
    {
        if(selectedHorpick == null) return;
        selectedHorpick.deleteSurface();
        pointsZ = null;
        pointsNull = null;
        surface.initializeSurface();
    }

    private void deleteArrays()
    {
        if(selectedHorpick == null) return;
        selectedHorpick.deleteArrays();
        patchIndexes = null;
        corCoefs = null;
    }

    private void pickSurface(float nPickF)
//    private void pickSurface(float nPickF, boolean stopSpiral)
    {
        boolean doBottom, doTop, doLeft, doRight;
        int i, j;

        // turn on the wait cursor
        StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);

        try
        {
/*
            windowSize = 2*(windowSize/2) + 1;
            maxWindowLength = windowSize + 2*StsMath.ceiling(maxPickDif) + 2*maxLag;
            nLags = 2*maxLag + 1;

            byteValues = new byte[maxWindowLength];
*/
            boolean forwardAdded = true, reverseAdded = true;

            int imin = iCenter;
            int imax = iCenter;
            int jmin = jCenter;
            int jmax = jCenter;

            int nSpiral = 0;
            while(forwardAdded || reverseAdded) // do forward and reverse loops until neither are adding points
            {
                forwardAdded = reverseAdded; // if reverse added points, do another forward loop
                while(forwardAdded)
//                while(forwardAdded || !stopSpiral)
                {
                    if(stopThread) return;

                    boolean terminate = true;
                    doBottom = imin > rowMin;
                    if(doBottom) { imin--; terminate = false; }
                    doTop = imax < rowMax;
                    if(doTop) { imax++; terminate = false; }
                    doLeft = jmin > colMin;
                    if(doLeft) { jmin--;  terminate = false; }
                    doRight = jmax < colMax;
                    if(doRight) { jmax++;  terminate = false; }
                    if(terminate) break;

                    if(debug) System.out.println("Forward spiral: " + nSpiral);
                    forwardAdded = false;

                    if(doBottom || doLeft) // lower-left corner
                        forwardAdded = addPoint(imin, jmin, StsSeismicVolume.YDIR, forwardBL1, forwardAdded);

                    if(doBottom)
                        for(j = jmin+1; j < jmax; j++)
                            forwardAdded = addPoint(imin, j, StsSeismicVolume.YDIR, forwardBottom, forwardAdded);

                    if(doBottom || doRight) // lower right corner
                        forwardAdded = addPoint(imin, jmax, StsSeismicVolume.XDIR, forwardBR1, forwardAdded);

                    if(doRight)
                        for(i = imin+1; i < imax; i++)
                            forwardAdded = addPoint(i, jmax, StsSeismicVolume.XDIR, forwardRight, forwardAdded);

                    if(doRight || doTop) // upper right corner
                        forwardAdded = addPoint(imax, jmax, StsSeismicVolume.YDIR, forwardTR1, forwardAdded);

                    if(doTop)
                        for(j = jmax-1; j > jmin; j--)
                            forwardAdded = addPoint(imax, j, StsSeismicVolume.YDIR, forwardTop, forwardAdded);

                    if(doTop || doLeft) // upper left corner
                            forwardAdded = addPoint(imax, jmin, StsSeismicVolume.XDIR, forwardTL1, forwardAdded);

                    if(doLeft)
                        for(i = imax-1; i > imin; i--)
                            forwardAdded = addPoint(i, jmin, StsSeismicVolume.XDIR, forwardLeft, forwardAdded);

                    if(doBottom || doLeft) // lower-left corner
                        forwardAdded = addPoint(imin, jmin, StsSeismicVolume.YDIR, forwardBL2, forwardAdded);

                    if(doBottom || doRight) // lower right corner
                        forwardAdded = addPoint(imin, jmax, StsSeismicVolume.XDIR, forwardBR2, forwardAdded);

                    if(doRight || doTop) // upper right corner
                        forwardAdded = addPoint(imax, jmax, StsSeismicVolume.YDIR, forwardTR2, forwardAdded);

                    if(doTop || doLeft) // upper left corner
                            forwardAdded = addPoint(imax, jmin, StsSeismicVolume.XDIR, forwardTL2, forwardAdded);

                    if(forwardAdded) reverseAdded = true; // if any points added in forward sweep, do a reverse sweep

                    nSpiral++;
                    setStatusMessage("Patch " + currentPatchNumber + " cor: " + minCorrel + " spiral: " + nSpiral);
                }
                int nForwardSpirals = nSpiral;
                StsMessageFiles.infoMessage("Patch " + currentPatchNumber + " cor: " + minCorrel + " ran  " + nForwardSpirals + " forward spirals");

//                if(mainDebug) System.out.println("Ran " + nSpiral + " forward spirals for patch: " + pickPatchIndex + " pointsAdded: " + forwardAdded);


                reverseAdded = forwardAdded;
                while(true) // if forward added points, do a reverse
//                while(reverseAdded) // if forward added points, do a reverse
//                while(reverseAdded || !stopSpiral) // if forward added points, do a reverse
                {
                    if(stopThread) return;

                    boolean terminate = true;
                    doBottom = imin < iCenter-1;
                    if(doBottom) { imin++; terminate = false; }
                    doTop = imax > iCenter+1;
                    if(doTop) { imax--; terminate = false; }
                    doLeft = jmin < jCenter-1;
                    if(doLeft) { jmin++;  terminate = false; }
                    doRight = jmax > jCenter+1;
                    if(doRight) { jmax--;  terminate = false; }
                    if(terminate) break;

                    if(debug) System.out.println("Reverse spiral: " + nSpiral);

                    reverseAdded = false;

                    if(doBottom || doLeft) // lower-left corner
                        reverseAdded = addPoint(imin, jmin, StsSeismicVolume.YDIR, reverseBL, reverseAdded);

                    if(doLeft)
                        for(i = imin+1; i < imax; i++)
                            reverseAdded = addPoint(i, jmin, StsSeismicVolume.XDIR, reverseLeft, reverseAdded);

                    if(doTop || doLeft) // upper left corner
                            reverseAdded = addPoint(imax, jmin, StsSeismicVolume.XDIR, reverseTL, reverseAdded);

                    if(doTop)
                        for(j = jmin+1; j < jmax; j++)
                            reverseAdded = addPoint(imax, j, StsSeismicVolume.YDIR, reverseTop, reverseAdded);

                    if(doRight || doTop) // upper right corner
                        reverseAdded = addPoint(imax, jmax, StsSeismicVolume.YDIR, reverseTR, reverseAdded);

                    if(doRight)
                        for(i = imax-1; i > imin; i--)
                            reverseAdded = addPoint(i, jmax, StsSeismicVolume.XDIR, reverseRight, reverseAdded);

                    if(doBottom || doRight) // lower right corner
                        reverseAdded = addPoint(imin, jmax, StsSeismicVolume.XDIR, reverseBR, reverseAdded);

                    if(doBottom)
                        for(j = jmax-1; j > jmin; j--)
                            reverseAdded = addPoint(imin, j, StsSeismicVolume.YDIR, reverseBottom, reverseAdded);

                    setStatusMessage("Patch " + currentPatchNumber + " cor: " + minCorrel + " spiral: " + nSpiral);
                    nSpiral--;
                }
                StsMessageFiles.infoMessage("Patch " + currentPatchNumber + " cor: " + minCorrel + " ran  " + nForwardSpirals + " reverse spirals");
//                if(mainDebug) System.out.println("Ran " + nSpiral + " reverse spirals for patch: " + pickPatchIndex + " pointsAdded: " + reverseAdded);

                // run only one forward and reverse spiral; may want to change this logic later
                break;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsAutoPickSurface.pickSurface() failed." +
                " for row: " + iCenter + " of " + rowMin + "-" + rowMax +
                " col: " + jCenter + " of " + colMin + "-" + colMax,
                e, StsException.WARNING);
        }
        finally
        {
			// turn off the wait cursor
			if(cursor != null) cursor.restoreCursor();
            floatValues = null;
        }
    }
    /**
     * Try to add a pickPoint at this ij using other traces at offsets from this ij
     *
     * @param i row of trace pick to add
     * @param j col of trace pick to add
     * @param dir direction of seismic plane containing this trace (XDIR:col YDIR:row)
     * @param traceOffsets row-col offsets to trace(s) to be used
     * @param pointsAdded boolean indicating points have been added; set to true if point added,
     * otherwise pass thru pointsAdded passed in
     * @return pointsAdded
     */
    private boolean addPoint(int i, int j, int dir, int[][] traceOffsets, boolean pointsAdded)
    {
        Trace prevTrace, nextTrace;
        byte surfaceNullValue = pointsNull[i][j];
        if(surfaceNullValue == SURF_NULL)
        {
            boolean isSamePatch = (patchIndexes[i*nCols + j] == (byte)pickPatchIndex);
            if(isSamePatch) return true;
            else            return pointsAdded;
        }

        // get possible range of picks
        int nTraceOffsets = traceOffsets.length;
        float[] prevNPickFs = new float[nTraceOffsets];
        boolean[] isPickGood = new boolean[nTraceOffsets];

        float minPickRangeF = largeFloat;
        float maxPickRangeF = -largeFloat;
        int nPrevPicks = 0;
        for(int n = 0; n < traceOffsets.length; n++)
        {
            int i0 = i + traceOffsets[n][0];
            int j0 = j + traceOffsets[n][1];
            if(!inRange(i0, j0)) continue;
            if(pointsNull[i0][j0] != SURF_NULL) continue;
            if(patchIndexes[i0*nCols + j0] != (byte)pickPatchIndex) continue;

            float prevNPickF;
            if(useTraces)
            {
                if(traces[i0][j0] != null)
                    prevNPickF = traces[i0][j0].nPickF;
                else
                {
                    prevNPickF = seismicVolume.getSliceCoor(pointsZ[i0][j0]);
                    traces[i0][j0] = new Trace(i0, j0, prevNPickF);
                }
            }
            else
                prevNPickF = seismicVolume.getSliceCoor(pointsZ[i0][j0]);

            isPickGood[n] = true;
            prevNPickFs[n] = prevNPickF;
            minPickRangeF = Math.min(prevNPickF, minPickRangeF);
            maxPickRangeF = Math.max(prevNPickF, maxPickRangeF);
            nPrevPicks++;
        }
        if(nPrevPicks == 0) return pointsAdded;
        if(maxPickRangeF - minPickRangeF > 2*maxPickDif) return pointsAdded;

        nextTrace = getSetTrace(i, j, minPickRangeF, maxPickRangeF);
//        setSampleRange(minPickRangeF)        nextTrace.classInitialize(i, j, dirNo);
        if(!nextTrace.computePossiblePicks(minPickRangeF, maxPickRangeF)) return pointsAdded;
        if(debug)
        {
            System.out.print("    Picks on new trace at row: " + i + " col: " + j);
            for(int n = 0; n < nPossiblePicks; n++)
                System.out.print(" " + possiblePicks[n]);
            System.out.println();
        }
        // check that other traces agree on a common pick: they all must be with + or - of MaxPickDif
        float currentPickF = nullValue;
        float nearestPickF;
        int nGoodPicks = 0;
        for(int n = 0; n < nTraceOffsets; n++)
        {
            if(!isPickGood[n]) continue;
            if(debug)
            {
                int i0 = i + traceOffsets[n][0];
                int j0 = j + traceOffsets[n][1];
                System.out.println("     Neighbor pick at row: " + i0 + " col: " + j0 + " is at: " + prevNPickFs[n]);
            }
            nearestPickF = findNearestPick(prevNPickFs[n]);
            if(nearestPickF == nullValue)
            {
                isPickGood[n] = false;
                continue;
            }

            float pickDif = Math.abs(nearestPickF - prevNPickFs[n]);
            if(pickDif > maxPickDif)
            {
                isPickGood[n] = false;
                continue;
            }

            // if other traces don't agree on common pick, return
            if(currentPickF != nullValue)
            {
                if(nearestPickF != currentPickF)
                {
                    if(debug) System.out.println("    Neighbor traces disagree on common pick at row: " + i + " col: " + j);
                    return pointsAdded;
                }
            }
            else
                currentPickF = nearestPickF;

            isPickGood[n] = true;
            nGoodPicks++;
            maxTestPickDif = Math.max(maxTestPickDif, pickDif);
            avgTestPickDif += Math.abs(pickDif);
            nTestPickDifs++;
        }

        if(nGoodPicks == 0)
        {
            if(debug) System.out.println("    Failed pick at row: " + i + " col: " + j + " no good neighbor picks ");
            return pointsAdded;
        }

        nextTrace.nPickF = currentPickF;
//        if(saveCenteredValues) nextTrace.computeCenteredValues();

        float bestCor = corCoefs[i][j];
        float newCor;
        for(int n = 0; n < nTraceOffsets; n++)
        {
            if(!isPickGood[n]) continue;
            int i0 = i + traceOffsets[n][0];
            int j0 = j + traceOffsets[n][1];
            float prevNPickF = prevNPickFs[n];

            prevTrace = null;
            if(useTraces) prevTrace = traces[i0][j0];
            if(prevTrace == null)
            {
                prevTrace = getSetTrace(i0, j0, minPickRangeF, maxPickRangeF);
                prevTrace.nPickF = prevNPickF;
//                if(saveCenteredValues) prevTrace.computeCenteredValues();
            }

            newCor = computeCrossCorrelation(prevTrace, nextTrace);
            if(!isTestPatch && newCor < minCorrel)
            {
                if(debug) System.out.println("    Failed pick at row: " + i + " col: " + j + " correl: " + newCor);
                return pointsAdded;
            }
            bestCor = Math.max(bestCor, newCor);
            corCoefs[i0][j0] = Math.max(corCoefs[i0][j0], newCor);
        }
        corCoefs[i][j] = bestCor;
//        corMin = Math.min(corMin, bestCor);
//        corMax = Math.max(corMax, bestCor);

        if(isTestPatch) return true;

/*
        int fn = Math.round(99*(currentPickF - (int)currentPickF));
        if(fn < 0 || fn > 99)
        {
            System.out.println("position stat: " + fn + " nPickF: " + currentPickF);
        }
        else
            posStat[fn] += 1;

        int cn = Math.round(99*bestCor);
        if(cn < 0 || cn > 99)
        {
            System.out.println("correl coef: " + bestCor);
            cn = StsMath.minMax(cn, 0, 99);
        }
        else
            corStat[cn] += 1;
*/
        // set this pick in arrays
        float z = seismicVolume.getZCoor(currentPickF);
//        if(StsMath.sameAsTol(z, pointsZ[i][j])) return pointsAdded;

        patchIndexes[i*nCols + j] = (byte)pickPatchIndex;
        pointsZ[i][j] = z;
        surface.interpolateNeighbors(i, j, z);
        pointsNull[i][j] = SURF_NULL;
        if(useTraces)
        {
            traces[i][j] = nextTrace;
//            if(saveCenteredValues) nextTrace.values = null;
            nAllocatedTraces++;
            nPatchTraces++;
            checkMemory(nextTrace);
        }

        if(debug) System.out.println("    Added pick at row: " + i + " col: " + j + " correl: " + corCoefs[i][j] + " z: " + pointsZ[i][j]);

        // check the stop criteria and return whether this patch point has been added to the surface
    /*
        if(surfaceNullValue == NULL_GAP_NOT_FILLED)
            return true;
        else if(stopCriteria == STOP_CRITERIA_REPLACE)
        {
            pointsNull[i][j] = NULL_GAP_NOT_FILLED;
            return true;
        }
        else if(stopCriteria == STOP_CRITERIA_IF_SAME_Z)
        {
            if(StsMath.sameAsTol(pointsZ[i][j], pointsZ[i][j], 1.0f)) return pointsAdded;
            pointsNull[i][j] = NULL_GAP_NOT_FILLED;
            return true;
        }
    */
        return true;
    }

/*
    class SpiralIterator implements Iterator
    {
        int rowCenter;
        int colCenter;
        int nSpirals;
        int nRows;
        int nCols;
        int row;
        int col;
        int nSpiral;
        int rowStart;
        int colStart;

        SpiralIterator(int rowCenter, int colCenter, int nSpirals, int nRows, int nCols)
        {
            this.rowCenter = rowCenter;
            this.colCenter = colCenter;
            this.nSpirals = nSpirals;
            this.nRows = nRows;
            this.nCols = nCols;
            row = rowCenter;
            col = colCenter;
            nSpiral = 0;
            rowStart = rowCenter;
            colStart = colCenter;
        }

        public boolean hasNext()
        {
            if(isAtStart())
            {
                nSpiral++;
                if(nSpiral == nSpirals) return false;
                rowStart =
        }
    }
*/
    private float findNearestSeedPick(int row, int col, float nPickF)
    {
        //float pick = nullValue;
//        Trace pickTrace = new Trace(row, col, nPickF, StsSeismicVolume.YDIR);

        //byte[] planeBytes = seismicVolume.readBytePlaneData(StsSeismicVolume.YDIR, seismicVolume.getYCoor(row));
//        float[] values = new float[nSamples];
//        getSignedValues(planeBytes, col*nSamples, nSamples, values);

        Trace pickTrace = new Trace(row, col, nPickF, nSamples);

        float minPickF = Math.max(nPickF - maxPossibleWindow/2, 0.0f);
        float maxPickF = Math.min(nPickF + maxPossibleWindow/2, nSamples-1);

        //int nStart = Math.min((int)minPickF - 1, 0);
        //int nEnd = Math.max((int)maxPickF + 1, nSamples-1);
        //int nValues = nEnd - nStart + 1;

        float nearestPick = nullValue;

        if(pickTrace.computePossiblePicks(minPickF, maxPickF))
        {
            float distance = largeFloat;
            for(int n = 0; n < nPossiblePicks; n++)
            {
                float newDistance = Math.abs(nPickF - possiblePicks[n]);
                if(distance <= newDistance) continue;
                distance = newDistance;
                nearestPick = possiblePicks[n];
            }
        }

        if(nearestPick == nullValue)
            new StsMessage(null, StsMessage.WARNING, "Failed to find a seed pick on this trace");
        else
        {
//            if(saveCenteredValues) pickTrace.computeCenteredValues();
            if(useTraces) traces[row][col] = pickTrace;
        }

        return nearestPick;
    }

    /** Given a pick and a window size (length), set the range of samples we want */
/*
    private void setSampleRange(float nPickF)
    {
        nFirstSample = StsMath.floor(nPickF) - maxWindowLength/2;
        nFirstSample = Math.max(0, nFirstSample);
        nLastSample = nFirstSample + maxWindowLength - 1;
        nLastSample = Math.min(nLastSample, nSamples - 1);
        length = nLastSample - nFirstSample + 1;
    }
*/

    private float findNearestMax(Trace trace)
    {
        float max = nullValue;
        float thisMax;
        int nFirstSample = trace.nFirstSample;
        int length = trace.length;
        float[] values = trace.getFloatValues();
        float distance = largeFloat;

        for(int n = 1; n < length-1; n++)
        {
            thisMax = nullValue;
            if(values[n-1] < values[n])
            {
                if(values[n] > values[n+1])
                {
                    thisMax = StsMath.findQuadraticMinMax(values[n-1], values[n], values[n+1], false);
                    if(thisMax == nullValue) continue;
                    thisMax += n + nFirstSample;
                }
                else if(values[n] == values[n+1])
                {
                    int n0 = n;
                    for(; n < length-1; n++)
                    {
                        if(values[n+1] != values[n])
                        {
                            if(values[n+1] < values[n]) thisMax = (float)(n + n0)/2 + nFirstSample; // values[n+1] < values[n]
                            break;
                        }
                    }
                }
            }
            if(thisMax != nullValue)
            {
                float newDistance =  Math.abs(thisMax - trace.nPickF);
                if(newDistance < distance)
                {
                    max = thisMax;
                    distance = newDistance;
                }
            }
        }
        return max;
    }

    private float findNearestMin(Trace trace)
    {
        float max = nullValue;
        float thisMax;
        int nFirstSample = trace.nFirstSample;
        int length = trace.length;
        float[] values = trace.getFloatValues();
        float distance = largeFloat;

        for(int n = 1; n < length-1; n++)
        {
            thisMax = nullValue;
            if(values[n-1] > values[n])
            {
                if(values[n] < values[n+1])
                {
                    thisMax = StsMath.findQuadraticMinMax(values[n-1], values[n], values[n+1], false);
                    if(thisMax == nullValue) continue;
                    thisMax += n + nFirstSample;
                }
                else if(values[n] == values[n+1])
                {
                    int n0 = n;
                    for(; n < length-1; n++)
                    {
                        if(values[n+1] != values[n])
                        {
                            if(values[n+1] > values[n]) thisMax = (float)(n + n0)/2 + nFirstSample; // values[n+1] < values[n]
                            break;
                        }
                    }
                }
            }
            if(thisMax != nullValue)
            {
                float newDistance =  Math.abs(thisMax - trace.nPickF);
                if(newDistance < distance)
                {
                    max = thisMax;
                    distance = newDistance;
                }
            }
        }
        return max;
    }

    private float findNearestMinusCrossing(Trace trace)
    {
        return findNearestMinusCrossing(trace, StsParameters.PLUS_AND_MINUS);
    }

    private float findNearestMinusCrossing(Trace trace, int searchDir)
    {
        int nFirstSample = trace.nFirstSample;
        int length = trace.length;
        float[] values = trace.getFloatValues();

        if(searchDir == StsParameters.PLUS)
        {
            int nStart = StsMath.ceiling(trace.nPickF);
            for(int n = nStart; n < values.length-1; n++)
            {
                if(values[n] > 0f && values[n+1] <= 0f)
                    return  n + nFirstSample + values[n]/(values[n] - values[n+1]);
            }
            return nullValue;
        }
        else if(searchDir == StsParameters.MINUS)
        {
            int nStart = StsMath.floor(trace.nPickF);
            for(int n = nStart; n > 0; n--)
            {
                if(values[n-1] > 0f && values[n] <= 0f)
                    return n + nFirstSample - 1 + values[n-1]/(values[n-1] - values[n]);
            }
            return nullValue;
        }
        else // both directions
        {
            float minusCrossing = findNearestMinusCrossing(trace, StsParameters.MINUS);
            float plusCrossing = findNearestMinusCrossing(trace, StsParameters.PLUS);
            if(minusCrossing == nullValue) return plusCrossing;
            else if(plusCrossing == nullValue) return minusCrossing;
            else
            {
                if(trace.nPickF - minusCrossing < plusCrossing - trace.nPickF)
                    return minusCrossing;
                else
                    return plusCrossing;
            }
        }
    }

    private float findNearestPlusCrossing(Trace trace)
    {
        return findNearestPlusCrossing(trace, StsParameters.PLUS_AND_MINUS);
    }

    private float findNearestPlusCrossing(Trace trace, int searchDir)
    {
        int nFirstSample = trace.nFirstSample;
        int length = trace.length;
        float[] values = trace.getFloatValues();

        if(searchDir == StsParameters.PLUS)
        {
            int nStart = StsMath.ceiling(trace.nPickF);
            for(int n = nStart; n < values.length-1; n++)
            {
                if(values[n] < 0f && values[n+1] >= 0f)
                    return n + nFirstSample + values[n]/(values[n] - values[n+1]);
            }
            return nullValue;
        }
        else if(searchDir == StsParameters.MINUS)
        {
            int nStart = StsMath.floor(trace.nPickF);
            for(int n = nStart; n > 0; n--)
            {
                if(values[n-1] < 0f && values[n] >= 0f)
                    return n + nFirstSample - 1 + values[n-1]/(values[n-1] - values[n]);
            }
            return nullValue;
        }
        else // both directions
        {
            float minusCrossing = findNearestPlusCrossing(trace, StsParameters.MINUS);
            float plusCrossing = findNearestPlusCrossing(trace, StsParameters.PLUS);
            if(minusCrossing == nullValue) return plusCrossing;
            else if(plusCrossing == nullValue) return minusCrossing;
            else
            {
                if(trace.nPickF - minusCrossing < plusCrossing - trace.nPickF)
                    return minusCrossing;
                else
                    return plusCrossing;
            }
        }
    }

    private void checkMemory(Trace trace)
    {
        if(memoryCheckFlag == MEMORY_CHECK_NO) return;

        if(memoryCheckFlag == MEMORY_CHECK_INITIAL)
        {
            // If we were unable to allocate all the slices, then memory is limited
            // and we won't save traces.  Notify the user that things will be slow.
        /*
            if(!trace.seismicSliceAllocationOK())
            {
                useTraces = false;
                traces = null;
                new StsMessage(model.win3d, StsMessage.WARNING, "Memory for seismic data is limited: process will be slow.");
                memoryCheckFlag = MEMORY_CHECK_YES;
                return;
            }
        */
            maxMemory = Runtime.getRuntime().maxMemory();

            long totalMemory = Runtime.getRuntime().totalMemory();
            long freeMemory = Runtime.getRuntime().freeMemory();
            long usedMemory = totalMemory - freeMemory;
            long availMemory = maxMemory - usedMemory;
            int traceSize = trace.getSize();

            if(debug) System.out.println("availMemory for trace array: " + availMemory);
            if(availMemory > 2*nRows*nCols*traceSize)
            {
                memoryCheckFlag = MEMORY_CHECK_NO;
                System.out.println("trace array will be used and memory won't be checked as availMemory is sufficient.");
            }
            else
            {
                memoryCheckFlag = MEMORY_CHECK_YES;
                traceCheckInterval = (int)(availMemory/10/traceSize);
                if(debug)
                {
                    int nPoints = nRows*nCols;
                    System.out.println("Horpick memory check will be made every " + traceCheckInterval + " out of " + nPoints);
                }
            }
        }
        else // MEMORY_CHECK_YES
        {
            if(nAllocatedTraces%traceCheckInterval == 0)
            {
            /*
                if(!trace.seismicSliceAllocationOK())
                {
                    useTraces = false;
                    traces = null;
                    new StsMessage(model.win3d, StsMessage.WARNING, "Memory for seismic data is limited: process will be slow.");
                    memoryCheckFlag = MEMORY_CHECK_NO;
                    return;
                }
            */
                long totalMemory = Runtime.getRuntime().totalMemory();
                long freeMemory = Runtime.getRuntime().freeMemory();
                long usedMemory = totalMemory - freeMemory;
                long availMemory = maxMemory - usedMemory;
                int traceSize = trace.getSize();

                // mainDebug testing
//                availMemory = (long)((0.5*nRows*nCols - nAllocatedTraces)*traceSize);

                long neededMemory = (nRows*nCols - nAllocatedTraces)*traceSize;
                if(neededMemory <= 0.8*availMemory) return;

                StsObjectRefList patches = selectedHorpick.getPatches();
                int nPatches = patches.getSize();
                if(nPatches == 1)
                {
                    if(useTraces == false)
                        new StsMessage(model.win3d, StsMessage.WARNING, "Getting low on memory and can't do anything about it.  Proceeding, but expect problems.");
                    else
                    {
                        new StsMessage(model.win3d, StsMessage.WARNING, "Dropping trace cache to reduce memory requirements.  Proceeding, but may have problems.");
                        useTraces = false;
                        deleteTraces();
                    }
                    memoryCheckFlag = MEMORY_CHECK_NO;
                    return;
                }

                int nFreeTracesNeeded = (int)(neededMemory/2/traceSize);
                int nTracesFreeable = 0;
                int currentPatchIndex = pickPatch.getIndex();
                int lowerTraceLimit = currentPatchIndex - 1;
                while(lowerTraceLimit >= 0 && nTracesFreeable < nFreeTracesNeeded)
                {
                    StsPickPatch otherPatch = (StsPickPatch)patches.getElement(lowerTraceLimit);
                    nTracesFreeable += otherPatch.getNPatchTraces();
                    lowerTraceLimit--;
                }
                lowerTraceLimit++;

                int upperTraceLimit = nPatches - 1;
                if(nTracesFreeable < nFreeTracesNeeded)
                {
                    while(upperTraceLimit > currentPatchIndex && nTracesFreeable < nFreeTracesNeeded)
                    {
                        StsPickPatch otherPatch = (StsPickPatch)patches.getElement(lowerTraceLimit);
                        nTracesFreeable += otherPatch.getNPatchTraces();
                        upperTraceLimit--;
                    }
                }
                upperTraceLimit++;

                // get rid of these traces
                int n = 0;
                for(int row = 0; row < nRows; row++)
                {
                    for(int col = 0; col < nCols; col++)
                    {
                        int patchIndex = StsMath.signedByteToUnsignedInt(patchIndexes[n]);
                        if(patchIndex == -1 || patchIndex == currentPatchIndex ||
                           patchIndex < lowerTraceLimit || patchIndex < upperTraceLimit)
                            continue;
                        traces[row][col] = null;
                        nAllocatedTraces--;
                        n++;
                    }
                }

                for(n = 0; n < nPatches; n++)
                {
                    if(n == currentPatchIndex || n < lowerTraceLimit || n < upperTraceLimit) continue;
                    StsPickPatch patch = (StsPickPatch)patches.getElement(n);
                    System.out.println("Removed " + patch.getNPatchTraces() + " traces from patch " + n);
                    patch.setNPatchTraces(0);
                }
            }
        }
    }

    final private boolean inRange(int i, int j)
    {
        return i >= 0 && i < nRows && j >= 0 && j < nCols;
    }

    /**
     * Given possible picks on the new trace at i-j, find the nearest one to the
     * previous pick at i0-j0
     * @return correlation coefficient
     */
    private float findNearestPick(float prevNPickF)
    {
        if(nPossiblePicks == 0) return nullValue;
        if(nPossiblePicks == 1) return possiblePicks[0];
        float dif = Math.abs(prevNPickF - possiblePicks[0]);
        float nearestPick = possiblePicks[0];
        for(int n = 1; n < nPossiblePicks; n++)
        {
            float newDif = Math.abs(prevNPickF - possiblePicks[n]);
            if(newDif < dif)
            {
                nearestPick = possiblePicks[n];
                dif = newDif;
            }
        }
        return nearestPick;
    }
/*
    private void crossCorrelate()
    {
        int center = length/2;
        int halfSize = windowSize/2;
        for(int lag = -maxLag; lag <= maxLag; lag++)
        {
            correl[lag+maxLag] = 0.0;
            for(int n = center-halfSize; n <= center+halfSize; n++)
                correl[lag+maxLag] += prevTrace.values[n]*nextTrace.values[n+lag];
        }
        for(int n = 0; n < nLags; n++)
            correl[n] /= windowSize*Math.sqrt(prevTrace.autoCovariance*nextTrace.autoCovariance);
    }
*/
    private float computeCrossCorrelation(Trace prevTrace, Trace nextTrace)
    {
        if(prevTrace == null || nextTrace == null) return 0.0f;

        float[] prevAdjustedValues, nextAdjustedValues;
//        if(!saveCenteredValues)
        {
            if(staticPrevAdjustedValues == null || prevTrace.length > staticPrevAdjustedValues.length)
                staticPrevAdjustedValues = new float[prevTrace.length];
            prevAdjustedValues = staticPrevAdjustedValues;
            prevTrace.computeCenteredValues(staticPrevAdjustedValues);

            if(staticNextAdjustedValues == null || nextTrace.length > staticNextAdjustedValues.length)
                staticNextAdjustedValues = new float[nextTrace.length];
            nextTrace.computeCenteredValues(staticNextAdjustedValues);
            nextAdjustedValues = staticNextAdjustedValues;
        }
/*
        else
        {
            prevAdjustedValues = prevTrace.adjustedValues;
            nextAdjustedValues = nextTrace.adjustedValues;
        }
*/
        int m1 = prevTrace.nFirstAdjustedValue;
        int p1 = prevTrace.nLastAdjustedValue;
        int c1 = prevTrace.nAdjustedCenter;
        int m2 = nextTrace.nFirstAdjustedValue;
        int p2 = nextTrace.nLastAdjustedValue;
        int c2 = nextTrace.nAdjustedCenter;

        int md = Math.max(m1-c1, m2-c2);
        int pd = Math.min(p1-c1, p2-c2);

        double cor = 0.0;
        int m = 0;
        double prevCovar = 0.0;
        double nextCovar = 0.0;
        int nValues = pd - md + 1;
        for(int n = md; n <= pd; n++, m++)
        {
            float coef = StsMath.simpsonCoefs(m, nValues);


            float prevValue = prevAdjustedValues[n+c1];
            float nextValue = nextAdjustedValues[n+c2];
            cor += coef*prevValue*nextValue;
            prevCovar += coef*prevValue*prevValue;
            nextCovar += coef*nextValue*nextValue;
        }

        if(nValues == 0) return 0;
        prevCovar /= nValues;
        nextCovar /= nValues;
        cor /= nValues*Math.sqrt(prevCovar*nextCovar);
/*
        if(!saveCenteredValues)
        {
            prevTrace.adjustedValues = null;
            nextTrace.adjustedValues = null;
        }
*/
        return (float)cor;
    }
/*
    private void getSignedValues(byte[] bytes, int offset, int length, float[] values)
    {
        int n0, n1, n;
        try
        {
            boolean hasNulls = false;
            for(n = 0; n < length; n++)
            {
                values[n] = (float)StsMath.signedByteToUnsignedInt(bytes[offset + n]);
                if(values[n] == 255f) hasNulls = true;
            }
            if(hasNulls)
            {
                for(n = 0; n < length; n++)
                {
                    if(values[n] == 255f)
                    {
                        boolean nullsFilled = false;
                        n0 = n;
                        n1 = length-1;
                        for(n = n+1; n < length; n++)
                        {
                            if(values[n] != 255f)
                            {
                                n1 = n - 1;
                                break;
                            }
                        }
                        fillNullValues(values, n0, n1);
                        n = n1;
                    }
                }
            }
            double avg = 0;
            for(n = 0; n < length; n++)
                avg += values[n];
            float avgF = (float)(avg/length);
            for(n = 0; n < length; n++)
                values[n] -= avgF;
        }
        catch(Exception e)
        {
            StsException.outputException("getSignedValues() failed.", e, StsException.WARNING);
        }
    }

    private void fillNullValues(float[] values, int nFirstNull, int nLastNull)
    {
        int nValues = values.length;

        if(nFirstNull > 0 && nLastNull < nValues-1)
        {
            float v0 = values[nFirstNull-1];
            float v1 = values[nLastNull+1];
            float df = 1.0f/(nLastNull - nFirstNull + 2);
            float dv = (v1 - v0);
            float f = df;
            for(int n = nFirstNull; n <= nLastNull; n++)
            {
                values[n] = v0 + f*dv;
                f += df;
            }
        }
        else if(nFirstNull == 0 && nLastNull < nValues-1)
        {
            for(int n = nFirstNull; n <= nLastNull; n++)
                values[n] = values[nLastNull+1];
        }
        else if(nFirstNull > 0 && nLastNull == nValues-1)
        {
            for(int n = nFirstNull; n <= nLastNull; n++)
                values[n] = values[nFirstNull-1];
        }
        else // entire interval is null
        {
            for(int n = 0; n < nValues; n++)
                values[n] = 0.0f;
        }
    }
*/
    private Trace getSetTrace(int row, int col, float minPickF, float maxPickF)
    {
        float avgPickF = (minPickF + maxPickF)/2;
        if(useTraces)
        {
            Trace trace = traces[row][col];
            if(trace != null && trace.values != null) return traces[row][col];
        }
        return new Trace(row, col, avgPickF);
    }

    class Trace
    {
        int row, col;
        int nFirstSample;
        int length;
        int dataStart;
        float nPickF;
        byte[] values;
//        float avgF;
//        byte[] adjustedValues = null;
//        float autoCovariance;
        int nAdjustedCenter;
        int nFirstAdjustedValue;
        int nLastAdjustedValue;

		/** pulls a complete trace out of a row-block. Used when getting the picked trace */
        Trace(int row, int col, float nPickF, int nSamples)
        {
            this.row = row;
            this.col = col;
            this.nPickF = nPickF;
            nFirstSample = 0;
            length = nSamples;
			values = new byte[length];
			seismicVolume.getRowBlockTraceValues(row, col, 0, nSamples-1, values);
			convertToSignedBytes();
        }

		/** pulls a partial trace out of a slice-block.  Used for the correlated traces as
		 *  we spiral out from the picked trace.  Using slice blocks limits the amount of
		 *  data which must be in memory.
		 */
		Trace(int row, int col, float nPickF)
        {
            this.row = row;
            this.col = col;
            this.nPickF = nPickF;

            nFirstSample = Math.max(0, StsMath.floor(nPickF) - maxWindowLength/2);
            int nLastSample = Math.min(nSamples-1, nFirstSample + maxWindowLength - 1);
            length = nLastSample - nFirstSample + 1;
            values = new byte[length];
            seismicVolume.getSliceBlockTraceValues(row, col, nFirstSample, nLastSample, values);
			convertToSignedBytes();
       }

		private void convertToSignedBytes()
		{
			for(int n = 0; n < length; n++)

				values[n] = (byte)(StsMath.signedByteToUnsignedInt(values[n]) - 127);
		}

        void computeAvg(byte[] values, int nSamples)
        {
            double avg = 0;
            int n;
            float avgF;

            try
            {
                for(n = 0; n < nSamples; n++)
                {
//                    float value = (float)(values[n] & 0xFF);
                    avg += values[n];
                }
                avgF = (float)(avg/nSamples);
            }
            catch(Exception e)
            {
                StsException.outputException("StsHorpickWizard.Trace.constructor() failed.",
                    e, StsException.WARNING);
            }
        }

        /** Check that all the planes needed for this trace are currently in memory.
         *  Generally called AFTER the trace is constructed to insure that we are
         *  still running in memory.
         */
	/*
        boolean seismicSliceAllocationOK()
        {
            return line2d.traceInMemory(nFirstSample, nLastSample);
        }
    */
        /** return the size of this trace object in bytes */
        public int getSize()
        {
            return (10 + length)*4;
        }
/*
        void computeAutoCovariance()
        {
            autoCovariance = 0.0f;
            for(int n = maxLag; n < windowSize + maxLag; n++)
            {
                float value = (float)(values[n] & 0xFF) - unsignedByteZero;
                autoCovariance += value*value;
            }
            autoCovariance /= windowSize;
        }
*/
        boolean computePossiblePicks(float minPickF, float maxPickF)
        {
            int nStart, nEnd;

            nStart = StsMath.floor(minPickF - 2*maxPickDif - nFirstSample);
            nStart = Math.max(nStart, 1);

            nEnd = StsMath.ceiling(maxPickF + 2*maxPickDif - nFirstSample);
            nEnd = Math.min(nEnd, length-2);

            boolean foundPicks = false;
            switch(pickType)
            {
                case PICK_MAX:
                    foundPicks = findMaxs(nStart, nEnd);
                    break;
                case PICK_MIN:
                    foundPicks = findMins(nStart, nEnd);
                    break;
                case PICK_ZERO_PLUS:
                    foundPicks = findPlusCrossings(nStart, nEnd);
                    break;
                case PICK_ZERO_MINUS:
                    foundPicks = findMinusCrossings(nStart, nEnd);
            }
            if(!foundPicks) return false;

            if(isTestPatch && nPossiblePicks > 1)
            {
                float waveLength = (possiblePicks[nPossiblePicks-1] - possiblePicks[0])/(nPossiblePicks - 1);
                minWaveLength = Math.min(minWaveLength, waveLength);
                maxWaveLength = Math.max(maxWaveLength, waveLength);
            }
            return true;
        }

        /** Called if we wish to save adjusted values in the trace;
         *  otherwise they will be computed as needed and saved in static arrays
         */
    /*
        void computeCenteredValues()
        {
            if(!saveCenteredValues) return;
            if(adjustedValues == null) adjustedValues = new float[length];
            computeCenteredValues(adjustedValues);
        }
    */
        /** Computes evenly spaced samples centered at the pick point
         *  Used in computing cross-correlations between traces.
         */
        void computeCenteredValues(float[] adjustedValues)
        {
            double f1, f2, f3;
            int n;

            try
            {
                nAdjustedCenter = (int)nPickF;
                f1 = nPickF - nAdjustedCenter;
                nAdjustedCenter -= nFirstSample;
                if(adjustedValues == null) adjustedValues = new float[length];
                nFirstAdjustedValue = nAdjustedCenter - windowSize/2;
                nFirstAdjustedValue = Math.max(1, nFirstAdjustedValue);
                nLastAdjustedValue = nAdjustedCenter + windowSize/2;
                nLastAdjustedValue = Math.min(length-3, nLastAdjustedValue);

                f2 = f1*f1;
                f3 = f2*f1;

                float[] floatValues = getFloatValues(nFirstAdjustedValue-1, nLastAdjustedValue+2);

                for(n = nFirstAdjustedValue; n <= nLastAdjustedValue; n++)
                {
                    double a = -floatValues[n-1]/2 + 3*floatValues[n]/2 - 3*floatValues[n+1]/2 + floatValues[n+2]/2;
                    double b = floatValues[n-1] - 5*floatValues[n]/2 + 2*floatValues[n+1] - floatValues[n+2]/2;
                    double c = (floatValues[n+1] - floatValues[n-1])/2;
                    double d = floatValues[n];
                    adjustedValues[n] = (float)(a*f3 + b*f2 + c*f1 + d);
                }
            }
            catch(Exception e)
            {
                StsException.outputException("StsHorpickWizard.computeCenteredValues() failed.", e, StsException.WARNING);
            }
        }

        float[] getFloatValues()
        {
            computeFloatValues(0, length-1);
            return floatValues;
        }

        float[] getFloatValues(int nStart, int nEnd)
        {
            computeFloatValues(nStart, nEnd);
            return floatValues;
        }

        void computeFloatValues(int nStart, int nEnd)
        {
            if(floatValues.length < nEnd+1) floatValues = new float[nEnd+1];
            for(int n = nStart; n <= nEnd; n++)
                floatValues[n] = (float)values[n];
        }

        boolean findMaxs(int nStart, int nEnd)
        {
            float f = 0.0f;
            nPossiblePicks = 0;
            if(possiblePicks == null) possiblePicks = new float[100];

            // first and last slopes same as neighboring, so skip using them
            float distance = largeFloat;
            float newPickF = -1.0f;

            float[] values = getFloatValues(nStart-1, nEnd+1);

            // if we are on a plateau, search for first value off of plateau if it exists
            if(values[nStart-1] == values[nStart])
            {
                nStart--;
                for(; nStart > 0; nStart--)
                    if(values[nStart-1] != values[nStart]) break;
                if(nStart == 0) nStart = 1;
            }

            for(int n = nStart; n <= nEnd; n++)
            {
                float thisMax = nullValue;
                if(values[n-1] < values[n])
                {
                    if(values[n] > values[n+1])
                    {
                        thisMax = StsMath.findQuadraticMinMax(values[n-1], values[n], values[n+1], false);
                        if(thisMax == nullValue) continue;
                        thisMax += n;
                    }
                    else if(values[n] == values[n+1])
                    {
                        int n0 = n;
                        for(; n <= nEnd; n++)
                        {
                            if(values[n+1] != values[n])
                            {
                                if(values[n+1] < values[n]) thisMax = (float)(n + n0)/2; // values[n+1] < values[n]
                                break;
                            }
                        }
                    }
                }
                if(thisMax != nullValue) possiblePicks[nPossiblePicks++] = thisMax + nFirstSample;
            }
            return nPossiblePicks > 0;
        }

        boolean findMins(int nStart, int nEnd)
        {
            float f = 0.0f;
            nPossiblePicks = 0;
            if(possiblePicks == null) possiblePicks = new float[100];

            // first and last slopes same as neighboring, so skip using them
            float distance = largeFloat;
            float newPickF = -1.0f;

            float[] values = getFloatValues(nStart-1, nEnd+1);

            for(int n = nStart; n <= nEnd; n++)
            {
                float thisMax = nullValue;
                if(values[n-1] > values[n])
                {
                    if(values[n] < values[n+1])
                    {
                        thisMax = StsMath.findQuadraticMinMax(values[n-1], values[n], values[n+1], false);
                        if(thisMax == nullValue) continue;
                        thisMax += n;
                    }
                    else if(values[n] == values[n+1])
                    {
                        int n0 = n;
                        for(; n < nEnd; n++)
                        {
                            if(values[n+1] != values[n])
                            {
                                if(values[n+1] > values[n]) thisMax = (float)(n + n0)/2; // values[n+1] < values[n]
                                break;
                            }
                        }
                    }
                }
                if(thisMax != nullValue) possiblePicks[nPossiblePicks++] = thisMax + nFirstSample;
            }
            return nPossiblePicks > 0;
        }
/*
        boolean isDifferent(float newNPickF)
        {
            for(int i = 0; i < nPossiblePicks; i++)
                 if(newNPickF == possiblePicks[i]) return false;
            return true;
        }
*/
        boolean findPlusCrossings(int nStart, int nEnd)
        {
            float f = 0.0f;
            nPossiblePicks = 0;
            if(possiblePicks == null) possiblePicks = new float[100];

            float distance = largeFloat;
            float newPickF = -1.0f;

            float[] values = getFloatValues(nStart, nEnd+1);

            for(int n = nStart; n <= nEnd; n++)
            {
                if(values[n] < 0f && values[n+1] >= 0f)
                {
                    float thisNPickF = n + values[n]/(values[n] - values[n+1]);
                    possiblePicks[nPossiblePicks++] = thisNPickF + nFirstSample;
                }
            }
            return nPossiblePicks > 0;
        }

        boolean findMinusCrossings(int nStart, int nEnd)
        {
            float f = 0.0f;
            nPossiblePicks = 0;
            if(possiblePicks == null) possiblePicks = new float[100];

            float distance = largeFloat;
            float newPickF = -1.0f;

            float[] values = getFloatValues(nStart, nEnd+1);

            for(int n = nStart; n <= nEnd; n++)
            {
                if(values[n] > 0f && values[n+1] <= 0f)
                {
                    float thisNPickF = n + values[n]/(values[n] - values[n+1]);
                    possiblePicks[nPossiblePicks++] = thisNPickF + nFirstSample;
                    float newDistance =  Math.abs(nPickF - thisNPickF);
                    if(newDistance < distance)
                    {
                        newPickF = thisNPickF;
                        distance = newDistance;
                    }
                }
            }
            return nPossiblePicks > 0;
        }
    }
/*
    public void smoothSurface()
    {
    }

    public void fillSurface()
    {
        surface.interpolateDistanceTransform();
    }
*/
/*
    public StsPickPatch setPatchSelected(Color color)
    {
        StsObjectRefList patches = selectedHorpick.getPatches();
        int nPatches = patches.getSize();
        for(int n = 0; n < nPatches; n++)
        {
            StsPickPatch patch = (StsPickPatch)patches.getElement(n);
            if(patch.getColor().equals(color))
            {
                pickPatch = patch;
                pickPatchIndex = pickPatch.index();
                return patch;
            }
        }
        return null;
    }
*/
    public Color[] getPatchColors()
    {
        if(selectedHorpick == null) return null;
//        if(selectedHorpick == null) return new Color[0];
        StsObjectRefList patches = selectedHorpick.getPatches();
        int nPatches = patches.getSize();
        Color[] patchColors = new Color[nPatches];
        for(int n = 0; n < nPatches; n++)
        {
            StsPickPatch patch = (StsPickPatch)patches.getElement(n);
            patchColors[n] = patch.getColor();
        }
        return patchColors;
    }

    public StsPickPatch getPickPatch() { return pickPatch; }
    public void setPickPatch(StsPickPatch pickPatch) { this.pickPatch = pickPatch; }

    public StsColorListItem getPickPatchColorItem()
    {
        if(pickPatch == null) return null;
        return new StsColorListItem(pickPatch);
    }

    public void setPickPatchColorItem(StsColorListItem pickPatchColorListItem)
    {
        pickPatch = (StsPickPatch)pickPatchColorListItem.getObject();
        if(runHorpick != null) runHorpick.panel.setPatchSelected(pickPatch);
        selectAllPatches = (pickPatch == null); // first item in list is "All" with a null patch associated
        String name = pickPatchColorListItem.getName();
        if(name.equals("All") || name.equals("None"))
            currentPatchNumber = -1;
        else
            currentPatchNumber = Integer.parseInt(name);
    }

    public void setStatusMessage(String message)
    {
        runHorpick.setStatusMessage(message);
    }

    public void setPickZ(float pickZ)
    {
        if(pickPatch == null) return;
        pickPatch.setPickZ(pickZ);
        runHorpick.panel.computeTestCorrelRange();
        model.win3dDisplay();
    }

    public float getPickZ()
    {
        if(pickPatch == null) return nullValue;
        return pickPatch.getPickZ();
    }

    public float getPickInline()
    {
        if(pickPatch == null) return nullValue;
        float y = pickPatch.getSeedPoint().getY();
        if(y == nullValue) return nullValue;
        return seismicVolume.getNearestBoundedRowNumFromY(y);
    }

    public void setPickInline(float rowNum)
    {
        if(pickPatch == null) return;
        float y = seismicVolume.getYFromRowNum(rowNum);
        pickPatch.getSeedPoint().setY(y);
        runHorpick.panel.computeTestCorrelRange();
        model.win3dDisplay();
    }

    public float getPickXline()
    {
        if(pickPatch == null) return nullValue;
        float x = pickPatch.getSeedPoint().getX();
        if(x == nullValue) return nullValue;
        return seismicVolume.getNearestBoundedColNumFromX(x);
    }

    public void setPickXline(float colNum)
    {
        if(pickPatch == null) return;
        float x = seismicVolume.getXFromColNum(colNum);
        pickPatch.getSeedPoint().setX(x);
        runHorpick.panel.computeTestCorrelRange();
        model.win3dDisplay();
    }

    public int getWindowSize()
    {
        if(pickPatch != null) windowSize = pickPatch.getWindowSize();
        return windowSize;
    }

    public void setWindowSize(int windowSize)
    {
        this.windowSize = windowSize;
        if(pickPatch == null) return;
        pickPatch.setWindowSize(windowSize);
    }

    public String[] getDisplayPropertyNames()
    {
        StsSeismicVolumeClass volumeClass = (StsSeismicVolumeClass)model.getCreateStsClass(StsSeismicVolume.class);
        String[] volumeNames = volumeClass.getNames();
        String[] propertyNames = new String[volumeNames.length + 3];
        propertyNames[0] = StsHorpick.displayPropertyNone;
        propertyNames[1] = StsHorpick.displayPropertyPatchColor;
        propertyNames[2] = StsHorpick.displayPropertyCorrelCoefs;
        for(int n = 3; n < propertyNames.length; n++)
            propertyNames[n] = volumeNames[n-3];
        return propertyNames;
    }

    public void initializeDisplayPropertyName()
    {
        setDisplayPropertyName(StsHorpick.displayPropertyPatchColor);
    }

    public void setDisplayPropertyName(String propertyName)
    {
        displayPropertyName = propertyName;
        if(selectedHorpick != null)
        {
            selectedHorpick.setDisplayPropertyName(propertyName);
            if(propertyName == StsHorpick.displayPropertyCorrelCoefs)
                finishHorpick.panel.corCoefsColorscaleEnabled(true);
            else
                finishHorpick.panel.corCoefsColorscaleEnabled(false);
        }
    }

    public String getDisplayPropertyName()
    {
        return displayPropertyName;
    }

    public void minCorSliderChanged(float value)
    {
        if(selectedHorpick == null) return;
        if(value == selectedHorpick.getMinCorrelFilter()) return;
        applyFilterCorrelValues(value);
        selectedHorpick.setMinCorrelFilter(value);
        model.win3dDisplay();
    }

    private void applyFilterCorrelValues(float minCorrelFilter)
    {

        if(corCoefs == null || pointsNull == null) return;

        nRows = surface.getNRows();
        nCols = surface.getNCols();

        // remove points whose correl coef is below min value
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
            {
                if(pointsNull[row][col] == StsSurface.SURF_GAP_SET && corCoefs[row][col] >= minCorrelFilter)
                    pointsNull[row][col] = StsSurface.SURF_PNT;
                else if(pointsNull[row][col] == StsSurface.SURF_PNT && corCoefs[row][col] < minCorrelFilter)
                    pointsNull[row][col] = StsSurface.SURF_GAP_SET;
            }

        // create 2d array of bytes and assign 1 to seed points
        // do diagonal sweeps and assign 1 to points which are NOT_NULL && one of its neighbors is 1

        byte[][] ok = new byte[nRows][nCols];

        StsObjectRefList patches = selectedHorpick.getPatches();
        int nPatches = patches.getSize();
        for(int n = 0; n < nPatches; n++)
        {
            StsPickPatch patch = (StsPickPatch)patches.getElement(n);
            StsGridPoint seedPoint = patch.getSeedPoint();
            ok[seedPoint.row][seedPoint.col] = 1;
        }

        boolean converged = false;
        int iter = 0;
        int maxIter = 10;
        boolean debug = false;

        if(debug)
        {
            System.out.println("minCorrelFilter: " + minCorrelFilter);
            System.out.println("    corCoefs[2][[10]: " + corCoefs[2][10]);
            System.out.println("    corCoefs[3][[11]: " + corCoefs[3][11]);
        }
        while(!converged)
        {
            converged = true;

            if(debug) System.out.println("Iteration: " + iter + " Sweep LL to UR");
            for(int row = 0; row < nRows; row++)
                for(int col = 0; col < nCols; col++)
                    if(addFilterPoint(row, col, 1, 1, ok)) converged = false;

            if(debug) System.out.println("Iteration: " + iter + " Sweep UR to LL");
            for(int row = nRows-1; row >= 0; row--)
                for(int col = nCols-1; col >= 0; col--)
                    if(addFilterPoint(row, col, -1, -1, ok)) converged = false;
/*
            if(mainDebug) System.out.println("Iteration: " + iter + " Sweep LR to UL");
            for(int row = 0; row < nRows; row++)
                for(int col = nCols-1; col >= 0; col--)
                    if(!addFilterPoint(row, col, ok)) converged = false;

            if(mainDebug) System.out.println("Iteration: " + iter + " Sweep UL to LR");
            for(int row = nRows-1; row >= 0; row--)
                for(int col = 0; col < nCols; col++)
                    if(!addFilterPoint(row, col, ok)) converged = false;
*/
            if(++iter >= maxIter) converged = true;
        }
        // any point flagged as 0 is an isolated point: remove it
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(ok[row][col] == 0)  pointsNull[row][col] = StsSurface.SURF_GAP_SET;

        surface.setTextureChanged();
//        selectedHorpick.deleteCorrelCoefsTextureDisplayable();
//        surface.setNewTextureDisplayable(selectedHorpick.getCorrelCoefsTextureDisplayable());
    }

    private boolean addFilterPoint(int row, int col, int dRow, int dCol, byte[][] ok)
    {
        boolean debug = false;
    /*
        boolean mainDebug = row == 2 && col == 10 || row == 3 && col == 11;
        if(mainDebug)
        {
            System.out.println(" ok: " + ok[row][col] + " pointsNull: " + pointsNull[row][col]);
            mainDebug = true;
        }
    */
        if(ok[row][col] == 1) return false;
        if(pointsNull[row][col] != StsSurface.SURF_PNT) return false;
        if(!filterPointOK(row, col, dRow, dCol, ok, debug)) return false;
        ok[row][col] = 1;
        return true;
    }

    private boolean filterPointOK(int row, int col, int dRow, int dCol, byte[][] ok, boolean debug)
    {
        if(neighborOK(row-dRow, col-dCol, ok, debug)) return true;
        if(neighborOK(row-dRow, col,      ok, debug)) return true;
        if(neighborOK(row-dRow, col+dCol, ok, debug)) return true;
        if(neighborOK(row,      col-dCol, ok, debug)) return true;
        return false;
    }

    // neighbor point has a 1
    private boolean neighborOK(int row, int col, byte[][] ok, boolean debug)
    {
        if(row < 0 || row >= nRows || col < 0 || col >= nCols) return false;
        if(debug) System.out.println("    neighbor at row: " + row + " col: " + col + " ok: " + ok[row][col]);
        return ok[row][col] == 1;
    }

    static public void main(String[] args)
    {
		StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsHorpickWizard HorpickWizard = new StsHorpickWizard(actionManager);
        HorpickWizard.start();
    }
}
