package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTimeDepthDisplayable;
import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.MVC.View3d.StsView3d;
import com.Sts.SeismicAttributes.StsSeismicVelocityConstructor;
import com.Sts.Types.*;
import com.Sts.UI.Beans.StsBooleanFieldBean;
import com.Sts.UI.Beans.StsFieldBean;
import com.Sts.UI.Beans.StsFloatFieldBean;
import com.Sts.UI.ObjectPanel.StsObjectPanel;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsToolkit;

import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

/** Here are the scenarios for velocityModel construction:
 *
 *  1. Input data: velocityVolume.  Output data: adjustedVelocityVolume.
 *     Algorithm: inputVolume is copied to outputVolume; depth added as zDomain if not already;
 *     well deviation path times, markers times and surface depths are adjusted using this velocityVolume.
 */
public class StsSeismicVelocityModel extends StsSeismicBoundingBox implements ActionListener, StsTreeObjectI
{
	/** Top datum time coresponding to top datum depth z0 */
	public float timeDatum;
	/** Top datum depth coresponding to top datum time t0 */
	public float depthDatum;
	/** Time min for volume to be defined; corresponds to t0. */
	public float timeMin;
	/** Time max for volume to be defined; corresponds to max time for project. */
	public float timeMax;
	/** Time increment for volume to be defined; corresponds to time inc for project. */
	public float timeInc;
    /** A list of StsModelSurfaces which have arrays for time (pointsZ) and depth (adjPointsZ)
     *  note that this could be depth and adjusted depth when we are using a depth-migrated volume or surfaces
     *  interpreted on a depth-migrated volume.
     */
	public float markerFactor;
    public int gridType;
	/** Indicates well markers will be used in velocity model. */
	public boolean useWellControl;
	/* Indicates that interval velocities are constant for each layer */
	public boolean hasConstantLayerIntervalVelocities = false;

    StsObjectRefList surfaces;
    /** Instances of StsVelocityGrid which describe the interval velocity above each surface */
	StsObjectRefList intervalVelocityGrids;
    /** number of surfaces used in velocity model construction */
    int nSurfaces = 0;
    /** input velocity model; if not null this will be calibrated using interval velocity data to create the output velocityVolume. */
    public StsSeismicVolume inputVelocityVolume;
    /** Output average velocity volume */
    public StsSeismicVolume velocityVolume = null;
	/** Output instantaneous velocity volume */
    StsSeismicVolume intervalvelocityVolume = null;
	/** row offset of velocityVolume in velocity model grid. */
	public int volumeRowOffset = 0;
	/** col offset of velocityVolume in velocity model grid. */
	public int volumeColOffset = 0;
	/** factor multiplying velocities in input velocity volume. Since seismic data and picked horizons are two-way,
	 *  if input velocity volume is one-way, factor is 0.5 to reduce one-way to two-way velocities. */
	public double scaleMultiplier = 1.0;
	/** this velocity model has been initialized */
	public transient boolean initialized = false;

	transient boolean changed = true;
	transient StsProgressPanel panel;
    /** transient versions of StsVelocityGrid(s). When successfully constructed, they are converted to the persistent invervalVelocityGrids list */
    transient public StsVelocityGrid[] intervalVelocityGridArray;
    /** transient versions of StsModelSurface(s). When successfully constructed, they are converted to the persistent surfaces list */
    transient StsModelSurface[] surfacesArray;
    /** array of constant interval velocities if specified; these override interval velocities derived from marker/surfaces and adjust the
     *  time values for these points. */
    transient float[] oneWayConstantIntervalVelocities;
	/** the interval velocity below the bottom surface; used if non-zero; otherwise, interval velocity in layer above will be continued down. */
	transient float oneWayBottomConstantIntervalVelocity;

    transient float depthMax = -StsParameters.largeFloat;

    transient InterpolatedTrace interpolatedTrace = null;

    static public final String VELOCITY = "velocity";
	static float maxAllowedError = 0.02f;

	static public final boolean debug = false;
	static public final int debugVolumeRow = 10;
	static public final int debugVolumeCol = 10;

	static public boolean interpolationDebug = false;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

	static public StsObjectPanel velocityVolumeObjectPanel = null;

    public String getGroupname()
    {
        return groupNone;
    }

    public StsSeismicVelocityModel()
	{
	}

	private StsSurface constructSurface(StsSeismicVolume v, int i, float[][] z)
	{
		byte zDomainByte = StsParameters.getZDomainFromString("Depth");
		byte vUnits = currentModel.getProject().getDepthUnits();
		byte hUnits = currentModel.getProject().getXyUnits();
        String Name = new Integer(i).toString();
		StsColor c = new StsColor(1.f,1.f,1.f);
		byte type = (byte)1;
		return StsSurface.constructSurface(name, c, type, v.getNCols(), v.getNRows(), v.xOrigin, v.yOrigin,
		 v.getXInc(), v.getYInc(), v.getXMin(), v.getYMin(), 0.f, z, false, StsParameters.nullValue, zDomainByte, vUnits, hUnits, null);
	}

	public StsSeismicVelocityModel( float timeDatum, float deptDatum, float velMin, float velMax, double caleMultiplier, double newTimeInc, StsProgressPanel panel) throws StsException
	{
		super(false);

	    StsSeismicVolume[] volumes = (StsSeismicVolume[])currentModel.getCastObjectList(StsSeismicVolume.class);
		StsSeismicVolume v = volumes[0];


		nSurfaces=0;
		surfacesArray = null;
		try
		{
			this.panel = panel;
			panel.appendLine("Constructing velocity model.");
			{
				setName(VELOCITY);
				stsDirectory = currentModel.getProject().getDataFullDirString();
			}
			initializeNumsToRowColRange();

			currentModel.disableDisplay();

			if (nSurfaces <= 0 || !buildVelocityModel(surfacesArray))
			{
				panel.appendLine("Seismic velocity model construction FAILED.");
				return;
			}

			velocityVolume.initializeColorscale();
			addToModel();
            currentModel.getProject().checkSetZDomain(StsProject.TD_DEPTH, StsProject.TD_TIME);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "constructor", e);
		}
		finally
		{
			currentModel.enableDisplay();
        }
	}
	public StsSeismicVelocityModel(StsModelSurface[] timeSurfacesArray, StsSeismicVolume inputVelocityVolume,
								   float timeDatum, float depthDatum, float velMin, float velMax, double scaleMultiplier,

								   float inputTimeInc, float[] oneWayConstantIntervalVelocities, float oneWayConstantBottomIntervalVelocity, boolean useWellControl, float markerFactor,
								   int gridType, StsProgressPanel panel) throws StsException

	{
        super(false);
        try
        {
            this.panel = panel;
            panel.appendLine("Constructing velocity model.");

            this.surfacesArray = timeSurfacesArray;
            this.oneWayConstantIntervalVelocities = oneWayConstantIntervalVelocities;
			this.oneWayBottomConstantIntervalVelocity = oneWayConstantBottomIntervalVelocity;
			this.useWellControl = useWellControl;
            if (surfacesArray == null)
                nSurfaces = 0;
            else
            {
                // remove surfaces which have no markers or constant interval velocities
                removeNonVelocitySurfaces();
                nSurfaces = surfacesArray.length;
            }

            this.timeDatum = timeDatum;
            this.depthDatum = depthDatum;

			this.markerFactor = markerFactor;
            this.gridType = gridType;
            this.inputVelocityVolume = inputVelocityVolume;
            zDomain = StsParameters.TD_TIME_STRING;

            // any surfaces which aren't model surfaces are converted to model surfaces so they can be isVisible in depth
    //        checkBuildModelSurfaces();
            // PostStack3d will be same size laterally as model.gridDef set in horizon construction wizard.
            StsGridDefinition gridDefinition = currentModel.getGridDefinition();
            initializeToBoundingBox(gridDefinition);

            // To be consistent with seismic for simplicity, make the velocity volume orthogonal in t,
            // so that in depth we will have a vertically distorted grid.
            // Vertically, adjust zMin of velocityModel to be congruent with existing project vertical scale.
            // Ignore this inputTimeInc if timeInc has already been specified by a previously loaded volume of some kind.
            checkChangeZRange(inputTimeInc);

            panel.appendLine("Verifying data....");
            if (inputVelocityVolume != null)
            {
                panel.appendLine("Using " + inputVelocityVolume.getName() + " as input velocity volume.");

                // if we have no surfaces, surfaces have no markers, or we don't want to use well control,
                // then surfaces and wells can't be used for calibrating this volume so velocityVolume is
                // equal to the inputVolume.
            /* tjl
				if (nSurfaces == 0 && !surfacesHaveMarkers() && !useWellControl)
				{
					panel.appendLine("Velocity model will be the same as the input: " + inputVelocityVolume.getName() + ".");
					velocityVolume = inputVelocityVolume;
					currentModel.getProject().checkSetZDomain(StsParameters.TD_TIME_DEPTH, StsParameters.TD_DEPTH);
					// jbwidunno
					//checkBuildModelSurfaces();
				}
		    */
				if (nSurfaces == 0 || (surfacesHaveMarkers() && !useWellControl) || !hasConstantLayerIntervalVelocities)
                {
					if(nSurfaces == 0)
						panel.appendLine("No surfaces selected");

					if(!surfacesHaveMarkers())
						panel.appendLine("No surfaces with markers provided");
					else if(!useWellControl)
						panel.appendLine("Ignoring markers provided");

					if(!hasConstantLayerIntervalVelocities)
						panel.appendLine("No layer interval velocities provided");

                    useWellControl = false; 
    //				new StsMessage(currentModel.win3d, StsMessage.ERROR,
    //						   "No surfaces with markers available, so a velocity model is same as input volume.");
                    panel.appendLine("Velocity model will be the same as the input: " + inputVelocityVolume.getName() + ".");
                    velocityVolume = inputVelocityVolume;
                    // jbw velocityVolume.setType(StsParameters.SAMPLE_TYPE_VEL_AVG);
                    setName(inputVelocityVolume.stemname);
                    StsProject project = currentModel.getProject();

                    project.checkSetZDomain(StsProject.TD_DEPTH, StsProject.TD_TIME);

                    //adjustWells();
                    //adjustSurfaces(); // initial adjustment of surfaces using inputVelocityVolume
                    //addToModel();
                }
                else
                {
                    setName(inputVelocityVolume.stemname + "." + VELOCITY);
                }
                stsDirectory = inputVelocityVolume.stsDirectory;
                inputVelocityVolume.initialize(currentModel); // initialize for reading
                dataMin = (float)scaleMultiplier*velMin;
                dataMax = (float)scaleMultiplier*velMax;
                this.scaleMultiplier = scaleMultiplier;
            }
			else if (nSurfaces == 0 || (surfacesHaveMarkers() && !useWellControl)  && !hasConstantLayerIntervalVelocities)
            {

				{
					panel.appendLine("Unable to build model with specified data.");
					new StsMessage(currentModel.win3d, StsMessage.ERROR,
								   "No velocity volume nor time surfaces with markers available, so a velocity model cannot be built.");
					throw new StsException(StsException.WARNING, "Can't build velocity model: no data available.");
				}
            }
			else if(hasConstantLayerIntervalVelocities)
			{
				setName(VELOCITY);
				stsDirectory = currentModel.getProject().getDataFullDirString();
				//addToModel();
				//return;

			}
            else // we have no inputVelocityVolume but we have surfaces; velocityModel can be built using surfaces with markers or constantIntervalVelocities or a combination
            {
                setName(VELOCITY);
                stsDirectory = currentModel.getProject().getDataFullDirString();
            }
            initializeNumsToRowColRange();

            currentModel.disableDisplay();

			if (nSurfaces > 0 && !buildVelocityModel(surfacesArray))
			{
				panel.appendLine("Seismic velocity model construction FAILED.");
				return;
			}
		/*
            if (nSurfaces > 0 && !hasConstantLayerIntervalVelocities && !buildVelocityModel(surfacesArray))
            {
                panel.appendLine("Seismic velocity model construction FAILED.");
                return;
            }
        */
            velocityVolume.initializeColorscale();
            addToModel();
            currentModel.getProject().checkSetZDomain(StsProject.TD_DEPTH, StsProject.TD_TIME);
            panel.appendLine("Adjusting well and surface ties...");
            adjustWellsAndSurfaces(useWellControl);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructor", e);
        }
        finally
        {
            currentModel.enableDisplay();
        }
    }

    public void checkSetDepthRange()
    {
        int nRows = velocityVolume.nRows;
        int nCols = velocityVolume.nCols;
        int nSlices = velocityVolume.nSlices;
        float maxAvgVelocity = 0.0f;
        for(int row = 0; row < nRows; row++)
        {
            float[] values = velocityVolume.readRowPlaneFloatData(row);
            int i = nSlices-1;
            for(int col = 0; col < nCols; col++, i += nSlices)
                maxAvgVelocity = Math.max(maxAvgVelocity, values[i]);
        }
        currentModel.getProject().checkSetDepthRange(maxAvgVelocity, depthDatum, timeDatum);
    }

    public String getTypeAsString() { return "Velocity"; }

    public String getUnitsString() { return currentModel.getProject().getVelocityUnits(); }

    private void checkChangeZRange(float newTimeInc)
    {
        // PostStack3d is native in time units, so set z range values to time.
        if(inputVelocityVolume != null)
        {
            zMin = inputVelocityVolume.zMin;
            zMax = inputVelocityVolume.zMax;
            zInc = inputVelocityVolume.zInc;
            timeMin = zMin;
            timeMax = zMax;
            timeInc = zInc;
            nSlices = Math.round( (timeMax - timeMin) / timeInc) + 1;
        }
        else
        {
            currentModel.getProject().checkChangeCongruentZMin(timeDatum, StsProject.TD_TIME);
            currentModel.getProject().timeIncLocked = true;
            timeInc = currentModel.getProject().getTimeInc();
            timeMin = timeDatum;
            zMin = this.timeMin;
            zInc = this.timeInc;
        }
    }

    private boolean surfacesHaveMarkers()
	{
		Iterator surfaceIterator = currentModel.getObjectIterator(StsModelSurface.class);
		ArrayList surfaceList = new ArrayList();
		while (surfaceIterator.hasNext())
		{
			StsModelSurface surface = (StsModelSurface) surfaceIterator.next();
			if (surface.getZDomainOriginal() == StsParameters.TD_TIME && surface.getMarker() != null) return true;
		}
		return false;
	}

    public boolean hasModelSurface(StsModelSurface surface)
	{
        if(surfaces == null) return false;
        int nSurfaces = surfaces.getSize();
        for(int n = 0; n < nSurfaces; n++)
        {
            StsModelSurface modelSurface = (StsModelSurface)surfaces.getElement(n);
            if(modelSurface == surface) return true;
        }
		return false;
	}

    private boolean buildVelocityModel(StsModelSurface[] surfacesArray) throws StsException
	{
//        depthSurfacesArray = new StsModelSurface[nSurfaces];
		intervalVelocityGridArray = new StsVelocityGrid[nSurfaces];

		// classInitialize depthSurfaces to same size and geometry as timeSurfaces
		for (int n = 0; n < nSurfaces; n++)
			surfacesArray[n].initializeDepthFromTime();

//		panel.initializeIntervals(nSurfaces);
/*
        if(intervalVelocities != null)
        {
            for(int n = 0; n < intervalVelocities.length; n++)
            {
                float velocity = intervalVelocities[n];
                if(velocity < dataMin) dataMin = velocity;
                if(velocity > dataMax) dataMax = velocity;
            }
        }
 */
        // classInitialize intervalVelocityGrids with time-depth objects above and below the interval
        float constantVelocity = 0.0f;
        if(oneWayConstantIntervalVelocities != null)
            constantVelocity = oneWayConstantIntervalVelocities[0]/2;

        if (nSurfaces > 0)
		{
          panel.appendLine("Creating velocity model grid for surface " + 1 + " of " + nSurfaces);
          intervalVelocityGridArray[0] = StsVelocityGrid.construct(this, timeDatum, depthDatum, surfacesArray[0], inputVelocityVolume, constantVelocity, panel);

//        progress = 10;
//        panel.setProgress(progress);
		  for (int n = 1, nn = 2; n < nSurfaces; n++, nn++)
		  {
            panel.appendLine("Creating velocity model grid for surface " + nn + " of " + nSurfaces);
            constantVelocity = 0.0f;
            if(oneWayConstantIntervalVelocities != null)
                constantVelocity = oneWayConstantIntervalVelocities[n]/2;
			intervalVelocityGridArray[n] = StsVelocityGrid.construct(this, surfacesArray[n - 1], surfacesArray[n], inputVelocityVolume, constantVelocity, panel);
		  }

		  if (!constructIntervalVelocities())
          {
              panel.appendLine("Failed to construct interval velocities.");
              throw new StsException(StsException.WARNING, "Failed to build velocity model.");
          }
			// jbw not sure this is right place or mechanism, but minmax's aren't being set
		  dataMin = intervalVelocityGridArray[0].velMin;
		  dataMax = intervalVelocityGridArray[0].velMax;
		  for(int n = 1; n < intervalVelocityGridArray.length; n++)
			{
				dataMin = Math.min(dataMin, intervalVelocityGridArray[n].velMin);
				dataMax = Math.max(dataMax, intervalVelocityGridArray[n].velMax);
			}
		  float scaleFactor = 254 / (dataMax - dataMin);
		  scaleOffset = -dataMin * scaleFactor;

		}
		else
		{
			// jbw need to init

			StsToolkit.copySubToSuperclass(inputVelocityVolume, this, StsSeismicBoundingBox.class, StsBoundingBox.class,true);
			copyInputVelocityVolume();
		}
        if (inputVelocityVolume == null)
        {
//            computeIntervalVelocityRange();
            setVelocityVolumeTimeMax();
        }

        if(!computeVelocityVolume()) return false;
		if(!computeintervalVelocityVolume()) return false;
		/* if (debug) */ debugCheck();
 /*
        if (velocityVolume != null)
		{
			currentModel.project.checkSetDepthMax(StsSeismicVelocityConstructor.depthMax);
		}
*/
        return true;
    }

    private void setVelocityVolumeTimeMax()
    {
        StsModelSurface bottomSurface = (StsModelSurface)surfaces.getLast();
		timeMax = getCurrentProject().getTimeMax();
        // timeMax = StsMath.intervalRoundUp(bottomSurface.getTimeMax(), timeMin, timeInc);
        zMax = timeMax;
        nSlices = Math.round( (timeMax - timeMin) / timeInc) + 1;
	/*
		if(hasBottomConstantIntervalVelocity())
		{
			timeMax += timeInc;
			nSlices++;
		}
	*/
        dbFieldChanged("timeMax", timeMax);
    }
/*
	private boolean hasBottomConstantIntervalVelocity()
	{
		return btmIntervalVelocity != 0.0f;
	}
*/
    private void removeNonVelocitySurfaces()
    {
        float[] newIntervalVelocities = null;

        int nSurfaces = surfacesArray.length;
        StsModelSurface[] newSurfaces = new StsModelSurface[nSurfaces];
        if(oneWayConstantIntervalVelocities != null)
            newIntervalVelocities = new float[nSurfaces+1];
        int nn = 0;
		// if interval velocities are constant for each layer, we won't build explicit velocity volumes
		hasConstantLayerIntervalVelocities = inputVelocityVolume == null && !useWellControl;
        for(int n = 0; n < nSurfaces; n++)
        {
            StsModelSurface surface = surfacesArray[n];
			boolean surfaceHasMarkers = surface.getMarker() != null;
			boolean hasConstantIntervalVelocity = oneWayConstantIntervalVelocities != null && oneWayConstantIntervalVelocities[n] != 0.0f;
            if (surfaceHasMarkers || hasConstantIntervalVelocity)
			{
                newSurfaces[nn] = surface;
				if(surfaceHasMarkers)
					hasConstantLayerIntervalVelocities = false;
			    if(hasConstantIntervalVelocity)
                    newIntervalVelocities[nn] = oneWayConstantIntervalVelocities[n];
                nn++;
            }
        }
        if(nn < nSurfaces)
		{
				surfacesArray = (StsModelSurface[])StsMath.trimArray(newSurfaces, nn);
			if(oneWayConstantIntervalVelocities != null)
				oneWayConstantIntervalVelocities = (float[])StsMath.trimArray(newIntervalVelocities, nn);
		}
    }

    private void adjustWellsAndSurfaces(boolean useWellControl)
	{
		if (velocityVolume != null)
           velocityVolume.setupReadRowFloatBlocks();
		adjustNonMarkerSurfaces(useWellControl);
		adjustWells();
		updatePlannedWells();
		checkAdjustModel();
//		if (velocityVolume != null)
//           velocityVolume.deleteRowFloatBlocks();
	}

	private boolean computeVelocityVolume()
	{
		velocityVolume = StsSeismicVelocityConstructor.constructor(currentModel, this, panel);
        return velocityVolume != null;
//		dbFieldChanged("velocityVolume", velocityVolume);
	}
	private boolean computeintervalVelocityVolume()
	{
		intervalvelocityVolume = StsSeismicVelocityConstructor.intervalconstructor(currentModel, this, panel);
		return intervalvelocityVolume != null;
//		dbFieldChanged("velocityVolume", velocityVolume);
	}

	private void debugCheck()
	{
		velocityVolume.setupReadRowFloatBlocks();
		if (inputVelocityVolume != null) inputVelocityVolume.setupReadRowFloatBlocks();

		if (debug)
		{
			debugCheckVolume(inputVelocityVolume);
			debugCheckVolume(velocityVolume);
		}

		for (int n = 0; n < nSurfaces; n++)
		{
			StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(n);
			velocityGrid.debugCheck(this);
		}
//		velocityVolume.deleteRowFloatBlocks();
//		if (inputVelocityVolume != null) inputVelocityVolume.deleteRowFloatBlocks();
	}

	private void debugCheckVolume(StsSeismicVolume velocityVolume)
	{
		if (velocityVolume == null)return;

		int debugModelRow = debugVolumeRow;
		int debugModelCol = debugVolumeCol;

		String pathname = velocityVolume.fileMapRowFloatBlocks.pathname;
		System.out.println("StsSeismicVelocityModel.debugCheck() for velocityVolume " + pathname + " volume row " +
						   debugVolumeRow + " volume col " +
						   debugVolumeCol);

		float avgVelocity = velocityVolume.fileMapRowFloatBlocks.getRowCubeFloatValueDebug(debugVolumeRow, debugVolumeCol, 0);
		System.out.println("    vol top " + " slice 0 avg Velocity " + avgVelocity);

		for (int s = 0; s < nSurfaces; s++)
		{
			StsModelSurface botSurface = (StsModelSurface) surfaces.getElement(s);
			float botTime = botSurface.pointsZ[debugModelRow][debugModelCol];
            float botDepth = botSurface.adjPointsZ[debugModelRow][debugModelCol];
            float computedBotDepth = getDepthFromVelocityVolumeRowCol(debugModelRow, debugModelCol, botTime);
//            float computedByteBotDepth = getDepthFromVelocityByteVolumeRowCol(debugModelRow, debugModelCol, botTime);
            System.out.println("    surface[" + debugModelRow + "][" + debugModelCol + "]" + botSurface.getName() + " surface depth " +
                    botDepth + " volume float depth " + computedBotDepth);  // + " volume byte depth " + computedByteBotDepth
        /*
            float botSliceF = velocityVolume.getBoundedSliceCoor(botTime);
			int sliceAbove = (int) botSliceF;
			avgVelocity = velocityVolume.fileMapRowFloatBlocks.getFloatValueDebug(debugVolumeRow, debugVolumeCol,
				sliceAbove);
			System.out.println("    surface " + botSurface.getName() + " sliceAbove " + sliceAbove +
							   " avg Velocity " + avgVelocity);
			int sliceBelow = sliceAbove + 1;
			avgVelocity = velocityVolume.fileMapRowFloatBlocks.getFloatValueDebug(debugVolumeRow, debugVolumeCol,
				sliceBelow);
			System.out.println("    surface " + botSurface.getName() + " sliceBelow " + sliceBelow +
							   " avg Velocity " + avgVelocity);
        */
		}
        /*
        int botVolumeSlice = velocityVolume.nCroppedSlices - 1;
		avgVelocity = velocityVolume.fileMapRowFloatBlocks.getFloatValueDebug(debugVolumeRow, debugVolumeCol,
			botVolumeSlice);
		System.out.println("    vol bot " + " slice " + botVolumeSlice + " avg Velocity " + avgVelocity);
		*/
	}

	/** turn off seismic and turn on velocity model on cursors */
	private void setTextureDisplays()
	{
		( (StsSeismicVolumeClass) currentModel.getStsClass(StsSeismicVolume.class)).setIsVisibleOnCursor(false);
		( (StsVirtualVolumeClass) currentModel.getStsClass(StsVirtualVolume.class)).setIsVisibleOnCursor(false);
		getStsClass().setIsVisible(true);
	}

	/** called by db read operation and unconditionally initializes velocityModel */
	public boolean initialize(StsModel model)
	{
		if (!initialized) initializeFixedObjects();
		initializeIntervalVelocities(); // may have changed if velocityModel fields are changed
		initialized = true;
		return initialized;
	}

	/** Called once to classInitialize fixed objects */
	private void initializeFixedObjects()
	{
 /*
        initializeSpectrum();
		if (colorscale != null) colorscale.addActionListener(this);
//            allocateVolumes();
		getStsClass().setIsVisible(false);
*/
	}

	/** In the db read, a methodCmd may be called using an argument object
	 *  which hasn't yet been initialized.  Check and classInitialize if necessary.
	 *  This problem goes away if in dbRead the object.classInitialize(model) is
	 *  called in sequence with the methodCmds rather than the objects being
	 *  initialized in one go at the end of the transaction read.  Fix some day...
	 */
	public void checkInitialize()
	{
		if (!initialized) initialize(currentModel);
	}

	/** If a velocityVolume has not been created, the velocity model
	 *  is an interval velocity model.  Build the interval velocities from
	 *  bounding time and depth surfaces pairs.
	 */

	public boolean initializeIntervalVelocities()
	{
		if (intervalVelocityGrids == null)return true;
		int nVelocityGrids = intervalVelocityGrids.getSize();
		for (int n = 0; n < nVelocityGrids; n++)
		{
			StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(n);
			if (velocityGrid == null)continue;
			if (!velocityGrid.initializeIntervalVelocities())return false;
		}
		return true;
	}

	public void convertDepthOnlyWells()
	{
		StsObjectList wells = currentModel.getInstances(StsWell.class);
		int nWells = wells.getSize();
		for (int n = 0; n < nWells; n++)
		{
			StsWell well = (StsWell) wells.getElement(n);
			if(well.zDomainSupported != StsProject.TD_TIME_DEPTH) well.adjustFromVelocityModel(this); // this skips wells already in time
        }
	}

	public void close()
	{
		if (inputVelocityVolume != null) inputVelocityVolume.close();
		if (velocityVolume != null) velocityVolume.close();
        interpolatedTrace = null;
    }

	/** If the timeSurfaces are different, delete this velocityModel, build a new one
	 *  and return it.  If only tMin or zMin have been changed, change values and
	 *  return this current velocityModel.  If no changes, simply return it.
	 */
	public boolean checkChangeVelocityModel(StsModelSurface[] timeSurfacesArray, StsSeismicVolume inputVelocityVolume,
											float timeDatum, float depthDatum, float velMin, float velMax, double scaleMultiplier,
											float inputTimeInc, float[] intervalVelocities, boolean useWellControl,
                                            StsProgressPanel panel) throws StsException
	{
        this.panel = panel;
		changed = false;
        panel.appendLine("Velocity model exists, modifying as requested...");

	    if(this.inputVelocityVolume != inputVelocityVolume)
		{
			changed = true;
			StsSeismicVolume.close(inputVelocityVolume);
			this.inputVelocityVolume = inputVelocityVolume;
		}

		this.surfacesArray = timeSurfacesArray;
		if (surfacesArray == null)
			nSurfaces = 0;
		else
        {
            removeNonVelocitySurfaces();
            nSurfaces = surfacesArray.length;
        }

        panel.appendLine("Verifying data....");
		if (inputVelocityVolume != null)
		{
            panel.appendLine("Using " + inputVelocityVolume.getName() + " as input velocity volume.");
			if (nSurfaces == 0 || (surfacesHaveMarkers() && !useWellControl) || !hasConstantLayerIntervalVelocities)
			{
                if(nSurfaces == 0)
                    panel.appendLine("No surfaces selected");

                if(!surfacesHaveMarkers())
					panel.appendLine("No surfaces with markers provided");
				else if(!useWellControl)
					panel.appendLine("Ignoring markers provided");

				if(!hasConstantLayerIntervalVelocities)
					panel.appendLine("No layer interval velocities provided");

//				copyInputVelocityVolume();
                panel.appendLine("Velocity model will be the same as the input: " + inputVelocityVolume.getName() + ".");
			    velocityVolume = inputVelocityVolume;
				adjustNonMarkerSurfaces(useWellControl);
				setName(inputVelocityVolume.stemname);
				return false;
			}
			else
			{
				setName(inputVelocityVolume.stemname + "." + VELOCITY);
			}
			stsDirectory = inputVelocityVolume.stsDirectory;
			inputVelocityVolume.initialize(currentModel);
			dataMin = (float)scaleMultiplier*velMin;
			dataMax = (float)scaleMultiplier*velMax;
			this.scaleMultiplier = scaleMultiplier;
		}
		else if (nSurfaces == 0 || (surfacesHaveMarkers() && !useWellControl)  && !hasConstantLayerIntervalVelocities)
		{
            panel.appendLine("Unable to build model with specified data.");
			new StsMessage(currentModel.win3d, StsMessage.ERROR,
						   "No velocity volume nor time surfaces with markers available, so a velocity model cannot be built.");
			throw new StsException(StsException.WARNING, "Can't build velocity model: no data available.");
		}
		else if(hasConstantLayerIntervalVelocities)
		{
			return false;
		}
		else
		{
			setName(VELOCITY);
			stsDirectory = currentModel.getProject().getDataFullDirString();
		}

		if (velocityVolume != null)
		{
			velocityVolume.close();
			velocityVolume = null;
		}

		if (this.timeDatum != timeDatum)
		{
			changed = true;
			this.timeDatum = timeDatum;
		}
		if (this.depthDatum != depthDatum)
		{
			changed = true;
			this.depthDatum = depthDatum;
		}
		if (velMin != dataMin)
		{
			changed = true;
			dataMin = velMin;
		}
		if (velMax != dataMax)
		{
			changed = true;
			dataMax = velMax;
		}
		if (this.scaleMultiplier != scaleMultiplier)
		{
			changed = true;
			this.scaleMultiplier = scaleMultiplier;
		}
        if(changed) return true;

        panel.appendLine("Velocity model input has not changed.");
        boolean runIt = StsYesNoDialog.questionValue(currentModel.win3d,"Velocity model input has not changed. Run anyways?");
        if (!runIt)
            return false;
       else
       {
           panel.appendLine("Re-computing velocity model.");
           return true;
       }
	}

	private boolean copyInputVelocityVolume()
	{
		try
		{
			velocityVolume = (StsSeismicVolume)StsToolkit.deepCopy(inputVelocityVolume);
			velocityVolume.setName(inputVelocityVolume.stemname + "." + VELOCITY);
			velocityVolume.setIndex(-1);
			velocityVolume.addToModel();
			velocityVolume.initialize();
			setName(inputVelocityVolume.stemname);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSeismicVelocityModel.copyInputVelocityVolume() failed.",
				e, StsException.WARNING);
			return false;
		}
	}

	public boolean delete()
	{
		try
		{
			if(surfaces != null) surfaces.delete(); // deletes just the list
			//           depthSurfaces.deleteAll(); // deletes the list and all depthSurfaces
			if (intervalVelocityGrids != null) intervalVelocityGrids.deleteAll(); // deletes the list and all velocityGrids
			super.delete();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsVelocityModel.delete() failed.", e, StsException.WARNING);
			return false;
		}
	}
/*
	private void minChanged(float t0, float z0)
	{
		if (this.timeDatum != t0)
		{
			this.timeDatum = t0;
			fieldChanged("t0", t0);
		}
		if (this.depthDatum != z0)
		{
			this.depthDatum = z0;
			fieldChanged("z0", z0);
		}
	}
*/
	public boolean isChanged()
	{
		return changed;
	}

	private boolean constructIntervalVelocities()
	{
		try
		{
//			if (inputVelocityVolume != null) inputVelocityVolume.setupReadRowFloatBlocks();
			for (int n = 0; n < nSurfaces; n++)
			{

				intervalVelocityGridArray[n].constructIntervalVelocities(markerFactor, gridType);

			}
//			if (inputVelocityVolume != null) inputVelocityVolume.deleteRowFloatBlocks();
			if (velocityVolume != null) velocityVolume.close();
			velocityVolume = null;
			changed = false;
			panel.finished();

			this.surfaces = StsObjectRefList.convertObjectsToRefList(surfacesArray, "surfaces", this);
			for (int n = 0; n < nSurfaces; n++)
			{
				intervalVelocityGridArray[n].addToModel();
			}

			this.intervalVelocityGrids = StsObjectRefList.convertObjectsToRefList(intervalVelocityGridArray,
				"intervalVelocityGrids", this);

			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsVelocityModel.constructModel() failed.",
										 e, StsException.WARNING);
			delete();
			return false;
		}
	}

	private void adjustSurfaces()
	{
        Iterator surfaceIterator = currentModel.getObjectIterator(StsModelSurface.class);
		ArrayList surfaceList = new ArrayList();
		while (surfaceIterator.hasNext())
		{
            StsModelSurface surface = (StsModelSurface) surfaceIterator.next();
            adjustSurface(surface);
        }
    }

    /** convert any surfaces which aren't currently model surfaces to model surfaces */
    private void checkBuildModelSurfaces()
    {
        StsObject[] nonModelSurfaces = currentModel.getObjectList(StsSurface.class);
        int nNonModelSurfaces = nonModelSurfaces.length;
        StsObject[] modelSurfaces = currentModel.getObjectList(StsModelSurface.class);
        StsGridDefinition gridDef = currentModel.getProject().getGridDefinition();
        StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass)currentModel.getStsClass(StsModelSurface.class);
        StsSurface[] buildSurfaces = new StsSurface[nNonModelSurfaces];
        int nBuildSurfaces = 0;
        for(int n = nNonModelSurfaces-1; n >= 0; n--)
        {
            StsSurface nonModelSurface = (StsSurface)nonModelSurfaces[n];
            if(!surfaceHasModelSurface(nonModelSurface, modelSurfaces))
                buildSurfaces[nBuildSurfaces++] = nonModelSurface;
        }
        if(nBuildSurfaces == 0) return;
        buildSurfaces = (StsSurface[])StsMath.trimArray(buildSurfaces, nBuildSurfaces);
        StsGridDefinition gridDefinition = currentModel.getGridDefinition();
        if(gridDefinition == null)
        {
            StsSurfaceClass surfaceClass = (StsSurfaceClass)currentModel.getStsClass(StsSurface.class);
            gridDefinition = surfaceClass.checkComputeUnionGrid(buildSurfaces);

        }
        if(gridDefinition == null) return;

        currentModel.getProject().setGridDefinition(gridDefinition);

        for(int n = 0; n < nBuildSurfaces; n++)
            StsModelSurface.constructModelSurfaceFromSurface(buildSurfaces[n], modelSurfaceClass, gridDefinition);
    }

    private void adjustNonMarkerSurfaces(boolean useWellControl)
	{
        StsObject[] modelSurfaces = currentModel.getObjectList(StsModelSurface.class);
        for(int n = 0; n < modelSurfaces.length; n++)
		{
			StsModelSurface modelSurface = (StsModelSurface)modelSurfaces[n];
            byte zDomainOriginal = modelSurface.getZDomainOriginal();
            StsMarker marker = modelSurface.getMarker();
            if ((zDomainOriginal == StsParameters.TD_TIME) && (marker == null || !useWellControl))
			{
                adjustSurface(modelSurface);
            }
        }
    }

    private boolean surfaceHasModelSurface(StsSurface nonModelSurface, StsObject[] modelSurfaces)
    {
        for(int n = 0; n < modelSurfaces.length; n++)
            if(((StsModelSurface)modelSurfaces[n]).originalSurface == nonModelSurface)
                return true;
        return false;
    }

    private void adjustSurface(StsModelSurface surface)
    {
        if(surfaces != null && surfaces.contains(surface)) return; // already adjusted
        float[][] times = surface.getPointsZ();
        float[][] depths = surface.getCreateAdjPointsZ();
        int nRows = surface.getNRows();
        int nCols = surface.getNCols();
        float y = surface.yMin;
        float yInc = surface.yInc;
        float xInc = surface.xInc;
        float depthMin = largeFloat;
        float depthMax = -largeFloat;
        for (int row = 0; row < nRows; row++, y += yInc)
        {
            float x = surface.xMin;
            for (int col = 0; col < nCols; col++, x += xInc)
            {
                float depth = (float) getZ(x, y, times[row][col]);
                depths[row][col] = depth;
                depthMin = Math.min(depthMin, depth);
                depthMax = Math.max(depthMax, depth);
            }
        }
        currentModel.getProject().checkChangeZRange(depthMin, depthMax, StsProject.TD_DEPTH);
    }

    /* jbw use well td-curves (and auxiliar vel funcs) to build vel model */
	private StsSeismicVolume doWellTimeDepth()
	{
		StsObjectList wells = currentModel.getInstances(StsWell.class);
		int nWells = wells.getSize();
		for (int n = 0; n < nWells; n++)
		{
			StsWell well = (StsWell) wells.getElement(n);
			well.addTimeDepthToVels(this);
		}
		return null;
    }
	private void adjustWells()
	{
		StsObjectList wells = currentModel.getInstances(StsWell.class);
		int nWells = wells.getSize();
		for (int n = 0; n < nWells; n++)
		{
			StsWell well = (StsWell) wells.getElement(n);
			well.checkAdjustWellTimes(this);
		}
	}

	private void updatePlannedWells()
	{
		StsObjectList wellPlanSets = currentModel.getInstances(StsWellPlanSet.class);
		int nWellPlanSets = wellPlanSets.getSize();
		for (int n = 0; n < nWellPlanSets; n++)
		{
			StsWellPlanSet wellPlanSet = (StsWellPlanSet)wellPlanSets.getElement(n);
			StsWellPlan plannedWell = wellPlanSet.copyLastWellPlan();
			if (plannedWell == null)continue;
			// if planned well is in depth, convert to time with velocity model
			// otherwise wells are planned in time, so update depth values
			plannedWell.adjustDepthPoints(this);
			plannedWell.setIndex(-1);
			plannedWell.addToModel();
		}
	}

	/** Add or adjust depth or time values to points when velocityModel is built or changed */
	private void checkAdjustModel()
	{
		StsArrayList classList = currentModel.getTimeDepthDisplayableClasses();
        Iterator classIterator = classList.iterator();
        while(classIterator.hasNext())
        {
            StsClass tdDisplayableClass = (StsClass)classIterator.next();
            StsObjectList objectList = tdDisplayableClass.getStsObjectList();
            StsObject[] objects = objectList.getObjectList();
			for(StsObject object : objects)
				((StsTimeDepthDisplayable)object).adjustTimeOrDepth(this);
		}
	}

	private boolean reconstructIntervalVelocities(float markerFactor, int gridType)
	{
		try
		{
			currentModel.disableDisplay();

			if (inputVelocityVolume != null) inputVelocityVolume.setupReadRowFloatBlocks();
			for (int n = 0; n < nSurfaces; n++)
			{
				StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(n);

				velocityGrid.reconstructIntervalVelocities(panel, markerFactor, gridType);
			}
			close();
			reconstructVelocityVolume();
			changed = false;
			currentModel.enableDisplay();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsVelocityModel.constructIntervalVelocities() failed.",
										 e, StsException.WARNING);
			delete();
			return false;
		}
	}

	private void computeIntervalVelocityVolume()
	{

	}

	private boolean reconstructVelocityVolume()
	{
//        classInitialize(currentModel);
		if(inputVelocityVolume == null)
         {
//            computeIntervalVelocityRange();
            setVelocityVolumeTimeMax();
        }

        if(!computeVelocityVolume()) return false;
		convertDepthOnlyWells(); // convert depth-only wells to time
        return true;
    }

	public boolean computeIntervalVelocityRange()
	{
		if (intervalVelocityGrids == null)return false;
//		dataMin = largeFloat;
//		dataMax = -largeFloat;
		for (int n = 0; n < nSurfaces; n++)
		{
			float[] velRange = ( (StsVelocityGrid) intervalVelocityGrids.getElement(n)).getVelocityRange();
			dataMin = Math.min(dataMin, velRange[0]);
			dataMax = Math.max(dataMax, velRange[1]);
		}
		if(velocityVolume != null && velocityVolume.colorscale != null)
            velocityVolume.colorscale.setEditRange(dataMin, dataMax);
		return true;
	}
/*
    public void adjustVelocityRange(float intervalVelocity)
    {
        if(intervalVelocity < dataMin)
            dataMin = intervalVelocity;
        if(intervalVelocity > dataMax)
            dataMax = intervalVelocity;
    }
*/
    public byte[] getTraceVelocityPlane(int dir, int nPlane, StsSeismicVelocityModel seismicVelocityModel)
	{
		float[][] timeGridlines = new float[nSurfaces][];
		double[][] intVelLines = new double[nSurfaces][];
		byte[] traceVelocityPlane;
		int sliceFirst, sliceLast;

		if (dir == StsCursor3d.XDIR)
		{
			int col = nPlane;
			for (int n = 0; n < nSurfaces; n++)
			{
				timeGridlines[n] = ( (StsModelSurface) surfaces.getElement(n)).getColGridlineValues(col);
				intVelLines[n] = ( (StsVelocityGrid) intervalVelocityGrids.getElement(n)).getColGridlineValues(col);
			}
			traceVelocityPlane = new byte[nRows * nSlices];
			int b = 0;
			for (int row = 0; row < nRows; row++)
			{
				int slice = 0;
				byte velocityByte = 0;
				for (int n = 0; n < nSurfaces; n++)
				{
					velocityByte = seismicVelocityModel.getByteValue(intVelLines[n][row]);
					float z = timeGridlines[n][row];
					sliceLast = (int) ( (z - zMin) / zInc);
					for (; slice <= sliceLast; slice++)
					{
						traceVelocityPlane[b++] = velocityByte;
					}
				}
				for (; slice < nSlices; slice++)
				{
					traceVelocityPlane[b++] = velocityByte;
				}
			}
			return traceVelocityPlane;
		}
		if (dir == StsCursor3d.YDIR)
		{
			int row = nPlane;
			for (int n = 0; n < nSurfaces; n++)
			{
				timeGridlines[n] = ( (StsModelSurface) surfaces.getElement(n)).getRowGridlineValues(row);
				intVelLines[n] = ( (StsVelocityGrid) intervalVelocityGrids.getElement(n)).getRowGridlineValues(row);
			}
			traceVelocityPlane = new byte[nCols * nSlices];
			int b = 0;
			for (int col = 0; col < nCols; col++)
			{
				int slice = 0;
				byte velocityByte = 0;
				for (int n = 0; n < nSurfaces; n++)
				{
					velocityByte = seismicVelocityModel.getByteValue(intVelLines[n][col]);
					float z = timeGridlines[n][col];
					sliceLast = (int) ( (z - zMin) / zInc);
					for (; slice <= sliceLast; slice++)
					{
						traceVelocityPlane[b++] = velocityByte;
					}
				}
				for (; slice < nSlices; slice++)
				{
					traceVelocityPlane[b++] = velocityByte;
				}
			}
			return traceVelocityPlane;
		}
		else
		{
			return null;
		}
	}

	public double getT(float[] xyzmt) throws StsException
	{
		checkInitialize();
		if (inputVelocityVolume == null)
			return getTimeFromIntervalVelocity(xyzmt);
		else
			return getTimeFromVelocityVolume(xyzmt);
	}

	public double getT(float x, float y, float z, float estimatedTime)
	{
		checkInitialize();
		if (inputVelocityVolume == null)
			return getTimeFromIntervalVelocity(x, y, z);
		else
			return getTimeFromVelocityVolume(x, y, z, estimatedTime);
	}

    public boolean adjustTimeOrDepthPoint(StsPoint point, boolean isOriginalDepth)
    {
        if(!isOriginalDepth) // point is in time and we need to set correct depth
        {
            float t = point.getT();
            float z = (float) getZ(point.getX(), point.getY(), t);
            // if this is original point, it has X, Y, T.  We need to extend vector to be
            // X, Y, Z, M, T. We then need to move t to new location and set new value of z.
            if(point.checkExtendVector(5))
                point.setT(t);
            point.setZ(z);
            return true;
        }
        else // point is in depth and we need to set the correct time
        {
            try
            {
                // compute t from velocity model given x, y, z
                float t = (float) getT(point.v);
                // if this is original point, it has X, Y, Z.  We need to extend vector to be
                // X, Y, Z, M, T. We then need to set t in new location.
                point.checkExtendVector(5);
                point.setT(t);
                return true;
            }
            catch(Exception e)
            {
                return false;
            }
        }
    }

    public float[] getTimesForDepthRange(float x, float y, float zMin, float zInc, int nValues)
    {
 		checkInitialize();
		float[] times = getTimesFromIntervalVelocity(x, y, zMin, zInc, nValues);
        return times;
//        if (inputVelocityVolume == null) return times;
//	    return getTimesFromVelocityVolume(x, y, zMin, zInc, nValues, times);
    }
    private double getTimeFromIntervalVelocity(float[] xyz)
	{
        return getTimeFromIntervalVelocity(xyz[0], xyz[1], xyz[2]);
    }

    public double getTimeFromIntervalVelocity(float x, float y, float z)
    {
        float rowF = this.getBoundedRowCoor(y);
		float colF = this.getBoundedColCoor(x);
		try
		{
			StsModelSurface topSurface = null; // timeSurface above z (if any)
			StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(0); // grid of interval velocities in interval containing z
			float topTime = timeDatum; // time at surface above z or at top of volume
			float topDepth = depthDatum; // depth at surface above z or at top of volume
			//           StsModelSurface botDepthSurface = (StsModelSurface) depthSurfaces.getElement(0); // depth surface below z
			StsModelSurface botSurface = (StsModelSurface) surfaces.getElement(0); // time surface below z
			float botDepth = botSurface.interpolateDepthNoNulls(rowF, colF); // depth at surface below z

			if (z > botDepth)
			{
				for (int n = 1; n < nSurfaces; n++)
				{
					topDepth = botDepth;
					topSurface = botSurface;
					topTime = topSurface.interpolateBilinearNoNulls(rowF, colF);
					botSurface = (StsModelSurface) surfaces.getElement(n);
					botSurface = (StsModelSurface) surfaces.getElement(n);
					botDepth = botSurface.interpolateDepthNoNulls(rowF, colF);
					velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(n);
					if (z <= botDepth)break;
					if (n == nSurfaces - 1)
					{
						topDepth = botDepth;
						topSurface = botSurface;
						topTime = topSurface.interpolateBilinearNoNulls(rowF, colF);
					}
				}
			}
			double intervalVelocity = velocityGrid.getVelocity(rowF, colF);
			//System.out.println("Interval Velocity=" + intervalVelocity);
			if (intervalVelocity == nullValue)return nullValue;
			return topTime + (z - topDepth) / intervalVelocity;
		}
		catch (Exception e)
		{
			StsException.systemError("StsSeismicVelocityModel.getTimeFromIntervalVelocity() failed to get time for row " + rowF + " col " + colF + " depth " + z);
			return nullValue;
		}
	}

    public float[] getTimesFromIntervalVelocity(float x, float y, float zMin, float zInc, int nValues)
    {
        float rowF = this.getBoundedRowCoor(y);
		float colF = this.getBoundedColCoor(x);
        float z = zMin;
        try
		{
			StsModelSurface topSurface = null; // timeSurface above z (if any)
			StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(0); // grid of interval velocities in interval containing z
			float topTime = timeDatum; // time at surface above z or at top of volume
			float topDepth = depthDatum; // depth at surface above z or at top of volume
			//           StsModelSurface botDepthSurface = (StsModelSurface) depthSurfaces.getElement(0); // depth surface below z
			StsModelSurface botSurface = (StsModelSurface) surfaces.getElement(0); // time surface below z
			float botDepth = botSurface.interpolateDepthNoNulls(rowF, colF); // depth at surface below z

            float[] times = new float[nValues];
            for(int i = 0; i < nValues; i++, z += zInc)
            {
                if (z > botDepth)
                {
                    for (int n = 1; n < nSurfaces; n++)
                    {
                        topDepth = botDepth;
                        topSurface = botSurface;
                        topTime = topSurface.interpolateBilinearNoNulls(rowF, colF);
                        botSurface = (StsModelSurface) surfaces.getElement(n);
                        botSurface = (StsModelSurface) surfaces.getElement(n);
                        botDepth = botSurface.interpolateDepthNoNulls(rowF, colF);
                        velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(n);
                        if (z <= botDepth)break;
                        if (n == nSurfaces - 1)
                        {
                            topDepth = botDepth;
                            topSurface = botSurface;
                            topTime = topSurface.interpolateBilinearNoNulls(rowF, colF);
                        }
                    }
                }
                double intervalVelocity = velocityGrid.getVelocity(rowF, colF);
                times[i] = (float)(topTime + (z - topDepth) / intervalVelocity);
            }
            return times;
        }
		catch (Exception e)
		{
			StsException.systemError("StsSeismicVelocityModel.getTimeFromIntervalVelocity() failed to get time for row " + rowF + " col " + colF + " depth " + z);
			return null;
		}
	}

	private float getTimeFromVelocityVolume(float[] xyzmt) throws StsException
	{
		if (xyzmt.length < 5)
		{
			throw new StsException(StsException.WARNING, "point.length < 5 for a time point.");
		}
		float x = xyzmt[0];
		float y = xyzmt[1];
		float z = xyzmt[2];
		float estimatedTime = xyzmt[4];
		return getTimeFromVelocityVolume(x, y, z, estimatedTime);
	}

	public float getTimeFromVelocityVolume(float x, float y, float z, float estimatedTime)
	{
		float v0, t0, z0, v1, t1 = 0.0f, z1 = 0.0f;
		float rowF = -1, colF = -1;
		try
		{
			rowF = velocityVolume.getBoundedRowCoor(y);
			colF = velocityVolume.getBoundedColCoor(x);

			if (interpolationDebug)
			{
				System.out.println("getTimeFromVelocityVolume mainDebug. row " + rowF + " col " + colF);
			}
			initializeVelocityTrace(velocityVolume, rowF, colF);
			int nMaxSlice = velocityVolume.nSlices - 1;
			float tMin = velocityVolume.zMin;
			float tMax = velocityVolume.zMax;
			float tInc = velocityVolume.zInc;
			v0 = getVelocityTraceValue(0);
			float zTop = v0 * (tMin - timeDatum) + depthDatum;
			if (z <= zTop)
			{
				return (z - depthDatum) / v0 + timeDatum;
			}
			v0 = getVelocityTraceValue(nMaxSlice);
			float zBot = v0 * (tMax - timeDatum) + depthDatum;
			if (z >= zBot)
			{
				return (z - depthDatum) / v0 + timeDatum;
			}
			float nSliceF = velocityVolume.getBoundedSliceCoor(estimatedTime);
			int nSlice = (int) nSliceF;
			if (nSlice >= nMaxSlice) nSlice = nMaxSlice - 1;

			v0 = getVelocityTraceValue(nSlice);
			t0 = tMin + nSlice * tInc;
			z0 = v0 * (t0 - timeDatum) + depthDatum;
			float f;
			if (z0 <= z)
			{
				z1 = z0;
				t1 = t0;
				v1 = v0;
				while (nSlice < nMaxSlice)
				{
					nSlice++;
					z0 = z1;
					t0 = t1;
					v0 = v1;
					v1 = getVelocityTraceValue(nSlice);
					t1 = t0 + tInc;
					z1 = v1 * (t1 - timeDatum) + depthDatum;
					if (z1 >= z)
					{
						f = (z - z0) / (z1 - z0);
						float t = t0 + f * tInc;
						if (interpolationDebug)
						{
							float v = v0 + f * (v1 - v0);
							System.out.println("getTimeFromVelocityVolume mainDebug. volumeRowF " + rowF + " volumeColF " +
											   colF + " volumeSliceF " + nSlice + f + " time " + t + " avg velocity " +
											   v);
						}
						return t;
					}
				}
				f = (z - z0) / (z1 - z0);
				return t0 + f * tInc;
			}
			else // z0 > z
			{
				if (nSlice == 0)
				{
					z1 = z0;
					t1 = t0;
					v0 = getVelocityTraceValue(0);
					t0 = t1 - tInc;
					z0 = v0 * (t0 - timeDatum) + depthDatum;
					f = (z - z0) / (z1 - z0);
					float t = t0 + f * tInc;
					if (interpolationDebug)
					{
						System.out.println("getTimeFromVelocityVolume mainDebug. volumeRowF " + rowF + " volumeColF " +
										   colF + " volumeSliceF " + nSlice + f + " time " + t + " avg velocity " + v0);
					}
					return t;
				}
				while (nSlice > 0)
				{
					nSlice--;
					z1 = z0;
					t1 = t0;
					v1 = v0;
					v0 = getVelocityTraceValue(nSlice);
					t0 = t1 - tInc;
					z0 = v0 * (t0 - timeDatum) + depthDatum;
					if (z0 <= z)
					{
						f = (z - z0) / (z1 - z0);
						float t = t0 + f * tInc;
						if (interpolationDebug)
						{
							float v = v0 + f * (v1 - v0);
							System.out.println("getTimeFromVelocityVolume mainDebug. volumeRowF " + rowF + " volumeColF " +
											   colF + " volumeSliceF " + nSlice + f + " time " + t + " avg velocity " +
											   v);
						}
						return t;
					}
				}
				f = (z - z0) / (z1 - z0);
				float t = t0 + f * tInc;
				if (interpolationDebug)
				{
					System.out.println("getTimeFromVelocityVolume mainDebug. volumeRowF " + rowF + " volumeColF " + colF +
									   " volumeSliceF " + nSlice + f + " time " + t + " avg velocity " + v0);
				}
				return t;
			}
		}
		catch (Exception e)
		{
			StsException.systemError("StsSeismicVolume.getTimeFromVelocityVolume() failed for row " + rowF + " col " +
									 colF);
			return estimatedTime;
		}
	}

    public double getZ(StsPoint point)
    {
        if(point.v.length < 5)
        {
            StsException.systemError(this, "getZ", "Can't get time from point: " + point.toString());
            return point.getZ();
        }
        return getZ(point.getX(), point.getY(), point.getT());
    }

    public double getZ(float x, float y, float t)
	{
		double value;
		if (inputVelocityVolume == null)
			return getDepthFromIntervalVelocity(x, y, t);
		else
		{
			value = getDepthFromVelocityVolume(x, y, t);
			if (value != nullValue)return value;
			return getDepthFromIntervalVelocity(x, y, t);
		}
	}

    public double getZFromRowCol(float rowF, float colF, float t)
	{
		double value;
		if (inputVelocityVolume == null)
			return getDepthFromIntervalVelocityRowCol(rowF, colF, t);
		else
		{
			value = getDepthFromVelocityVolumeRowCol(rowF, colF, t);
			if (value != nullValue)return value;
			return getDepthFromIntervalVelocityRowCol(rowF, colF, t);
		}
	}
	/*
			public double getZ(float x, float y, float t)
			{
				double value;
				if(velocityVolume == null || !velocityVolume.isInsideXY(x, y))
					return getDepthFromIntervalVelocity(x, y, t);
				else
				{
					value = getDepthFromVelocityVolume(x, y, t);
					if(value != nullValue) return value;
					return getDepthFromIntervalVelocity(x, y, t);
				}
			}
	 */
	public double getDepthFromIntervalVelocity(float x, float y, float t)
	{
		float rowF = this.getBoundedRowCoor(y);
		float colF = this.getBoundedColCoor(x);
        return getDepthFromIntervalVelocityRowCol(rowF, colF, t);
    }

	public double getDepthFromIntervalVelocityRowCol(float rowF, float colF, float t)
	{
		float topTime = timeDatum; // time at surface above z or at top of volume
		float topDepth = depthDatum; // depth at surface above z or at top of volume
		StsModelSurface botSurface = (StsModelSurface) surfaces.getElement(0);
		StsVelocityGrid velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(0);
		float botTime = botSurface.interpolateBilinearNoNulls(rowF, colF);
		float botDepth = botSurface.interpolateDepthNoNulls(rowF, colF);
		int nSurface = 0;

		while(t <= botTime && nSurface < nSurfaces)
		{
			topTime = botTime;
			topDepth = botDepth;
			botSurface = (StsModelSurface) surfaces.getElement(nSurface);
			botTime = botSurface.interpolateBilinearNoNulls(rowF, colF);
			botDepth = botSurface.interpolateDepthNoNulls(rowF, colF);
			velocityGrid = (StsVelocityGrid) intervalVelocityGrids.getElement(nSurface);
			nSurface++;
		}
		double intervalVelocity = nullValue;
		if (nSurface < nSurfaces || oneWayBottomConstantIntervalVelocity == 0.0f) // null if z is above first surface
			intervalVelocity = velocityGrid.getVelocity(rowF, colF);
		else
			intervalVelocity = oneWayBottomConstantIntervalVelocity;
		return topDepth + (t - topTime) * intervalVelocity;
	}

	private float getDepthFromVelocityVolume(float x, float y, float t)
	{
        float rowF = velocityVolume.getBoundedRowCoor(y);
        float colF = velocityVolume.getBoundedColCoor(x);
        return getDepthFromVelocityVolumeRowCol(rowF, colF, t);
	}

    private float getDepthFromVelocityVolumeRowCol(float rowF, float colF, float t)
	{
		try
		{
			initializeVelocityTrace(velocityVolume, rowF, colF);
			float nSliceF = velocityVolume.getBoundedSliceCoor(t);
			int nSlice = (int) nSliceF;
			float v = getVelocityTraceValue(nSlice);
			return v*(t - timeDatum) + depthDatum;
		}
		catch (Exception e)
		{
			StsException.systemError("getDepthFromVelocityVolume.getValue() failed for row " + rowF + " col " + colF);
			return nullValue;
		}
	}

    public float getDepthFromVelocityByteVolumeRowCol(int row, int col, float t)
	{
		try
		{
			float v = velocityVolume.getScaledByteValue(row, col, t);
			return v*(t - timeDatum) + depthDatum;
		}
		catch (Exception e)
		{
			StsException.systemError("getDepthFromVelocityByteVolumeRowCol.getValue() failed for row " + row + " col " + col);
			return nullValue;
		}
	}

    public boolean initializeVelocityTrace(StsSeismicVolume velocityVolume, float rowF, float colF)
    {
        if(velocityVolume == null) return false;
        if (interpolatedTrace == null || interpolatedTrace.velocityVolume != velocityVolume)
            interpolatedTrace = new InterpolatedTrace(velocityVolume);
        return interpolatedTrace.initialize(rowF, colF);
    }

    public float getVelocityTraceValue(int nSlice)
    {
        return interpolatedTrace.getValue(nSlice);
    }

    public boolean getVelocityTraceValues(StsSeismicVolume velocityVolume, float rowF, float colF, float[] values)
    {
        if(!initializeVelocityTrace(velocityVolume, rowF, colF)) return false;
        return interpolatedTrace.getValues(rowF, colF, values);
    }


    public float[] getDepthTraceFromIntervalVelocities(float timeMin, float timeInc, int nSlices, float x, float y)
    {
        try
        {
            int row = getNearestBoundedRowCoor(y);
            int col = getNearestBoundedColCoor(x);
            int nSurfaces = getNSurfaces();
            float[] traceTimeValues = getSurfaceTimes(row, col);
            float[] traceIntervalVelocities = getIntervalVelocities(row, col);
            float[] traceDepthValues = getSurfaceDepths(row, col);
            int nBotSurface = 0;
            float topTime = timeDatum;
            float botTime = traceTimeValues[0];
            float topDepth = depthDatum;
            float botDepth = traceDepthValues[0];
            float intervalVelocity = traceIntervalVelocities[0];
            float[] depths = new float[nSlices];
            float t = timeMin;
            for (int n = 0; n < nSlices; n++, t += timeInc)
            {
                while(t > botTime && nBotSurface < nSurfaces-1)
                {
                    topTime = botTime;
                    topDepth = botDepth;
                    nBotSurface++;
                    botTime = traceTimeValues[nBotSurface];
                    botDepth = traceDepthValues[nBotSurface];
                    intervalVelocity = traceIntervalVelocities[nBotSurface];
                }
                depths[n] = topDepth + (t - topTime)*intervalVelocity;
            }
            return depths;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getDepthTraceFromIntervalVelocities", e);
            return null;
        }
    }

    public class InterpolatedTrace
    {
        StsSeismicVolume velocityVolume;
        int row = -1, col = -1;
        float dx, dy;
        int index;
        float[] trace00, trace01, trace10, trace11;
        int type = 0;

        InterpolatedTrace(StsSeismicVolume velocityVolume)
        {
            this.velocityVolume = velocityVolume;
            velocityVolume.setupReadRowFloatBlocks();
            trace00 = new float[nSlices];
            trace01 = new float[nSlices];
            trace10 = new float[nSlices];
            trace11 = new float[nSlices];
        }

        boolean initialize(float rowF, float colF)
        {
            rowF = StsMath.minMax(rowF, 0.0f, velocityVolume.nRows-1);
            row = (int)rowF;
            dy = rowF - row;
            if(dy < 0.01f)
                dy = 0.0f;
            else if(dy > 0.99f)
            {
                dy = 0.0f;
                row++;
            }

            colF = StsMath.minMax(colF, 0.0f, velocityVolume.nCols-1);
            col = (int)colF;
            dx = colF - col;
            if(dx < 0.01f)
                dx = 0.0f;
            else if(dx > 0.99f)
            {
                dx = 0.0f;
                col++;
            }

            if (dy == 0.0f)
            {
                if (dx == 0.0f)
                {
                    if(!getTraceValues(row, col, 0, nSlices - 1, YDIR, false, trace00)) return false;
                    type = 0;
                }
                else
                { // dy == 0.0f, point is along bottom row
                    if(!getTraceValues(row, col, 0, nSlices - 1, YDIR, false, trace00)) return false;
                    if(!getTraceValues(row, col + 1, 0, nSlices - 1, YDIR, false, trace01)) return false;
                    type = 1;
                }
            }
            else
            { // dy != 0.0f
                if (dx == 0.0f)
                { // point is along left column
                    if(!getTraceValues(row, col, 0, nSlices - 1, YDIR, false, trace00)) return false;
                    if(!getTraceValues(row + 1, col, 0, nSlices - 1, YDIR, false, trace10)) return false;

                    type = 3;
                }
                else
                { // point is in the middle
                    if(!getTraceValues(row, col, 0, nSlices - 1, YDIR, false, trace00)) return false;
                    if(!getTraceValues(row, col + 1, 0, nSlices - 1, YDIR, false, trace01)) return false;
                    if(!getTraceValues(row + 1, col, 0, nSlices - 1, YDIR, false, trace10)) return false;
                    if(!getTraceValues(row + 1, col + 1, 0, nSlices - 1, YDIR, false, trace11)) return false;
                    type = 4;
                }
            }
            return true;
        }

        boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, boolean useByteCubes, float[] floatData)
	    {
            if(!velocityVolume.getTraceValues(row, col, sliceMin, sliceMax, dir, useByteCubes, floatData)) return false;
            return floatData[0] != StsParameters.nullValue;
        }

        boolean getValues(float rowF, float colF, float[] values)
        {
            if(!initialize(rowF, colF)) {
				System.out.println("error getting row "+row+" "+col);
				System.out.println("error getting row "+nRows+" "+nCols);
				return false;
			}
            if (type == 0)
            { // point is in lower-left corner
                System.arraycopy(trace00, 0, values, 0, nSlices);
            }
            else if (type == 1)
            { // point is along bottom row
                return StsMath.interpolate(trace00, trace01, nSlices, dx, values);
            }
            else if (type == 3)
            { // point is along left column
                return StsMath.interpolate(trace00, trace10, nSlices, dy, values);
            }
            else
            { // type == 4, point is in the middle
                return StsMath.interpolateBilinear(trace00, trace01, trace10, trace11, nSlices, dx, dy, values);
            }
            return true;
        }

        float getValue(int nSlice)
        {
            nSlice = StsMath.minMax(nSlice, 0, nSlices - 1);
            if (type == 0)
            { // point is in lower-left corner
                return trace00[nSlice];
            }
            else if (type == 1)
            { // point is along bottom row
                return (1.0f - dx) * trace00[nSlice] + dx * trace01[nSlice];
            }
            else if (type == 3)
            { // point is along left column
                return (1.0f - dy) * trace00[nSlice] + dy * trace10[nSlice];
            }
            else
            { // type == 4, point is in the middle
                float value1 = (1.0f - dx) * trace00[nSlice] +
                    dx * trace01[nSlice];
                float value2 = (1.0f - dx) * trace10[nSlice] +
                    dx * trace11[nSlice];
                return (1.0f - dy) * value1 + dy * value2;
            }
        }
    }


    public float[] getSurfaceTimes(int row, int col)
	{
		float[] surfaceTimes = new float[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
            StsSurface surface = (StsSurface)surfaces.getElement(n);
            surfaceTimes[n] = surface.pointsZ[row][col];
		}
		return surfaceTimes;
	}


    public float[] getSurfaceDepths(int row, int col)
	{
		float[] surfaceDepths = new float[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
            StsSurface surface = (StsSurface)surfaces.getElement(n);
            surfaceDepths[n] = surface.adjPointsZ[row][col];
		}
		return surfaceDepths;
	}

	public float[] getIntervalVelocities(int row, int col)
	{
		float[] velocities = new float[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
            StsVelocityGrid velocityGrid = (StsVelocityGrid)intervalVelocityGrids.getElement(n);
            velocities[n] = (float)velocityGrid.velocities[row][col];
		}
		return velocities;
	}

	/** get constant interval velocities for each interval: datum to first surface, between surfaces, and below last surface.
	 *  length is nSurfaces+1.
	 * @return array of intervalVelocities
	 */
	public float[] getOneWayConstantIntervalVelocities()
	{
		return oneWayConstantIntervalVelocities;
	}

	public float getOneWayConstantBottomIntervalVelocity()
	{
		return oneWayBottomConstantIntervalVelocity;
	}

    public float getDataMin()
	{
		return super.getDataMin();
	}

	public void setDataMin(float dataMin)
	{
		this.dataMin = dataMin;
	}

    public void setDataRange(float dataMin, float dataMax)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
    }

    public float getDataMax()
	{
		return super.getDataMax();
	}

	public void setDataMax(float dataMax)
	{
		this.dataMax = dataMax;
	}

	public float getTimeMin()
	{
		return timeMin;
	}

	public float getTimeMax()
	{
		return timeMax;
	}

	public float getTimeInc()
	{
		return timeInc;
	}

	public StsSeismicVolume getInputVelocityVolume()
	{
		return inputVelocityVolume;
	}

	public StsSeismicVolume getVelocityVolume()
	{
		return velocityVolume;
	}

	public StsVelocityGrid[] getVelocityGrids()
	{
		if (intervalVelocityGrids == null) return null;
		return (StsVelocityGrid[]) intervalVelocityGrids.getCastList();
	}

	public StsModelSurface[] getTimeSurfaces()
	{
		if (surfaces == null) return null;
		return (StsModelSurface[]) surfaces.getCastList();
	}

	public int getNSurfaces()
	{
		return nSurfaces;
	}

    public float getMinDepthAtTime(float time)
    {
        double minDepth = getZFromRowCol(0.0f, 0.0f, time);
        minDepth = Math.min(minDepth, getZFromRowCol(0.0f, (float)(nCols-1), time));
        minDepth = Math.min(minDepth, getZFromRowCol((float)(nRows-1), (float)(nCols-1), time));
        minDepth = Math.min(minDepth, getZFromRowCol((float)(nRows-1), 0.0f, time));
        minDepth = Math.min(minDepth, getZFromRowCol((float)(nRows-1), 0.0f, time));
        minDepth = Math.min(minDepth, getZFromRowCol((float)(nRows-1)/2, (float)(nCols-1)/2, time));
        return (float)minDepth;
    }

    public float getMaxDepthAtTime(float time)
    {
        double maxDepth = getZFromRowCol(0.0f, 0.0f, time);
        maxDepth = Math.max(maxDepth, getZFromRowCol(0.0f, (float)(nCols-1), time));
        maxDepth = Math.max(maxDepth, getZFromRowCol((float)(nRows-1), (float)(nCols-1), time));
        maxDepth = Math.max(maxDepth, getZFromRowCol((float)(nRows-1), 0.0f, time));
        maxDepth = Math.max(maxDepth, getZFromRowCol((float)(nRows-1), 0.0f, time));
        maxDepth = Math.max(maxDepth, getZFromRowCol((float)(nRows-1)/2, (float)(nCols-1)/2, time));
        return (float)maxDepth;
    }
 

	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
		{
		    new StsBooleanFieldBean(StsSeismicVelocityModel.class, "isVisible", "Visible")
	    };
        return displayFields;
	}

	public StsFieldBean[] getPropertyFields()
	{
        if(propertyFields != null) return propertyFields;

	    propertyFields = new StsFieldBean[]
        {
            new StsFloatFieldBean(StsSeismicVelocityModel.class, "dataMin", 0.0f, 10.0f, "Data Min"),
            new StsFloatFieldBean(StsSeismicVelocityModel.class, "dataMax", 0.0f, 10.0f, "Data Max")
//            colorscaleBean
	    };
        return propertyFields;
	}

	public StsObjectPanel getObjectPanel()
	{
		if (velocityVolumeObjectPanel == null)
		{
			velocityVolumeObjectPanel = StsObjectPanel.constructor(this, true);
		}
		return velocityVolumeObjectPanel;
	}

	public void treeObjectSelected()
	{
		StsSeismicVelocityModelClass velocityVolumeClass = (StsSeismicVelocityModelClass) currentModel.getStsClass("com.Sts.DBTypes.StsSeismicVelocityModel");
		velocityVolumeClass.selected(this);
		currentModel.getGlPanel3d().checkAddView(StsView3d.class);
		currentModel.win3dDisplayAll();
	}

	public boolean anyDependencies()
	{
		return false;
	}

	public Object[] getChildren()
	{
		return new Object[0];
	}
/*
	public boolean allocateVolumes()
	{
		try
		{
			dataCubes = new StsByteDataCube[3];
			dataCubes[XDIR] = new StsByteDataCube(XDIR, nRows, nCols, nCroppedSlices, null);
			dataCubes[YDIR] = new StsByteDataCube(YDIR, nRows, nCols, nCroppedSlices, null);
			dataCubes[ZDIR] = new StsByteDataCube(ZDIR, nRows, nCols, nCroppedSlices, null);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsVirtualVolume.allocateVolumes() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}
*/
/*
    public void initializeSpectrum()
	{
		try
		{
			StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
			colorscale = new StsColorscale("SeismicVelocity", spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_RWB), dataMin, dataMax);
			colorscale.setEditRange(dataMin, dataMax);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVelocityModel.initializeSpectrum() failed.", e, StsException.WARNING);
		}
	}

	public void setSpectrumDialog(StsSpectrumDialog spectrumDialog)
	{
		this.spectrumDialog = spectrumDialog;
	}
*/
	public StsSeismicVelocityModelClass getSeismicVelocityVolumeClass()
	{
		return (StsSeismicVelocityModelClass) currentModel.getCreateStsClass(getClass());
	}
/*
	public boolean getIsVisibleOnCursor()
	{
		return isVisible && getSeismicVelocityVolumeClass().getIsVisibleOnCursor();
	}
*/
    /** velocity can always be drawn when corresponding velocityVolume (an StsSeismicVolume) is selected for display. So don't display it again here. */
	public Object readPlaneData(int dir, float dirCoordinate)
	{
		return null;
//		if(velocityVolume == null) return null;
//		return velocityVolume.readPlaneData(dirNo, dirCoordinate);
	}
/*
	public Object readPlaneData(int dirNo, int nPlane)
	{
		byte[] plane = null;

		if (dirNo == XDIR)
		{
			int col = nPlane;
			plane = getTraceVelocityPlane(dirNo, col, this);
		}
		else if (dirNo == YDIR)
		{
			int row = nPlane;
			plane = getTraceVelocityPlane(dirNo, row, this);
		}
		return plane;
	}
*/
/*
   public void actionPerformed(ActionEvent e)
   {
	   if (e.getSource() instanceof StsColorscale)
	   {
		   currentModel.displayIfCursor3dObjectChanged(this);
	   }
	   else
	   {
		   System.out.println("unknown action "+e+"\n   on "+this);
		   //this.fireActionPerformed(null);
	   }
	   return;
   }
*/
/*
	public void itemStateChanged(ItemEvent e)
	{
		if (e.getItem() instanceof StsColorscale)
		{
			currentModel.displayIfCursor3dObjectChanged(this);
		}
		else
		{
			ItemEvent event = null;
			fireItemStateChanged(event);
		}
		return;
	}
*/
}
