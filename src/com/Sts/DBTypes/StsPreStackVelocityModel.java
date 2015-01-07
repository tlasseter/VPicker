package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.nio.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 15, 2007
 * Time: 10:13:35 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class StsPreStackVelocityModel extends StsSeismicBoundingBox implements StsDialogFace
{
    /** volume this velocityModel belongs to; initialized by volume */
    public StsPreStackLineSet lineSet;
    /** optional input velocity volume (3d) or lineset (2d) */
    StsSeismic seismicVelocityVolume = null;
    /** velocity is defined by this array of profiles */
    StsObjectRefList velocityProfiles;//	boolean isComplete = false;
    /** current local X coordinate */
    // float currentX;
    /** current local Y coordinate */
    // float currentY;//	float xInc;
    /** multiplier to convert imported velocities to project units */
    float multiplier = 1.0f;
    /** current units for velocity */
    String velUnits = StsParameters.VEL_FT_PER_MSEC;
    /** velocity type being displayed/exported */
    byte velocityType = StsParameters.SAMPLE_TYPE_VEL_RMS;

    /** colorscale for stack velocities */
    // TODO: currently creating four colorscales; instead should only create as needed and add to refList. TJL 1/16/07
    // TODO: Should imbed colorlist as transient in colorscale, and "manage" a bit cleaner. TJL 1/16/07
    protected StsObjectRefList velocityColorscales; // one for each of the velocityTypes: interval, rms, instantaneous, average

    /** this is the current profile being referenced or used.  Only for convenience so we don't have to search the entire list all the time. */
//    transient StsVelocityProfile currentVelocityProfile = null;

    transient StsColorList[] velocityColorLists = null; // one for each of the velocityTypes: interval, rms, instantaneous, average
    transient StsColorscale currentColorscale = null;
    transient StsColorList currentColorList = null;

    transient StsVelocityProfile[][] currentVelocityProfiles;
    //    transient private StsPoint[] minVelocityProfile = null;
    //    transient private StsPoint[] maxVelocityProfile = null;
    /** minVelocities for all profiles */
    transient float[] minVelocities = null;
    /** maxVelocities for all profiles */
    transient float[] maxVelocities = null;
    /** minimum profile time/depth */
    transient float corridorZMin = largeFloat;
    /** maximum profile time/depth */
    transient float corridorZMax = -largeFloat;
    /** index of profileZMin in range zMin, zMax with intervals zInc */
    transient int corridorZMinIndex = 0;
    /** index of profileZMaxin range zMin, zMax with intervals zInc */
    transient int corridorZMaxIndex = 0;
    /** this is the velocity shift in from minVelocities and maxVelocities based on corridorPercentage defining the corridor */
    transient float[] corridorAdjustments = null;

    transient public StsRadialInterpolation interpolation = null;

    transient protected ExportPanelTabbedPane tabbedPane = null;

    static public final String nameExtension = "-velModel";

    static final boolean debug = false;

    public static final float DefaultVelMin = 5f;
    public static final float DefaultVelMax = 20f;
    public static final String[] VelocitySpectrumNames = StsSpectrumClass.velocitySpectrums;

    static StsComboBoxFieldBean velocityTypeBean;

//    abstract public void initializeSuperGather();

    abstract protected void constructInterpolatedVelocityPoints();

    abstract public float[] getVolumeVelocities(int row, int col);

    abstract public float[][] getVolumeVelocityProfile(int row, int col);

    public StsPreStackVelocityModel()
    {
    }

    public StsPreStackVelocityModel(boolean persistent)
    {
        super(persistent);
    }

    public StsPreStackVelocityModel(StsPreStackLineSet lineSet, String name)
    {
        this(lineSet);
        setName(new String(name)); // make copy in case we want to change it
    }

    public StsPreStackVelocityModel(StsPreStackLineSet lineSet)
    {
        super(false);
        this.lineSet = lineSet;
        initializeToBoundingBox(lineSet);
        setName(lineSet.getName() + nameExtension);
        // need to make a copy, since user might change this name independent of seismic volume name
        this.isDataFloat = true;
        /*
            currentX = lineSet.getXMin();
            currentY = lineSet.getYMin();
        */
        dataMin = StsObject.getCurrentModel().getProject().convertVelocityToProjectUnits(DefaultVelMin, StsParameters.VEL_FT_PER_MSEC);
        dataMax = StsObject.getCurrentModel().getProject().convertVelocityToProjectUnits(DefaultVelMax, StsParameters.VEL_FT_PER_MSEC);
        dataMax = (float) StsMath.niceNumber(dataMax, false);
        createColorscales();
        addToModel();
    }

    public void initializeVelocityProfiles()
    {
        velocityProfiles = StsObjectRefList.constructor(100, 100, "velocityProfiles", this);
    }

    public void initializeInterpolation()
    {
        interpolation = null;
    }

    public StsPoint[] getVelocityProfilePoints(int row, int col)
    {
        StsVelocityProfile velocityProfile = getComputeVelocityProfile(row, col);
        if(velocityProfile == null) return null;
        return velocityProfile.getProfilePoints();
    }

    public StsPoint[] getInterpolatedVelocityProfilePoints(int row, int col)
    {
        StsVelocityProfile interpolatedVelocityProfile = getInterpolatedVelocityProfile(row, col);
        if(interpolatedVelocityProfile == null) return null;
        return interpolatedVelocityProfile.getProfilePoints();
    }

    public StsPoint getInterpolatedTopMute(int row, int col)
    {
        StsVelocityProfile interpolatedVelocityProfile = getComputeVelocityProfile(row, col);
        if(interpolatedVelocityProfile == null) return null;
        return interpolatedVelocityProfile.topMute;
    }

    public StsPoint getInterpolatedBottomMute(int row, int col)
    {
        StsVelocityProfile interpolatedVelocityProfile = getComputeVelocityProfile(row, col);
        if(interpolatedVelocityProfile == null) return null;
        return interpolatedVelocityProfile.bottomMute;
    }

    public void setPreStackLineSet(StsPreStackLineSet lineSet)
    {
        fieldChanged("lineSet", lineSet);
    }

    public void addToModel()
    {
        super.addToModel();
//        velocityProfiles = StsObjectRefList.constructor(100, 100, "velocityProfiles", this);
//		refreshObjectPanel();
    }

    public boolean initialize(StsModel model)
    {
        initializeColorscales();
        initializeCurrentVelocityProfile();
        if(velocityProfiles != null)
        {
            checkInitializeInterpolation();
//            interpolation.run();
        }
        isDataFloat = true;
        return true;
    }

    private void initializeCurrentVelocityProfile()
    {
        if(velocityProfiles == null || velocityProfiles.getSize() == 0 || lineSet == null) return;
//        int row = lineSet.getNearestBoundedRowCoor(currentY);
//        int col = lineSet.getNearestBoundedColCoor(currentX);
//        setCurrentVelocityProfile(row, col);
    }

    protected void createColorscales()
    {
        try
        {
            if(velocityColorscales != null) return;
            String[] velocityNames = StsParameters.VEL_STRINGS;
            int nColorscales = velocityNames.length;
            velocityColorscales = StsObjectRefList.constructor(nColorscales, 2, "velocityColorscales", this);
            for(int n = 0; n < nColorscales; n++)
            {
                StsSpectrum spectrum = currentModel.getSpectrum(VelocitySpectrumNames[n]);
                StsColorscale colorscale = new StsColorscale(velocityNames[n], spectrum, dataMin, dataMax);
                velocityColorscales.add(colorscale);
            }
            initializeColorscales();
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackLineSet3d.initializeColorscale() failed.", e, StsException.WARNING);
        }
    }

    protected void initializeColorscales()
    {
        // JKF -- added this because createColorscales() was only being called in a constructor not used when loading db.
        if(velocityColorscales == null)
            createColorscales();

        int nColorscales = velocityColorscales.getSize();
        velocityColorLists = new StsColorList[nColorscales];
        for(int n = 0; n < nColorscales; n++)
        {
            StsColorscale colorscale = (StsColorscale) velocityColorscales.getElement(n);
            velocityColorLists[n] = new StsColorList(colorscale);
            colorscale.addActionListener(this);
        }
        currentColorscale = (StsColorscale) velocityColorscales.getElement(velocityType);
        currentColorList = velocityColorLists[velocityType];
    }

    private StsColorscale addColorscale(String name, float velocityMin, float velocityMax)
    {
        StsSpectrum spectrum = currentModel.getSpectrum(name);
        StsColorscale colorscale = new StsColorscale(name, spectrum, velocityMin, velocityMax);
        setCurrentColorscale(colorscale);
        colorscale.addActionListener(this);
        return colorscale;
    }

    public void checkCurrentVelocityProfile(StsWin3dBase window)
    {
        StsVelocityProfile velocityProfile = getCurrentVelocityProfile(window);
        if(velocityProfile != null) velocityProfile.checkCurrentVelocityProfile();
    }

    /**
     * the line or xline has changed, so set the current velocityProfile for this gather.
     * If we have an existing profile, commit it if changed.  If the velocityVolume has
     * already created a profile (in the objectRefList), set it.  If not,  set null.
     */
/*
    public void setCurrentVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile = getComputeVelocityProfile(row, col);
        if(velocityProfile == null)
        {
            if(debug) System.out.println("No velocity profile for row " + row + " col " + col);
            return;
        }
        if(currentVelocityProfile != null && !currentVelocityProfile.sameRowCol(velocityProfile))
		{
            currentVelocityProfile.checkCurrentVelocityProfile();
        }
        currentVelocityProfile = velocityProfile;
    }
*/
/*
    public void initializeCurrentVelocityProfile(int row, int col)
    {
        if(currentVelocityProfile != null)
        {
            if(currentVelocityProfile.row == row && currentVelocityProfile.col == col)
            {
                if(debug)
                    System.out.println("Continuing use of velocity profile " + currentVelocityProfile.toString());
                return;
            }
            if(currentVelocityProfile.changeType != StsVelocityProfile.CHANGE_NONE)
            {
                if(debug)
                    System.out.println("Commiting changes to velocity profile " + currentVelocityProfile.toString());
//                commitCurrentVelocityProfile();
            }
        }

		  // if a plane of data is already available in the current display mode then don't reset the
		  // mode back to velocity. If no data is available then set it back.
		  if (!(lineSet.isPlaneOK(lineSet.lineSetClass.getDisplayMode(), StsRotatedGridBoundingBox.XDIR, col) &&
				  lineSet.isPlaneOK(lineSet.lineSetClass.getDisplayMode(), StsRotatedGridBoundingBox.YDIR, row)))
		  {
			  lineSet.lineSetClass.setDisplayMode(StsPreStackLineSetClass.DISPLAY_MODE_VELOCITY);
		  }
        if(velocityProfiles == null) return;

        for(int n = 0; n < velocityProfiles.getSize(); n++)
        {
            currentVelocityProfile = (StsVelocityProfile)velocityProfiles.getElement(n);
            if(currentVelocityProfile.row == row && currentVelocityProfile.col == col)
            {
                if(debug)
                    System.out.println("Moved to existing velocity profile " + currentVelocityProfile.toString());
                return;
            }
        }
        if(debug) System.out.println("No velocity profile for row " + row + " col " + col);
        currentVelocityProfile = null;
    }
*/
/*
    public void initializeCurrentVelocityProfile(StsVelocityProfile profile)
    {
        currentVelocityProfile = profile;
    }
*/
    public boolean hasVelocityProfile(int row, int col)
    {
        return (interpolation != null) && interpolation.isDataPoint(row, col);
    }

    public StsVelocityProfile getExistingVelocityProfile(int row, int col)
    {
        if(!hasVelocityProfile(row, col)) return null;
        return currentVelocityProfiles[row][col];
    }

    public int getNProfiles()
    {
        if(velocityProfiles == null) return 0;
        return velocityProfiles.getSize();
    }

    public float[] getVelocities(int row, int col)
    {
        int nVolumeSlices = lineSet.getNSlices(row);
        float volumeZMin = lineSet.getZMin(row);
        float volumeZInc = lineSet.getZInc(row);

        StsVelocityProfile velocityProfile = getComputeVelocityProfile(row, col);
//        StsVelocityProfile velocityProfile = getExistingVelocityProfile(row, col);
        if(velocityProfile == null || velocityProfile.getProfilePoints().length < 2) return null;

        float[] velocities = new float[nVolumeSlices];
        float z = volumeZMin;
        for(int n = 0; n < nVolumeSlices; n++, z += volumeZInc)
            velocities[n] = velocityProfile.getVelocityFromTzero(z);
        return velocities;
    }

    public float getMultiplier()
    {
        return multiplier;
    }

    public void setMultiplier(float mult)
    {
        multiplier = mult;
    }

    public StsObjectRefList getVelocityProfiles()
    {
        return velocityProfiles;
    }

    public void setVelocityProfiles(StsObjectRefList velocityProfiles)
    {
        //       if(this.velocityProfiles != null) velocityProfiles.delete();
        //       this.velocityProfiles = velocityProfiles;  // this is done by the fieldChanged method below/ TJL 9/10/07.
        fieldChanged("velocityProfiles", velocityProfiles);
    }

    public void addProfile(StsVelocityProfile profile)
    {
//        currentVelocityProfile = profile
        if(!profile.isPersistent())
        {
            currentModel.add(profile);
            currentModel.commit();
        }
        profile.checkInitializeMutes(lineSet);
//        profile.interpolated = false;
        int row = profile.row;
        int col = profile.col;
//        interpolatedVelocityProfiles[row][col] = profile;
        insertProfileInList(profile);
//        velocityProfiles.sort();
        // add this point to collection of new points (should be only this one in the list)
        interpolation.addDataPoint(profile, row, col);
        // redo interpolation with this new point added
//        interpolation.run();
    }

    private void insertProfileInList(StsVelocityProfile velocityProfile)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles == 0)
        {
            velocityProfiles.add(velocityProfile);
            return;
        }
        int rowColIndex = lineSet.getVolumeRowColIndex(velocityProfile.row, velocityProfile.col);
        StsVelocityProfile nextProfile = (StsVelocityProfile) velocityProfiles.getElement(0);
        int nextRowColIndex = -1;
        for(int n = 0; n < nProfiles; n++)
        {
            int prevRowColIndex = nextRowColIndex;
            nextProfile = (StsVelocityProfile) velocityProfiles.getElement(n);
            nextRowColIndex = lineSet.getVolumeRowColIndex(nextProfile.row, nextProfile.col);
            if(rowColIndex > prevRowColIndex && rowColIndex < nextRowColIndex)
            {
                velocityProfiles.insertBefore(n, velocityProfile);
                return;
            }
        }
        velocityProfiles.add(velocityProfile);
    }

    /*
        public boolean interpolateVelocityVolume(StsVelocityProfile velocityProfile)
        {
            if(!checkInitializeInterpolation()) return false;
            interpolation.addDataPoint(velocityProfile, velocityProfile.row, velocityProfile.col);
            interpolation.run();
    //        interpolation.runUpdate(velocityProfile, velocityProfile.row, velocityProfile.col, false);
            return true;
        }
    */
    public void updateVelocityProfile(StsVelocityProfile velocityProfile, int indexPicked)
    {
        if(velocityProfile == null) return;
        if(velocityProfile.changeType == StsVelocityProfile.CHANGE_NONE) return;
        updateInterpolation(velocityProfile);
        StsPoint[] oldPoints = velocityProfile.getInitialProfilePoints();
        adjustProfileCorridor(velocityProfile, indexPicked, oldPoints);
        corridorAdjustments = null;
    }

    public void updateInterpolation(StsVelocityProfile velocityProfile)
    {
        if(velocityProfile == null) return;
        if(!checkRunInterpolation()) return;
        int row = velocityProfile.row;
        int col = velocityProfile.col;
        interpolation.updateDataPoint(velocityProfile, row, col);
        lineSet.updatePreStackVolumes(interpolation);
//        interpolation.updateDataPoint(currentVelocityProfile, currentVelocityProfile.row, currentVelocityProfile.col);
    }

    private void updatePreStackVolumes()
    {
        lineSet.updatePreStackVolumes(interpolation);
    }

    /** Don't call this directly; call velocityProfile.delete() which calls this method. */
    public void deleteVelocityProfile(StsVelocityProfile velocityProfile)
    {
        velocityProfiles.delete(velocityProfile);
        int row = velocityProfile.row;
        int col = velocityProfile.col;
        currentVelocityProfiles[row][col] = null;
        interpolation.initialize();
        initializeInterpolateVelocityProfiles();
//        interpolation.run();
        // replace this actual profile with an interpolated one
        velocityProfile = getComputeVelocityProfile(row, col);
//        currentVelocityProfile = velocityProfile;
        recomputeCorridor();
    }


    private int[] getNextProfileRowCol(int row, int col)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles == 0) return null;
        int currentRowColIndex = lineSet.getVolumeRowColIndex(row, col);
        int nextRowColIndex = -1;
        for(int n = 0; n < velocityProfiles.getSize(); n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            int prevRowColIndex = nextRowColIndex;
            nextRowColIndex = lineSet.getVolumeRowColIndex(profile.row, profile.col);
            if(currentRowColIndex >= prevRowColIndex && currentRowColIndex < nextRowColIndex)
                return new int[]
                        {profile.row, profile.col};
        }
        // must be > last, so return the first
        StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(0);
        return new int[]
                {profile.row, profile.col};
    }

    private int[] getPrevProfileRowCol(int row, int col)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles == 0) return null;
        int currentRowColIndex = lineSet.getVolumeRowColIndex(row, col);
        int prevRowColIndex = Integer.MAX_VALUE;
        for(int n = nProfiles - 1; n >= 0; n--)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            int nextRowColIndex = prevRowColIndex;
            nextRowColIndex = lineSet.getVolumeRowColIndex(profile.row, profile.col);
            if(currentRowColIndex > nextRowColIndex && currentRowColIndex <= prevRowColIndex)
                return new int[]
                        {profile.row, profile.col};
        }
        // must be < first, so return the last
        StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(nProfiles - 1);
        return new int[]
                {profile.row, profile.col};
    }

    /**
     * get the subset of velocitProfiles whose row-col is within dRow&dCol of the given row-col.
     * Don't include the current profile.  Then reorder the list so it continues up from the current profile
     * and wraps around from the beginning.  The picker then picks the first one encountered under the mouse,
     * guaranteeing that it will cycle thru the picks so the user can see all picks at that same mouse location.
     */
    public StsVelocityProfile[] getNeighborProfiles(int row, int col, int dRow, int dCol)
    {
        int nVelocityProfiles = velocityProfiles.getSize();
        if(nVelocityProfiles == 0) return new StsVelocityProfile[0];
        StsVelocityProfile[] neighborProfiles = new StsVelocityProfile[(2 * dRow + 1) * (2 * dCol + 1)];
        int nNeighbors = 0;
        int minRow = Math.max(0, row - dRow);
        int maxRow = Math.min(nRows - 1, row + dRow);
        int minCol = Math.max(0, col - dCol);
        int maxCol = Math.min(nCols - 1, col + dCol);

        int nCurrentProfile = 0; // position of current profile in list
        for(int n = 0; n < nVelocityProfiles; n++)
        {
            StsVelocityProfile velocityProfile = (StsVelocityProfile) velocityProfiles.getElement(n);
            int profileRow = velocityProfile.row;
            int profileCol = velocityProfile.col;
            if(profileRow < minRow) continue;
            if(profileRow > maxRow) break;
            if(profileRow == row && profileCol == col)
            {
                nCurrentProfile = nNeighbors;
                continue;
            }
            if(profileCol >= minCol && profileCol <= maxCol)
                neighborProfiles[nNeighbors++] = velocityProfile;
        }
        StsVelocityProfile[] wrappedNeighborProfiles = new StsVelocityProfile[nNeighbors];
        int m = 0;
        for(int n = nCurrentProfile; n < nNeighbors; n++, m++)
            wrappedNeighborProfiles[m] = neighborProfiles[n];
        for(int n = 0; n < nCurrentProfile; n++, m++)
            wrappedNeighborProfiles[m] = neighborProfiles[n];
        return wrappedNeighborProfiles;
    }

    public boolean canExport()
    {
        return true;
    }

    public boolean export()
    {
        try
        {
            /*
                exportStartInline = line2d.getRowNumMin();
                exportEndInline = line2d.getRowNumMax();
                exportStartCrossline = line2d.getColNumMin();
                exportEndCrossline = line2d.getColNumMax();
                exportIncCrossline = line2d.getColNumInc();
                exportIncInline = line2d.getRowNumInc();
                exportStartZ = line2d.getZMin();
                exportEndZ = line2d.getZMax();
                exportIncZ = line2d.getZInc();
            */
            StsProcessDismissDialog exportDialog = new StsProcessDismissDialog(currentModel.win3d, this, "Velocity Model Export", true);
            return true;
        }
        catch(Exception e)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Prestack velocity model export failed.");
            return false;
        }
    }

    public float getMinDepthAtTime(float time)
    {
        int nVelocityProfiles = this.getNProfiles();
        float minVelocity = StsParameters.largeFloat;
        for(int n = 0; n < nVelocityProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            float velocity = profile.getVelocityFromTzero(time);
            minVelocity = Math.min(minVelocity, velocity);
        }
	    return time*minVelocity;
    }

    public float getMaxDepthAtTime(float time)
    {
        float maxVelocity = 0;
        int nVelocityProfiles = this.getNProfiles();
        for(int n = 0; n < nVelocityProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            float velocity = profile.getVelocityFromTzero(time);
            maxVelocity = Math.max(maxVelocity, velocity);
        }
        return time*maxVelocity;
    }

    public void dialogSelectionType(int type)
    {
        StsDialogFace dialogFace = tabbedPane.getDialogFace();
        dialogFace.dialogSelectionType(type);
    }

    public StsPreStackLineSet getSeismicVolume()
    {
        return lineSet;
    }

    public boolean hasProfiles()
    {
        if(velocityProfiles == null) return false;
        return this.velocityProfiles.getSize() > 0;
    }

    public String getVelocityUnitsString()
    {
        return velUnits;
    }

    public void setVelocityUnitsString(String option)
    {
        velUnits = option;
        setMultiplier(currentModel.getProject().calculateVelScaleMultiplier(velUnits));
        dbFieldChanged("velUnits", velUnits);
    }

    public boolean recomputeCorridor()
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles == 0) return false;

        computeProfileZRange();
        //       profileZMinIndex = this.getNearestSliceCoor(profileZMin);
        //       profileZMaxIndex = getNearestSliceCoor(profileZMax);
        int nProfileSamples = corridorZMaxIndex - corridorZMinIndex + 1;
        int nGoodProfiles = 0;
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            if(profile.getNProfilePoints() >= 2)
            {
                nGoodProfiles++;
                if(nGoodProfiles >= 2) break;
            }
        }
        //        if(nGoodProfiles < 2) return false;

        minVelocities = new float[nProfileSamples];
        maxVelocities = new float[nProfileSamples];
        corridorAdjustments = null;
        int p = 0;
        // initializeSuperGather min and max velocity profiles to first profile that has minimum of 2 points
        for(p = 0; p < nProfiles; p++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(p);
            StsPoint[] points = profile.getProfilePoints();
            int nPoints = points.length;
            if(nPoints < 2) continue;
            float[] velocities = getProfileVelocityFloats(profile, corridorZMin, zInc, nProfileSamples);
            if(velocities == null) continue;
            for(int n = 0; n < nProfileSamples; n++)
            {
                minVelocities[n] = velocities[n];
                maxVelocities[n] = velocities[n];
            }
            break;
        }
        for(; p < nProfiles; p++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(p);

            StsPoint[] points = profile.getProfilePoints();
            int nPoints = points.length;
            if(nPoints < 2) continue;
            float[] velocities = getProfileVelocityFloats(profile, corridorZMin, zInc, nProfileSamples);
            if(velocities == null) continue;
            for(int n = 0; n < nProfileSamples; n++)
            {
                if(velocities[n] < minVelocities[n]) minVelocities[n] = velocities[n];
                else if(velocities[n] > maxVelocities[n]) maxVelocities[n] = velocities[n];
            }
        }
        return true;
    }

    private void computeProfileZRange()
    {
        int nProfiles = velocityProfiles.getSize();
        corridorZMin = largeFloat;
        corridorZMax = -largeFloat;
        boolean initialized = false;
        for(int p = 0; p < nProfiles; p++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(p);
            StsPoint[] points = profile.getProfilePoints();
            int nPoints = points.length;
            if(nPoints < 2) continue;
            corridorZMin = Math.min(corridorZMin, points[0].v[1]);
            corridorZMax = Math.max(corridorZMax, points[nPoints - 1].v[1]);
            initialized = true;
        }
        if(initialized)
        {
            corridorZMin = StsMath.intervalRoundDown(corridorZMin, zMin, zInc);
            corridorZMax = StsMath.intervalRoundUp(corridorZMax, zMin, zInc);
        }
        else
        {
            corridorZMin = zMin;
            corridorZMax = zMax;
        }
        corridorZMinIndex = getNearestBoundedSliceCoor(corridorZMin);
        corridorZMaxIndex = getNearestBoundedSliceCoor(corridorZMax);
        corridorZMin = getZCoor(corridorZMinIndex);
        corridorZMax = getZCoor(corridorZMaxIndex);
    }

    private float[] getProfileVelocityFloats(int nProfile)
    {
        StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(nProfile);
        return getProfileVelocityFloats(profile, zMin, zInc, nSlices);
    }

    private float[] getProfileVelocityFloats(StsVelocityProfile profile, float zMin, float zInc, int nValues)
    {
        StsPoint[] profilePoints = profile.getProfilePoints();
        if(profilePoints == null) return null;
        return StsMath.interpolateValues(profilePoints, nValues, zMin, zInc, 1, 0);
    }

    public void addProfileToCorridor(StsVelocityProfile newProfile, float zMin, float zInc, int nValues)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles <= 1) return;
        if(minVelocities == null)
        {
            recomputeCorridor();
            return;
        }
        float[] newVelocities = getProfileVelocityFloats(newProfile, corridorZMin, zInc, nValues);
        for(int n = 0; n < nValues; n++)
        {
            int where = StsMath.belowBetweenAbove(newVelocities[n], minVelocities[n], maxVelocities[n]);
            if(where == StsMath.BELOW)
                minVelocities[n] = newVelocities[n];
            else if(where == StsMath.ABOVE)
                maxVelocities[n] = newVelocities[n];
        }
    }

    static int MOVE_INSIDE = 0;
    static int MOVE_ABOVE = 1;
    static int MOVE_BELOW = -1;

    public void adjustProfileCorridor(StsVelocityProfile adjustedProfile, int indexChanged, StsPoint[] originalPoints)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles <= 1) return;

        if(minVelocities == null)
        {
            recomputeCorridor();
            return;
        }

        StsPoint[] points = adjustedProfile.getProfilePoints();
        int nPoints = points.length;
        int nOriginalPoints = originalPoints.length;
        if(nPoints < 2) return;

        if(checkRecomputeCorridor(points, indexChanged)) return;
        int pointIndexMin = Math.max(0, indexChanged - 1);
        int pointIndexMax = Math.min(nOriginalPoints - 1, indexChanged + 1);
        int indexMin, indexMax;
        if(pointIndexMin == 0)
            indexMin = 0;
        else
        {
            float intervalZMin = points[pointIndexMin].v[1];
            float pointZMin = StsMath.intervalRoundDown(intervalZMin, zMin, zInc);
            indexMin = Math.round((pointZMin - corridorZMin) / zInc);
        }
        if(pointIndexMax == nOriginalPoints - 1)
            indexMax = corridorZMaxIndex - corridorZMinIndex;
        else
        {
            float intervalZMax = points[pointIndexMax].v[1];
            float pointZMax = StsMath.intervalRoundUp(intervalZMax, zMin, zInc);
            indexMax = Math.round((pointZMax - corridorZMin) / zInc);
        }

        if(mustAdjustProfileCorridor(points, originalPoints, indexMin, indexMax))
            adjustProfileCorridor(indexMin, indexMax);
    }

    /**
     * If adjusted point z is outside of current corridor range, recompute the entire corridor.
     * TODO It should be possible to just extrapolate the existing corridor to cover this Z, and then adjust the corridor in the affected range.
     *
     * @return true if it does requiring that complete corridor be recomputed (see TODO).
     */
    private boolean checkRecomputeCorridor(StsPoint[] points, int indexChanged)
    {
        int nPoints = points.length;
        if(nPoints == 2)
        {
            recomputeCorridor();
            return true;
        }
        if(indexChanged >= points.length)
        {  //user deleted last vertex
            recomputeCorridor();
            return true;
        }
        float adjustedZ = points[indexChanged].v[1];
        if(adjustedZ < corridorZMin || adjustedZ > corridorZMax)
        {
            recomputeCorridor();
            return true;
        }
        return false;
    }

    private boolean mustAdjustProfileCorridor(StsPoint[] newPoints, StsPoint[] oldPoints, int indexMin, int indexMax)
    {
        if(oldPoints == null) return true;
        if(minVelocities == null) return true;
        float[] newVelocities = getVelocitiesOverProfileInterval(newPoints, indexMin, indexMax);
        float[] oldVelocities = getVelocitiesOverProfileInterval(oldPoints, indexMin, indexMax);
        for(int i = 0, n = indexMin; n <= indexMax; n++, i++)
        {
            int whereNew = StsMath.belowBetweenAbove(newVelocities[i], minVelocities[n], maxVelocities[n]);
            int whereOld = StsMath.belowBetweenAbove(oldVelocities[i], minVelocities[n], maxVelocities[n]);
            if(whereNew != whereOld)
                return true;
        }
        return false;
    }

    // TODO use more efficient method since points are monotonically increasing in z
    private float[] getVelocitiesOverProfileInterval(StsPoint[] points, int indexMin, int indexMax)
    {
        int nPoints = points.length;
        if(nPoints < 2) return null;
        int nVelocities = indexMax - indexMin + 1;
        float[] velocities = new float[nVelocities];

        float z = corridorZMin + indexMin * zInc;
        int nInterval = getInterval(points, z);
        float z0 = points[nInterval].v[1];
        float v0 = points[nInterval].v[0];
        float z1 = points[nInterval + 1].v[1];
        float v1 = points[nInterval + 1].v[0];
        double dvdz = (v1 - v0) / (z1 - z0);
        double a = v0 - z0 * dvdz;
        double b = dvdz;
        for(int i = 0, n = indexMin; n <= indexMax; n++, i++, z += zInc)
        {
            if(z > z1 && nInterval < nPoints - 2)
            {
                nInterval = getInterval(points, z);
                z0 = points[nInterval].v[1];
                v0 = points[nInterval].v[0];
                z1 = points[nInterval + 1].v[1];
                v1 = points[nInterval + 1].v[0];
                dvdz = (v1 - v0) / (z1 - z0);
                a = v0 - z0 * dvdz;
                b = dvdz;
            }
            velocities[i] = (float) (a + b * z);
        }
        return velocities;
    }

    private int getInterval(StsPoint[] points, float z)
    {
        int nIntervals = points.length - 1;
        for(int n = 0; n < nIntervals; n++)
            if(z < points[n + 1].v[1])
                return n;
        return nIntervals - 1;
    }

    private void adjustProfileCorridor(int indexMin, int indexMax)
    {
        for(int n = indexMin; n <= indexMax; n++)
        {
            minVelocities[n] = largeFloat;
            maxVelocities[n] = -largeFloat;
        }

        int nProfiles = velocityProfiles.getSize();
        for(int p = 0; p < nProfiles; p++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(p);
            StsPoint[] points = profile.getProfilePoints();
            float[] velocities = this.getVelocitiesOverProfileInterval(points, indexMin, indexMax);
            if(velocities == null) continue;
            for(int i = 0, n = indexMin; n <= indexMax; n++, i++)
            {
                float v = velocities[i];
                if(v < minVelocities[n]) minVelocities[n] = v;
                if(v > maxVelocities[n]) maxVelocities[n] = v;
            }
        }
    }

    public void clearVelocityProfiles()
    {
        minVelocities = null;
        maxVelocities = null;
        corridorAdjustments = null;
    }

    public void drawProfileCorridor(GL gl, float percentage)
    {
        if(minVelocities == null)
        {
            if(!recomputeCorridor())
                return;
        }

        // Apply percentage
        checkComputeCorridorAdjustments(percentage);

        if(minVelocities == null) return;
        StsGLDraw.drawDottedLineStrip2d(gl, StsColor.CYAN, minVelocities, corridorAdjustments, -1f, corridorZMin, zInc, 2);

        if(maxVelocities == null) return;
        StsGLDraw.drawDottedLineStrip2d(gl, StsColor.CYAN, maxVelocities, corridorAdjustments, 1f, corridorZMin, zInc, 2);
    }

    private boolean checkComputeCorridorAdjustments(float percentage)
    {
        if(corridorAdjustments != null) return true;
        if(minVelocities == null) return false;
        int nProfileSamples = corridorZMaxIndex - corridorZMinIndex + 1;
        corridorAdjustments = new float[nProfileSamples];
        for(int n = 0; n < nProfileSamples; n++)
        {
            float velDiff = maxVelocities[n] - minVelocities[n];
            corridorAdjustments[n] = velDiff * ((100 - percentage) / 100.0f) * 0.5f;
        }
        return true;
    }

    public void drawInterpolatedProfile(StsGLPanel3d glPanel3d, StsPoint[] points, StsColor stsColor, float percentage)
    {
        GL gl = glPanel3d.getGL();

        float lineWidth = 2.0f;

        int nPoints = points.length;

        try
        {
            boolean corridorAdjust = (nPoints >= 2) && checkComputeCorridorAdjustments(percentage);
            if(corridorAdjust)
            {
                int nVelocityPoints = corridorZMaxIndex - corridorZMinIndex + 1;
                float[] velocities = StsMath.interpolateValues(points, nVelocityPoints, corridorZMin, zInc, 1, 0);
                if(velocities == null) return;
                drawCorridorLimitedProfile(gl, velocities, nVelocityPoints, stsColor, lineWidth, true);
            }
            else
            {
                StsGLDraw.drawDottedLineStrip2d(gl, stsColor, points, 2);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawInterpolatedProfile", e);
        }
    }

    private void drawCorridorLimitedProfile(GL gl, float[] velocities, int nVelocityPoints, StsColor stsColor, float lineWidth, boolean stippled)
    {
        try
        {
            float z = corridorZMin;
            if(minVelocities == null) return;

            gl.glDisable(GL.GL_LIGHTING);
            if(stippled)
            {
                gl.glLineStipple(1, StsGraphicParameters.dottedLine);
                gl.glEnable(GL.GL_LINE_STIPPLE);
            }
            stsColor.setGLColor(gl);
            gl.glLineWidth(lineWidth);

            gl.glBegin(GL.GL_LINE_STRIP);
            for(int n = 0; n < nVelocityPoints; n++, z += zInc)
            {
                float velocity = StsMath.minMax(velocities[n], minVelocities[n] + corridorAdjustments[n], maxVelocities[n] - corridorAdjustments[n]);
                gl.glVertex2f(velocity, z);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawCorridorLimitedProfile", e);
        }
        finally
        {
            gl.glEnd();
            gl.glLineWidth(1.0f);
            if(stippled) gl.glDisable(GL.GL_LINE_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public void drawProfile(StsGLPanel3d glPanel3d, StsPoint[] points, StsColor stsColor, float percentage, boolean drawVertices, boolean active, boolean drawLabels, int pickedIdx, boolean drawBackground)
    {
        GL gl = glPanel3d.getGL();

        int width;
        if(active)
            width = 4;
        else
            width = 2;

        int nPoints = points.length;

        try
        {
            int nProfiles = velocityProfiles.getSize();
            boolean corridorAdjust = (nProfiles > 1) && (nPoints >= 2) && checkComputeCorridorAdjustments(percentage);
            if(corridorAdjust)
            {
                int nVelocityPoints = corridorZMaxIndex - corridorZMinIndex + 1;
                float[] velocities = StsMath.interpolateValues(points, nVelocityPoints, corridorZMin, zInc, 1, 0);
                if(velocities == null) return;
                drawCorridorLimitedProfile(gl, velocities, nVelocityPoints, StsColor.BLACK, width, false);
                drawCorridorLimitedProfile(gl, velocities, nVelocityPoints, stsColor, width / 2, false);
            }
            else
            {
                StsGLDraw.drawLineStrip2d(gl, StsColor.BLACK, points, 2, width);
                StsGLDraw.drawLineStrip2d(gl, stsColor, points, 2, width / 2);
            }

            if(!drawVertices) return;
            //            StsColor.WHITE.setGLColor(gl);
            StsVelocityProfile.drawVertices(points, gl, pickedIdx);
            if(drawLabels) StsVelocityProfile.drawLabels(points, gl, drawBackground);
        }
        catch(Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
    }

    public String getName()
    {
        if(lineSet == null) return "nullVolume" + nameExtension;
        else return super.getName();
    }

    public void setVelocityTypeString(String typeString)
    {
        byte velocityType = StsParameters.getVelocityTypeFromString(typeString);
        if(this.velocityType == velocityType) return;
        this.velocityType = velocityType;
        currentModel.viewObjectChangedAndRepaint(this, lineSet);
    }

    public String getVelocityTypeString()
    {
        return StsParameters.VEL_STRINGS[velocityType];
    }

    public void setVelocityType(byte type)
    {
        this.velocityType = type;
        //StsPreStackVelocityModel3d.velocityTypeBean.setValueObject(StsParameters.VEL_STRINGS[type]);
        if(velocityTypeBean != null) velocityTypeBean.setValueObject(StsParameters.VEL_STRINGS[type]);
        currentColorscale = (StsColorscale) velocityColorscales.getElement(velocityType);
        currentColorList = velocityColorLists[velocityType];
    }

    public byte getVelocityType()
    {
        return velocityType;
    }

    public void actionPerformed(ActionEvent e)
    {
        // actionEvent fired when colorscale edit closed
        if(!(e.getSource() instanceof StsColorscale)) return;
        StsColorscale colorscale = (StsColorscale) e.getSource();
        int nColorscales = velocityColorscales.getSize();
        for(int n = 0; n < nColorscales; n++)
        {
            StsColorscale velocityColorscale = (StsColorscale) velocityColorscales.getElement(n);
            if(velocityColorscale == colorscale)
            {
                velocityColorLists[n].setColorListChanged(true);
                currentModel.viewObjectRepaint(this, lineSet);
                currentModel.displayIfCursor3dObjectChanged(lineSet);
                return;
            }
        }
    }

    public void setCurrentColorscale(StsColorscale colorscale)
    {
        if(currentColorscale == colorscale) return;
        currentColorscale = colorscale;
        if(StsPreStackVelocityModel3d.colorscaleBean != null)
            StsPreStackVelocityModel3d.colorscaleBean.setValueObject(colorscale);
        setCurrentColorList(colorscale);
    }

    public StsColorscale getCurrentColorscale()
    {
        return currentColorscale;
    }

    public void setCurrentColorList(StsColorscale colorscale)
    {
        if(velocityColorLists != null)
        {
            for(int n = 0; n < velocityColorLists.length; n++)
            {
                if(velocityColorLists[n].colorscale == colorscale)
                {
                    currentColorList = velocityColorLists[n];
                    return;
                }
            }
        }
        addColorList(colorscale);
    }

    public StsColorList addColorList(StsColorscale colorscale)
    {
        currentColorList = new StsColorList(colorscale);
        velocityColorLists = (StsColorList[]) StsMath.arrayAddElement(velocityColorLists, currentColorList);
        return currentColorList;
    }

    public StsColorList getVelocityColorList()
    {
        return currentColorList;
    }

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }

    /** @return true if initialized or able to initialize; return false if cannot initialize */
    public boolean checkInitializeInterpolation()
    {
        try
        {
            if(interpolation != null) return true;
            int nMinPoints = 1;
            int nMaxPoints = 10;
            StsPreStackLineSetClass lineSetClass = (StsPreStackLineSetClass) lineSet.getStsClass();
            float neighborRadius = lineSetClass.getStackNeighborRadius();
            byte stackOption = lineSetClass.getStackOption();
            interpolation = StsRadialInterpolation.constructor(lineSet, nMinPoints, nMaxPoints, neighborRadius);
            StsMethod updateMethod = new StsMethod(this, "clearInterpolatedVelocityProfile", new Class[]
                    {Integer.class, Integer.class});
            interpolation.setUpdateMethod(updateMethod);
            constructInterpolatedVelocityPoints();
            initializeInterpolateVelocityProfiles();
            //			if(!volumeInterpolation.complete) volumeInterpolation.run();
            if(debug) System.out.println("Constructed velocity interpolator.");
            //            interpolation.run();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackVelocityModel3d:initializeInterpolation() failed.", e, StsException.WARNING);
            return false;
        }
    }

    /**
     * initialize interpolator if not already initialized, then run interpolation if not up to date
     *
     * @return true if interpolation is up to date or interpolator has been successfully run
     */
    public boolean checkRunInterpolation()
    {
        if(!checkInitializeInterpolation()) return false;
        if(interpolation.upToDate) return true;
        interpolation.run();
        return true;
    }

    public void initializeInterpolateVelocityProfiles()
    {
        int nProfiles = velocityProfiles.getSize();
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            interpolation.addDataPoint(profile, profile.row, profile.col);
            currentVelocityProfiles[profile.row][profile.col] = profile;
        }
        interpolation.run();
    }

    public void clearInterpolatedVelocityProfile(Integer rowInteger, Integer colInteger)
    {
        int row = rowInteger.intValue();
        int col = colInteger.intValue();
        this.currentVelocityProfiles[row][col] = null;
        // JKF This is causing an issue with reloading planOK statuses. Commenting out for now. Discuss with Tom.
        //        if(lineSet instanceof StsPreStackLineSet3d)
        //            ((StsPreStackLineSet3d)lineSet).setPlaneFalse(row, col);
    }

    /*
        public StsVelocityProfile getInterpolatedVelocityProfile()
        {
            StsSuperGather gather = lineSet.getGather();
            return getExistingOrInterpolatedVelocityProfile(gather.superGatherRow, gather.superGatherCol);
        }

        public StsVelocityProfile getExistingOrInterpolatedVelocityProfile(int row, int col)
        {
            if(currentVelocityProfiles == null) return null;
            return currentVelocityProfiles[row][col];
        }
    */
    public StsVelocityProfile getInterpolatedVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile = null;
        if(!checkInitializeInterpolation()) return null;
        velocityProfile = currentVelocityProfiles[row][col];
        if(!velocityProfile.isInterpolated()) return null;
        return getComputeVelocityProfile(row, col);
    }

    public StsVelocityProfile getCurrentVelocityProfile(StsWin3dBase window)
    {
        StsSuperGather gather = lineSet.getSuperGather(window);
        if(gather == null) return null;
        return getVelocityProfile(gather.superGatherRow, gather.superGatherCol);
//        return getComputeVelocityProfile(gather.superGatherRow, gather.superGatherCol);
//        currentVelocityProfile = getComputeVelocityProfile(gather.superGatherRow, gather.superGatherCol);
//        return currentVelocityProfile;
    }
/*
    public StsVelocityProfile getComputeVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile = getExistingVelocityProfile(row, col);
        if(velocityProfile != null)
        {
            if(debug) System.out.println("Found existing velocityProfile at row " + row + " col " + col);
            return velocityProfile;
        }
        velocityProfile = getInterpolatedVelocityProfile(row, col);
        if(velocityProfile != null)
        {
            if(debug) System.out.println("Interpolated velocityProfile at row " + row + " col " + col);
            return velocityProfile;
        }
        else
            return null;
    }
*/

    public StsVelocityProfile getVelocityProfile(int row, int col)
    {
        if(currentVelocityProfiles == null) return null;
        if(!this.isInsideRowColRange(row, col)) return null;
        return currentVelocityProfiles[row][col];
    }

    public StsVelocityProfile getComputeVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile = null;
        if(!checkInitializeInterpolation()) return null;
        if(currentVelocityProfiles == null) return null;
        if(row >= currentVelocityProfiles.length || col >= currentVelocityProfiles[row].length) return null;
        velocityProfile = currentVelocityProfiles[row][col];
        if(velocityProfile != null)
        {
            if(debug) System.out.println("Existing velocityProfile at row " + row + " col " + col);
            return velocityProfile;
        }
        velocityProfile = computeInterpolatedVelocityProfile(row, col);
        if(velocityProfile != null)
        {
            if(debug) System.out.println("Interpolated velocityProfile at row " + row + " col " + col);
            return velocityProfile;
        }
        else
        {
            if(debug) System.out.println("Initialized first velocityProfile at row " + row + " col " + col);
            return computeInitialVelocityProfile(row, col);
        }
    }


    public boolean isActualVelocityProfile(int row, int col)
    {
        return interpolation.isDataPoint(row, col);
    }

    public StsVelocityProfile getActualVelocityProfile(int row, int col)
    {
        if(!interpolation.isDataPoint(row, col)) return null;
        return currentVelocityProfiles[row][col];
    }

    public StsVelocityProfile computeInterpolatedVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile;

        if(velocityProfiles.getSize() == 0) return null;
        if(!checkRunInterpolation()) return null;

        StsRadialInterpolation.Weights dataWeights = interpolation.getWeights(row, col);
        if(dataWeights == null) return null;

        int nWeights = dataWeights.nWeights;
        if(nWeights == 0) return null;
        Object[] dataObjects = dataWeights.dataObjects;
        // if we aren't interpolating (row-col same) just set this profile in the currentProfiles 2d array
        // if not the same row-col, but a single weight (only one profile exists), copy it and set in array
        if(nWeights == 1)
        {
            StsVelocityProfile otherProfile = (StsVelocityProfile) dataObjects[0];
            if(otherProfile.row == row && otherProfile.col == col)
                velocityProfile = otherProfile;
            else
                velocityProfile = new StsVelocityProfile(row, col, otherProfile);
        }
        else
        {
            double[] objectWeights = dataWeights.weights;
            TreeSet sortedTreeSet = new TreeSet(StsPoint.getComparator(1));
            StsVelocityProfile[] profiles = new StsVelocityProfile[nWeights];
            int nGoodProfiles = 0;
            double[] weights = new double[nWeights];
            for(int n = 0; n < nWeights; n++)
            {
                StsVelocityProfile profile = (StsVelocityProfile) dataObjects[n];
                StsPoint[] points = profile.getProfilePoints();
                if(points == null) continue;
                int nPoints = points.length;
                if(points.length == 0) continue;
                profiles[nGoodProfiles] = profile;
                weights[nGoodProfiles] = objectWeights[n];
                nGoodProfiles++;
                StsPoint[] newPoints = new StsPoint[nPoints];
                System.arraycopy(points, 0, newPoints, 0, nPoints);
                for(int p = 0; p < nPoints; p++)
                    sortedTreeSet.add(newPoints[p]);
            }
            if(nGoodProfiles < nWeights)
            {
                profiles = (StsVelocityProfile[]) StsMath.trimArray(profiles, nGoodProfiles);
                weights = (double[]) StsMath.trimArray(weights, nGoodProfiles);
                nWeights = nGoodProfiles;
            }
            for(int w = 0; w < nWeights; w++)
                weights[w] /= dataWeights.weightSum;

            int nSortedPoints = sortedTreeSet.size();
            if(nSortedPoints == 0) return null;
            StsPoint[] sortedPoints = new StsPoint[nSortedPoints];
            Iterator iter = sortedTreeSet.iterator();
            int p = 0;
            while (iter.hasNext())
                sortedPoints[p++] = new StsPoint((StsPoint) iter.next());
            for(p = 0; p < nSortedPoints; p++)
            {
                float v = 0.0f;
                float t = sortedPoints[p].v[1];
                for(int w = 0; w < nWeights; w++)
                {
                    float vv = profiles[w].getVelocityFromTzero(t);
                    v += weights[w] * vv;
                }
                //            v /= dataWeights.weightSum;
                if(Float.isNaN(v))
                {
                    //				System.out.println("way bad juju");
                    return null;
                }
                sortedPoints[p].v[0] = v;
            }
            velocityProfile = new StsVelocityProfile(row, col);
            velocityProfile.setProfilePoints(sortedPoints);
            computeInterpolatedMutes(velocityProfile, profiles, weights);
        }
        currentVelocityProfiles[row][col] = velocityProfile;
        return velocityProfile;
    }

    private StsVelocityProfile computeInitialVelocityProfile(int row, int col)
    {
        StsVelocityProfile velocityProfile = new StsVelocityProfile(row, col);
        velocityProfile.initializeMutes(lineSet);
        if(debug) System.out.println("Creating first velocity profile " + velocityProfile.toString());
        if(!checkInitializeInterpolation()) return null;
        currentVelocityProfiles[row][col] = velocityProfile;
        return velocityProfile;
        /*
            currentVelocityProfile = new StsVelocityProfile(row, col);
            currentVelocityProfile.initializeMutes(lineSet);
            if(debug) System.out.println("Creating first velocity profile " + currentVelocityProfile.toString());
            if(!checkInitializeInterpolation()) return null;
            currentVelocityProfiles[row][col] = currentVelocityProfile;
            return currentVelocityProfile;
        */
    }

    public boolean computeInterpolatedMutes(StsVelocityProfile interpolatedProfile)
    {
        if(!checkRunInterpolation()) return false;
        int row = interpolatedProfile.row;
        int col = interpolatedProfile.col;
        if(interpolation.isDataPoint(row, col)) return true;
        StsRadialInterpolation.Weights dataWeights = interpolation.getWeights(row, col);
        if(dataWeights == null) return false;
        int nWeights = dataWeights.nWeights;
        if(nWeights == 0) return false;
        double[] weights = dataWeights.weights;
        Object[] dataObjects = dataWeights.dataObjects;
        StsVelocityProfile[] profiles = new StsVelocityProfile[nWeights];
        for(int n = 0; n < nWeights; n++)
            profiles[n] = (StsVelocityProfile) dataObjects[n];
        computeInterpolatedMutes(interpolatedProfile, profiles, weights);
        return true;
    }

    /** interpolate mute from other profiles with given weights.  We are assuming sum of weights == 1.0 */
    public boolean computeInterpolatedMutes(StsVelocityProfile interpolatedProfile, StsVelocityProfile[] profiles, double[] weights)
    {
        int nWeights = weights.length;
        StsPoint topMute = null, bottomMute = null;
        for(int w = 0; w < nWeights; w++)
        {
            double weight = weights[w];

            StsPoint profileTopMute = profiles[w].topMute;
            if(profileTopMute != null)
            {
                if(topMute == null)
                {
                    topMute = new StsPoint(profileTopMute);
                    topMute.multiply(weight);
                }
                else
                {
                    topMute.addScaledVector(profileTopMute, (float) weight);
                }
            }
            StsPoint profileBottomMute = profiles[w].bottomMute;
            if(profileBottomMute != null)
            {
                if(bottomMute == null)
                {
                    bottomMute = new StsPoint(profileBottomMute);
                    bottomMute.multiply(weight);
                }
                else
                {
                    bottomMute.addScaledVector(profileBottomMute, (float) weight);
                }
            }
        }
        interpolatedProfile.topMute = topMute;
        interpolatedProfile.bottomMute = bottomMute;
        return true;
    }

    public boolean computeInterpolatedByteVelocityProfile(int row, int col, byte[] byteVelocities)
    {
        float[] floatData = computeInterpolatedFloatVelocityProfile(row, col, nSlices);
        double scale = 254 / (dataMax - dataMin);
        for(int n = 0; n < nSlices; n++)
            byteVelocities[n] = StsMath.unsignedIntToUnsignedByte((int) ((floatData[n] - dataMin) * scale));
        return true;
    }

    public boolean scaleVelocityTrace(float[] floatData, byte[] byteVelocities)
    {
        double scale = 254 / (dataMax - dataMin);
        for(int n = 0; n < floatData.length; n++)
            byteVelocities[n] = StsMath.unsignedIntToUnsignedByte((int) ((floatData[n] - dataMin) * scale));
        return true;
    }

    public float[] computeInterpolatedFloatVelocityProfile(int row, int col, int nSlices)
    {
        StsVelocityProfile interpolatedProfile = getComputeVelocityProfile(row, col);
        if(interpolatedProfile == null) return null;
        StsPoint[] points = interpolatedProfile.getProfilePoints();
        if(points == null || points.length < 2) return null;

        float volumeTMin = lineSet.getZMin();
        float volumeTInc = lineSet.getZInc();
        float t = volumeTMin;

        float[] vz0 = points[0].v;
        float[] vz1 = points[1].v;
        int nNextPoint = 2;
        int nPoints = points.length;
        float v0 = vz0[0];
        float v1 = vz1[0];
        int n = 0;
        float[] velocities = new float[nSlices];
        // rms velocity is constant down to the first point and is
        // linearly extrapolated beyond the last point
        for(; n < nSlices; n++, t += volumeTInc)
        {
            if(t > vz0[1]) break;
            velocities[n] = vz0[0];
        }
        for(; n < nSlices; n++, t += volumeTInc)
        {
            while (t > vz1[1] && nNextPoint < nPoints)
            {
                vz0 = vz1;
                vz1 = points[nNextPoint++].v;
            }
            v0 = v1;
            float f = (t - vz0[1]) / (vz1[1] - vz0[1]);
            v1 = vz0[0] + f * (vz1[0] - vz0[0]);
            if(velocityType == StsParameters.SAMPLE_TYPE_VEL_RMS)
                velocities[n] = v1;
            else if(velocityType == StsParameters.SAMPLE_TYPE_VEL_INTERVAL)
                // velocities[n] = (float)Math.sqrt((v1 * v1 * (t + volumeTInc) - v0 * v0 * t) / volumeTInc); //this code lets interval velocity increase with depth
                velocities[n] = (float) Math.sqrt((vz1[0] * vz1[0] * vz1[1] - vz0[0] * vz0[0] * vz0[1]) / (vz1[1] - vz0[1])); //this code holds interval velocity constant between picks
        }
        return velocities;
    }

    public boolean computeInterpolatedByteVelocityProfile(int row, int col, float velocityMin, float velocityMax, byte[] velocities, int nSlices, int position)
    {
        try
        {
            float[] floatVelocities = this.computeInterpolatedFloatVelocityProfile(row, col, nSlices);
            if(floatVelocities == null)
                System.arraycopy(StsPreStackLineSet.byteTransparentTrace, 0, velocities, position, nSlices);
            else
            {
                double scale = 254 / (velocityMax - velocityMin);
                for(int n = 0; n < nSlices; n++, position++)
                    velocities[position] = StsMath.unsignedIntToUnsignedByte((int) ((floatVelocities[n] - velocityMin) * scale));
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackVelocityModel3d.computeInterpolatedByteVelocityProfile() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean computeInterpolatedByteBufferVelocityProfile(int row, int col, float velocityMin, float velocityMax, ByteBuffer buffer, int nSlices, boolean neighborsOnly)
    {
        try
        {
            if(neighborsOnly && !interpolation.isNeighbor[row][col])
            {
                buffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nSlices);
                return true;
            }
            float[] floatVelocities = this.computeInterpolatedFloatVelocityProfile(row, col, nSlices);

            if(floatVelocities == null)
                buffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nSlices);
            else
            {
                double scale = 254 / (velocityMax - velocityMin);
                for(int n = 0; n < nSlices; n++)
                {
                    byte velocity = StsMath.unsignedIntToUnsignedByte((int) ((floatVelocities[n] - velocityMin) * scale));
                    buffer.put(velocity);
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackVelocityModel.computeInterpolatedByteVelocityProfile() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean computeInterpolatedFloatBufferVelocityProfile(int row, int col, FloatBuffer velocities, int nSlices)
    {
        float[] traceVelocities = computeInterpolatedFloatVelocityProfile(row, col, nSlices);
        if(traceVelocities == null) return false;
        velocities.put(traceVelocities);
        return true;
    }

    /** prestackVelocityModel is not currently using external files for data storage; on export, row planes are computed on the fly */
    public boolean setupReadRowFloatBlocks()
    {
        return true;
    }

    public boolean setupRowByteBlocks()
    {
        return true;
    }

    public void debugBuildProfiles()
    {
        Iterator rowColIterator = lineSet.getRowColIterator(1);
        while (rowColIterator.hasNext())
        {
            int[] rowCol = (int[]) rowColIterator.next();
            if(rowCol == null) return;
            int row = rowCol[0];
            int col = rowCol[1];
            if(!lineSet.isNTracesBelowThreshold(row, col))
                buildDebugProfile(row, col);
        }
        velocityProfiles.sort();
        checkInitializeInterpolation();
    }

    private void buildDebugProfile(int row, int col)
    {
        StsVelocityProfile profile = this.getComputeVelocityProfile(row, col);
        if(profile == null) return;
        if(lineSet.isNTracesBelowThreshold(row, col)) return;
        if(profile.isInterpolated())
        {
            currentModel.add(profile);
            velocityProfiles.add(profile);
            currentVelocityProfiles[row][col] = profile;
        }
    }

    public void setSeismicVelocityVolume(StsSeismic seismicVolume)
    {
        this.seismicVelocityVolume = seismicVolume;
        dbFieldChanged("seismicVelocityVolume", seismicVelocityVolume);
    }

    public StsSeismic getSeismicVelocityVolume()
    {
        return seismicVelocityVolume;
    }

    public boolean getDisplayAnalysisPoints()
    {
        StsPreStackVelocityModelClass velocityModelClass = (StsPreStackVelocityModelClass) getStsClass();
        return velocityModelClass.getDisplayAnalysisPoints();
    }

    class ExportPanelTabbedPane extends JTabbedPane
    {
        StsDialogFace[] panels = null;

        ExportPanelTabbedPane()
        {
        }

        public Component add(String title, Component component)
        {
            return super.add(title, component);
        }

        StsDialogFace getDialogFace()
        {
            return (StsDialogFace) getSelectedComponent();
        }
    }

    /**
     * set minimum and maximum velocity for display purposes
     * (resets colorbar)
     *
     * @param velocityMax
     * @param velocityMin
     */
    public void setVelocityMaxMin(float velocityMax, float velocityMin)
    {
        dataMin = velocityMin;
        dataMax = velocityMax;
        setColorScaleRange(velocityMin, velocityMax, true);
    }

    public void setColorScaleRange(float velocityMin, float velocityMax, boolean changed)
    {
        if(currentColorList != null)
        {
            currentColorList.colorscale.setRange(velocityMin, velocityMax);
            if(changed) currentColorList.colorscale.colorsChanged();
        }
    }

    public float[] computeMinMax()
    {
        float min = Float.MAX_VALUE;
        float max = 0 - min;
        if(velocityProfiles == null) return null;
        int nProfiles = velocityProfiles.getSize();
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            if(profile != null)
            {
                if(velocityType == StsParameters.SAMPLE_TYPE_VEL_INTERVAL)
                    profile = profile.getIntervalVelocityProfile();
                StsPoint[] points = profile.getProfilePoints();
                for(int j = 0; j < points.length; j++)
                {
                    min = Math.min(min, points[j].v[0]);
                    max = Math.max(max, points[j].v[0]);
                }
            }
        }
        return new float[]{min, max};
    }

    public abstract StsFieldBean[] getDisplayFields();

    public boolean computeInterpolatedFloatVelocityProfile(int row, int col, float[] velocities, int position, int nSlices)
    {
        float[] traceVelocities = computeInterpolatedFloatVelocityProfile(row, col, nSlices);
        if(traceVelocities == null) return false;
        System.arraycopy(traceVelocities, 0, velocities, position, nSlices);
        return true;
    }

    public float[] computeFloatVelocityPlane(int dir, int nPlane)
    {
        int row = -1, col = -1;
        float[] velocities = null;
        ByteBuffer buffer = null;

        if(dir == ZDIR) return null;
        StsPreStackLineSet.checkTransparentTrace(nSlices);
        // nCroppedSlices = lineSet.getNSlices(); // temporary!  remove when new db built!!
        StsPreStackLineSetClass lineSetClass = lineSet.lineSetClass;
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
//      float velocityMin = rangeProperties.velocityMin;
//      float velocityMax = rangeProperties.velocityMax;
        int position = 0;
        if(dir == XDIR)
        {
            velocities = new float[nRows * nSlices];
            col = nPlane;
            for(row = 0; row < nRows; row++, position += nSlices)
                computeInterpolatedFloatVelocityProfile(row, col, velocities, position, nSlices);
        }
        else if(dir == YDIR)
        {
            velocities = new float[nCols * nSlices];
            row = nPlane;
            for(col = 0; col < nCols; col++, position += nSlices)
                computeInterpolatedFloatVelocityProfile(row, col, velocities, position, nSlices);
        }
        return velocities;
    }
}