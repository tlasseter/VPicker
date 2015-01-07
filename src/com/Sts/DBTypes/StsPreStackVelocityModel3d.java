package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.nio.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

public class StsPreStackVelocityModel3d extends StsPreStackVelocityModel implements Cloneable, Serializable, StsTreeObjectI
{
//    transient public StsSpiralRadialInterpolation volumeInterpolation = null;
    /** user specified picking mode: create or edit
		 boolean pickMode = CREATE; */
    /** convenience copy of distance increment between rows */
//	float yInc;
	/** convenience copy of distance increment between cols */

/*
        static public final boolean CREATE = false;
        static public final boolean EDIT = true;
     */

	static public StsFieldBean[] displayFields = null;

    static StsComboBoxFieldBean colorscalesBean;
	static StsEditableColorscaleFieldBean colorscaleBean;
	static public StsFieldBean[] propertyFields = null;

	static public StsObjectPanel objectPanel = null;

    public String getGroupname()
    {
        return StsSeismicBoundingBox.group3d;
    }

    public StsPreStackVelocityModel3d()
	{
	}

	public StsPreStackVelocityModel3d(boolean persistent)
	{
		super(persistent);
	}

	public StsPreStackVelocityModel3d(StsPreStackLineSet lineSet, String name)
	{
		super(lineSet, name);
	}

	public StsPreStackVelocityModel3d(StsPreStackLineSet lineSet)
	{
        super(lineSet);
	}

    static public void setStsObjectCopierCloneFlags(StsObjectCopier copier)
	{
		copier.setFieldCloneStsObject("lineSet");
	}


//	public boolean isComplete() { return isComplete; };
//	public void setIsComplete(boolean isComplete) { this.isComplete = isComplete; }

    public float[] getVolumeVelocities(int row, int col)
	{
		if(seismicVelocityVolume == null)return null;
		int nVolumeSlices = lineSet.getNSlices();
		float volumeZMin = lineSet.getZMin();
		float volumeZInc = lineSet.getZInc();
		int nVelocitySlices = seismicVelocityVolume.getNSlices();
		float[] velocityVolumeVelocities = new float[nVelocitySlices];
		if(!seismicVelocityVolume.getTraceValues(row, col, 0, nVelocitySlices - 1, YDIR, true, velocityVolumeVelocities))return null;
		float velocityZMin = seismicVelocityVolume.getZMin();
		float velocityZInc = seismicVelocityVolume.getZInc();
		float[] velocities = new float[nVolumeSlices];
		float z = volumeZMin;
		for(int n = 0; n < nVolumeSlices; n++, z += volumeZInc)
		{
			float indexF = (z - velocityZMin) / velocityZInc;
			int i = StsMath.minMax((int)indexF, 0, nVelocitySlices - 2);
			float f = indexF - i;
			float v0 = velocityVolumeVelocities[i];
			float v1 = velocityVolumeVelocities[i + 1];
			float v = v0 + f * (v1 - v0);
			velocities[n] = v * multiplier;
		}
		return velocities;
	}

	public float[][] getVolumeVelocityProfile(int row, int col)
	{
		if(seismicVelocityVolume == null)return null;
		int nVolumeSlices = lineSet.getNSlices();
		float volumeZMin = lineSet.getZMin();
		float volumeZInc = lineSet.getZInc();
		int nVelocitySlices = seismicVelocityVolume.getNSlices();
		float[] velocityVolumeVelocities = new float[nVelocitySlices];
		if(!seismicVelocityVolume.getTraceValues(row, col, 0, nVelocitySlices - 1, YDIR, true, velocityVolumeVelocities))return null;
		float velocityZMin = seismicVelocityVolume.getZMin();
		float velocityZInc = seismicVelocityVolume.getZInc();
		float[][] velocities = new float[nVolumeSlices][2];
		float z = volumeZMin;
		for(int n = 0; n < nVolumeSlices; n++, z += volumeZInc)
		{
			float indexF = (z - velocityZMin) / velocityZInc;
			int i = StsMath.minMax((int)indexF, 0, nVelocitySlices - 2);
			float f = indexF - i;
			float v0 = velocityVolumeVelocities[i];
			float v1 = velocityVolumeVelocities[i + 1];
			float v = v0 + f * (v1 - v0);
			velocities[n][0] = v * multiplier;
			velocities[n][1] = z;
		}
		return velocities;
	}
/*
    public StsVelocityProfile getCreateCurrentVelocityProfile(float rowNum, float colNum)
    {
        int row = lineSet.getRowFromRowNum(rowNum);
        int col = lineSet.getColFromColNum(colNum);
        StsVelocityProfile velocityProfile = getComputeVelocityProfile(row, col);
        if(velocityProfile == null) return null;

        currentVelocityProfile = velocityProfile;
        if(debug)
        {
            if(velocityProfile.isInterpolated())
                System.out.println("Using interpolated profile " + currentVelocityProfile.toString());
            else
                System.out.println("Using existing velocity profile " + velocityProfile.toString());
        }
        return currentVelocityProfile;
    }
*/
    public void updateVelocityProfile(StsVelocityProfile adjustedProfile, int indexPicked)
    {
        if(adjustedProfile == null)return;
        if(adjustedProfile.isInterpolated()) return;
        super.updateVelocityProfile(adjustedProfile, indexPicked);
        currentModel.clearDisplayTextured3dCursors(lineSet);
    }

	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null) return displayFields;
        velocityTypeBean = new StsComboBoxFieldBean(StsPreStackVelocityModel3d.class, "velocityTypeString", "Velocity type", StsParameters.VEL_STRINGS);
	    displayFields = new StsFieldBean[]
        {
 //           new StsBooleanFieldBean(StsPreStackVelocityModel3d.class, "isVisible", "Visible"),
            velocityTypeBean
        };
        return displayFields;
	}

	public StsFieldBean[] getPropertyFields()
	{
        if(propertyFields != null) return propertyFields;
        colorscalesBean = new StsComboBoxFieldBean(StsPreStackVelocityModel3d.class, "currentColorscale", "Colorscales");
	    colorscaleBean = new StsEditableColorscaleFieldBean(StsPreStackVelocityModel3d.class, "currentColorscale");

	    propertyFields = new StsFieldBean[]
        {
            new StsFloatFieldBean(StsPreStackVelocityModel3d.class, "dataMin", false, "Velocity Min:"),
            new StsFloatFieldBean(StsPreStackVelocityModel3d.class, "dataMax", false, "Velocity Max:"),
            new StsStringFieldBean(StsPreStackVelocityModel3d.class, "velocityUnitsString", false, "Input Velocity Units:"),
            colorscalesBean,
            colorscaleBean
        };
		colorscalesBean.setListItems(velocityColorscales.getTrimmedList());
		colorscalesBean.setSelectedItem(currentColorscale);

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
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		 tabbedPane = new ExportPanelTabbedPane();

		 StsSeismicExportPanel volExportPanel = new StsSeismicExportPanel3d(currentModel, this, "Export Velocities as SEGY", true);
		 tabbedPane.add("Velocity PostStack3d", volExportPanel);
		 StsPreStackHandVelExportPanel handVelExportPanel = new StsPreStackHandVelExportPanel3d(currentModel, this);
		 tabbedPane.add("Hand Vels", handVelExportPanel);
		 return tabbedPane;
	}

	public void treeObjectSelected()
	{
		StsClass velocityVolumeClass = currentModel.getStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d");
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

    public float[] computeInterpolatedDepthFloatVelocityProfile(int row, int col)
	{
        StsVelocityProfile interpolatedProfile = getComputeVelocityProfile(row, col);
        if(interpolatedProfile == null) return null;
        StsPoint[] points = interpolatedProfile.getProfilePoints();
		if(points == null)return null;

		float volumeTMin = lineSet.getZMin();
		float volumeTInc = lineSet.getZInc();
		float t = volumeTMin;

		float[] vz0 = points[0].v;
		float[] vz1 = points[1].v;
		int nNextPoint = 2;
		int nPoints = points.length;
		float v1 = vz1[0];
		float v0 = v1;
		int n = 0;
		float[] velocities = new float[nSlices];
		// rms velocity is constant down to the first point and is
		// linearly extrapolated beyond the last point
		for(; n < nSlices; n++, t += volumeTInc)
		{
			if(t > vz1[1])break;
			velocities[n] = v0;
		}
		for(; n < nSlices; n++, t += volumeTInc)
		{
			while(t > vz1[1] && nNextPoint < nPoints)
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
				velocities[n] = (float)Math.sqrt((v1 * v1 * (t + volumeTInc) - v0 * v0 * t) / volumeTInc);
		}
		return velocities;
	}

    /*
        public boolean computeInterpolatedVelocityProfilesPoints()
        {
            if(!checkInitializeInterpolation()) return false;
            if(interpolatedVelocityPoints != null)return true;
            interpolatedVelocityPoints = new StsPoint[nRows][nCols][];
            for(int row = 0; row < nRows; row++)
                for(int col = 0; col < nCols; col++)
                    interpolatedVelocityPoints[row][col] = computeInterpolatedVelocityProfilePoints(row, col);
            return true;
        }
    */
    public float[] readRowPlaneFloatData(int nPlane)
	{
		return computeFloatVelocityPlane(YDIR, nPlane);
	}

	public byte[] readRowPlaneByteData(int nPlane)
	{
		return computeByteVelocityPlane(YDIR, nPlane);
	}

	public FloatBuffer readRowPlaneFloatBufferData(int nPlane)
	{
		return computeFloatBufferVelocityPlane(YDIR, nPlane);
	}

	public ByteBuffer readRowPlaneByteBufferData(int nPlane)
	{
		return computeByteBufferVelocityPlane(YDIR, nPlane, false);
	}

	public byte[] computeByteVelocityPlane(int dir, int nPlane)
	{
		int row = -1, col = -1;
		byte[] velocities = null;
		ByteBuffer buffer = null;

		if(dir == ZDIR)return null;
		StsPreStackLineSet.checkTransparentTrace(nSlices);
//		nCroppedSlices = lineSet.getNSlices(); // temporary!  remove when new db built!!
//		StsPreStackLineSet3dClass lineSet = line2d.lineSet;
//		StsSemblanceRangeProperties rangeProperties = line2d.semblanceRangeProperties;
//		float velocityMin = rangeProperties.velocityMin;
//		float velocityMax = rangeProperties.velocityMax;
		int position = 0;
		if(dir == XDIR)
		{
			velocities = new byte[nRows * nSlices];
			col = nPlane;
			for(row = 0; row < nRows; row++, position += nSlices)
				computeInterpolatedByteVelocityProfile(row, col, dataMin, dataMax, velocities, nSlices, position);
		}
		else if(dir == YDIR)
		{
			velocities = new byte[nCols * nSlices];
			row = nPlane;
			for(col = 0; col < nCols; col++, position += nSlices)
				computeInterpolatedByteVelocityProfile(row, col, dataMin, dataMax, velocities, nSlices, position);
		}
		return velocities;
	}

	public ByteBuffer computeByteBufferVelocityPlane(int dir, int nPlane, boolean neighborsOnly)
	{
		int row = -1, col = -1;
		byte[] velocities = null;
		ByteBuffer buffer = null;

		if(dir == ZDIR)return null;
		StsPreStackLineSet.checkTransparentTrace(nSlices);
//		nCroppedSlices = lineSet.getNSlices(); // temporary!  remove when new db built!!
		StsPreStackLineSetClass lineSetClass = lineSet.lineSetClass;
		StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
//		float velocityMin = rangeProperties.velocityMin;
//		float velocityMax = rangeProperties.velocityMax;
		velocities = new byte[nSlices];
		int position = 0;
		if(dir == XDIR)
		{
			buffer = ByteBuffer.allocateDirect(nRows * nSlices);
			col = nPlane;
			for(row = 0; row < nRows; row++, position++)
				computeInterpolatedByteBufferVelocityProfile(row, col, dataMin, dataMax, buffer, nSlices, neighborsOnly);
		}
		else if(dir == YDIR)
		{
			buffer = ByteBuffer.allocateDirect(nCols * nSlices);
			row = nPlane;
			for(col = 0; col < nCols; col++)
				computeInterpolatedByteBufferVelocityProfile(row, col, dataMin, dataMax, buffer, nSlices, neighborsOnly);
		}
		buffer.rewind();

        // Filter buffer

		return buffer;
	}

	public ByteBuffer computeByteBufferVelocityArbLine(int[] rows, int[]cols, boolean neighborsOnly)
	{
		int row = -1, col = -1;
		byte[] velocities = null;
		ByteBuffer buffer = null;


		StsPreStackLineSet.checkTransparentTrace(nSlices);
//		nCroppedSlices = lineSet.getNSlices(); // temporary!  remove when new db built!!
		StsPreStackLineSetClass lineSetClass = lineSet.lineSetClass;
		StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
		float velocityMin = rangeProperties.velocityMin;
		float velocityMax = rangeProperties.velocityMax;
		velocities = new byte[nSlices];
		int position = 0;
		StsStatusArea statusArea = currentModel.win3d.statusArea;
		statusArea.setStatus("Computing velocity display...");
		statusArea.addProgress();

		statusArea.setMaximum(rows.length);

		{
			buffer = ByteBuffer.allocateDirect(rows.length*nSlices);

			for(int i = 0; i < rows.length; i++)
			{
				row = rows[i];
				col = cols[i];
				computeInterpolatedByteBufferVelocityProfile(row, col, velocityMin, velocityMax, buffer, nSlices, neighborsOnly);
				statusArea.setProgress(i);
			}
		}

		buffer.rewind();
		statusArea.removeProgress();
		statusArea.clearStatus();

		return buffer;
	}

	public FloatBuffer computeFloatBufferVelocityPlane(int dir, int nPlane)
	{
		int row = -1, col = -1;
		float[] velocities = null;
		FloatBuffer buffer = null;

		if(dir == ZDIR)return null;
		StsPreStackLineSet.checkTransparentTrace(nSlices);
//		nCroppedSlices = lineSet.getNSlices(); // temporary!  remove when new db built!!
//		StsPreStackLineSet3dClass lineSet = line2d.lineSet;
//		StsSemblanceRangeProperties rangeProperties = line2d.semblanceRangeProperties;
//		float velocityMin = rangeProperties.velocityMin;
//		float velocityMax = rangeProperties.velocityMax;
		if(dir == XDIR)
		{
			buffer = FloatBuffer.allocate(nRows * nSlices);
			col = nPlane;
			for(row = 0; row < nRows; row++)
				computeInterpolatedFloatBufferVelocityProfile(row, col, buffer, nSlices);
		}
		else if(dir == YDIR)
		{
			buffer = FloatBuffer.allocate(nCols * nSlices);
			row = nPlane;
			for(col = 0; col < nCols; col++)
				computeInterpolatedFloatBufferVelocityProfile(row, col, buffer, nSlices);
		}
		return buffer;
	}

    public boolean initialize()
	{
		xInc = lineSet.getXInc();
		yInc = lineSet.getYInc();
//		currentX = lineSet.getXMin();
//		currentY = lineSet.getYMin();
		StsWindowFamily windowFamily = currentModel.getWindowFamily(currentModel.win3d);
//		windowFamily.adjustCursor(YDIR, currentY);
//		windowFamily.adjustCursor(XDIR, currentX);
        return true;
	}

    protected void constructInterpolatedVelocityPoints()
    {
        currentVelocityProfiles = new StsVelocityProfile[nRows][nCols];
    }
    
    @Override
    public boolean initialize(StsModel model)
    {
        if (!super.initialize(model)) return false;
        StsWindowFamily windowFamily = currentModel.getWindowFamily(currentModel.win3d);
 //       windowFamily.adjustCursor(YDIR, currentY);
 //       windowFamily.adjustCursor(XDIR, currentX);
        return true;
    }
}
