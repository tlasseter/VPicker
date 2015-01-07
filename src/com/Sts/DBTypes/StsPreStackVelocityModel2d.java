package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
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
public class StsPreStackVelocityModel2d extends StsPreStackVelocityModel implements Cloneable, StsTreeObjectI, StsViewable
{
    /** velocity data is stored in a floatBuffer and accessed/rebuilt as needed */
    transient FloatBuffer velocityFloatBuffer;
    /** For each gather, if true, an actual profile is in the neighborhood */
    transient public boolean[][] isNeighbor;

	static public StsFieldBean[] displayFields = null;

    static StsComboBoxFieldBean colorscalesBean;
	static StsEditableColorscaleFieldBean colorscaleBean;
	static public StsFieldBean[] propertyFields = null;

	static public StsObjectPanel objectPanel = null;

    public String getGroupname()
    {
        return StsSeismicBoundingBox.group2d;
    }

    public StsPreStackVelocityModel2d()
	{
	}

     public StsPreStackVelocityModel2d(boolean persistent)
	{
		super(persistent);
	}

	public StsPreStackVelocityModel2d(StsPreStackLineSet lineSet, String name)
	{
		super(lineSet, name);
        initialize();
    }

	public StsPreStackVelocityModel2d(StsPreStackLineSet lineSet)
	{
        super(lineSet);
        initialize();
    }

    public boolean initialize()
    {
        return true;
    }

    public int getColumnsInRow(int row)
	{
		return lineSet.lines[row].nCols;
	}

    public double[][] getLineXYCDPs(int row)
	{
        return lineSet.getLineXYCDPs(row);
	}

	public int getNSlices(int row)
	{
		return lineSet.getNSlices(row);
	}

    public int getNCols(int row)
    {
        return lineSet.getNColsForRow(row);
    }

    public float getZInc(int row)
    {
        return lineSet.getZInc(row);
    }

	public StsSeismicBoundingBox getLineBoundingBox(int row)
	{
		return lineSet.getLineBoundingBox(row);
	}

    protected void constructInterpolatedVelocityPoints()
    {
        StsPreStackLine[] lines = lineSet.lines;
        currentVelocityProfiles = new StsVelocityProfile[nRows][];
        for(int row = 0; row < nRows; row++)
            currentVelocityProfiles[row] = new StsVelocityProfile[lines[row].nCols];
    }

    public float[] getVolumeVelocities(int row, int col)
    {
        return null;
    }

    public float[][] getVolumeVelocityProfile(int row, int col)
    {
        return null;
    }
/*
   public void setCurrentVelocityProfile(int row, int col)
    {
        if(currentVelocityProfile != null)
        {
            if(currentVelocityProfile.row == row && currentVelocityProfile.col == col)
            {
                if(debug) System.out.println("Continuing use of velocity profile " + currentVelocityProfile.toString());
                return;
            }
            checkCurrentVelocityProfile();
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
    private void initializeNeighborhood()
    {
        isNeighbor = new boolean[nRows][];
        for(int row = 0; row < nRows; row++)
        {
            int nGathers = lineSet.lines[row].nGathers;
            isNeighbor[row] = new boolean[nGathers];
        }
        float neighborRadius = lineSet.lineSetClass.getStackNeighborRadius();
        float radiusSq = neighborRadius*neighborRadius;
        int nProfiles = velocityProfiles.getSize();
        for(int n = 0; n < nProfiles; n++)
                addProfileToNeighborhood((StsVelocityProfile)velocityProfiles.getElement(n), radiusSq);
    }

    private void  addProfileToNeighborhood(StsVelocityProfile profile, float radiusSq)
    {
        int profileRow = profile.row;
        int profileCol = profile.col;
        StsPreStackLine2d line = (StsPreStackLine2d)lineSet.lines[profileRow];
        int nLineCols = line.nGathers;
        float x0 = line.cdpX[profileCol];
        float y0 = line.cdpY[profileCol];
        isNeighbor[profileRow][profileCol] = true;
        for(int col = profileCol+1; col < nLineCols; col++)
        {
            float x = line.cdpX[col];
            float y = line.cdpY[col];
            float distSq = (x - x0)*(x - x0) + (y - y0)*(y - y0);
            if(distSq > radiusSq) break;
            isNeighbor[profileRow][col] = true;
        }
        for(int col = profileCol-1; col >= 0; col--)
        {
            float x = line.cdpX[col];
            float y = line.cdpY[col];
            float distSq = (x - x0)*(x - x0) + (y - y0)*(y - y0);
            if(distSq > radiusSq) break;
            isNeighbor[profileRow][col] = true;
        }
    }
*/
    /** TODO use different data structure so velocityProfiles are organized by lines */
    public StsVelocityProfile[] getVelocityProfiles2d(int nLine)
    {
        int nProfiles = velocityProfiles.getSize();
        StsVelocityProfile[] lineProfiles = new StsVelocityProfile[nProfiles];
        int nLineProfiles = 0;
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile)velocityProfiles.getElement(n);
            int profileLine = profile.row;
            if(profileLine > nLine) continue;
            if(profileLine == nLine)
                lineProfiles[nLineProfiles++] = profile;
        }
        return (StsVelocityProfile[]) StsMath.trimArray(lineProfiles, nLineProfiles);
    }
/*
    public void updateVelocityProfile(StsVelocityProfile adjustedProfile, int indexPicked, StsPoint[] originalPoints)
    {
        if(adjustedProfile == null)return;
        super.updateVelocityProfile(adjustedProfile, indexPicked, originalPoints);
        StsPreStackLine2d line = (StsPreStackLine2d)lineSet.lines[adjustedProfile.row];
        line.textureChanged();
        line.lineStackOk = false;
    }
*/
    public Buffer computeLineVelocityBuffer(StsPreStackLine2d line, boolean neighborsOnly, boolean outputFloats)
    {
        Buffer buffer = null;

        int row = line.lineIndex;
        int nGathers = line.nCols;

        int nLineSlices = line.nSlices;
        StsPreStackLineSet.checkTransparentTrace(nLineSlices);

        if (outputFloats)
        {
            buffer = ByteBuffer.allocateDirect(nGathers * nLineSlices * 4).asFloatBuffer();
            int bufferIndex = 0;
            for(int col = 0; col < nGathers; col++, bufferIndex += nLineSlices)
                computeInterpolatedFloatBufferVelocityProfile(row, col, (FloatBuffer)buffer, nLineSlices);
        }
        else
        {
            float floatMin = dataMin;
            float floatMax = dataMax;
            boolean autoSaturate = lineSet.semblanceRangeProperties.autoSaturateColors;
            float[] floatPlane = computeFloatVelocityPlane(row);
            if (floatPlane == null) return null;
            if (autoSaturate)
            {
                float[] minMax = StsMath.minMax(floatPlane);
                floatMin = minMax[0];
                floatMax = minMax[1];
                for (int i=0; i<lineSet.nRows; i++) //search min/max of other lines so all are isVisible to same scale
                {
                    if ( i != row)
                    {
                        float[] otherPlane = computeFloatVelocityPlane(i); 
                        if (otherPlane != null)
                        {
                            minMax = StsMath.minMax(otherPlane);
                            floatMin = Math.min(minMax[0], floatMin);
                            floatMax = Math.max(minMax[1], floatMax);
                        }
                    }
                }
                setColorScaleRange((float)dataMin, (float)dataMax, false); //this should update velocity colorbar to show new velocity ranges
            }

            double scale = 253.0/(floatMax - floatMin);  //using 254 was causing the max velocity to come out transparent - very obvious in interval-velocity mode SWC 11/17/09
            ByteBuffer buffer2 = ByteBuffer.allocateDirect(nGathers * nLineSlices);
            byte velocity = 0;
            for(int n = 0; n < nGathers * nLineSlices; n++)
            {
                velocity = StsMath.unsignedIntToUnsignedByte((int)((floatPlane[n] - floatMin) * scale));
                buffer2.put(velocity);
            }
            buffer2.rewind();
            return buffer2;
            
            /*
            buffer = ByteBuffer.allocateDirect(nGathers * nLineSlices);
            //int bufferIndex = 0;
            for(int col = 0; col < nGathers; col++)
                computeInterpolatedByteBufferVelocityProfile(row, col, floatMin, floatMax, (ByteBuffer)buffer, nLineSlices, neighborsOnly);
                */
        }
        buffer.rewind();
        return buffer;
    }
/*
	public float[] computeInterpolatedFloatVelocityProfile(int row, int col)
	{
		StsPoint[] points = this.getComputeVelocityProfile(row, col).profilePoints;
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
		float[] velocities = new float[nCroppedSlices];
		// rms velocity is constant down to the first point and is
		// linearly extrapolated beyond the last point
		for(; n < nCroppedSlices; n++, t += volumeTInc)
		{
			if(t > vz1[1])break;
			velocities[n] = v0;
		}
		for(; n < nCroppedSlices; n++, t += volumeTInc)
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
*/
    public ByteBuffer computeVelocity2dBytes(StsPreStackLine2d line)
    {
        if(!hasProfiles(line)) return null;
 //       StsVelocityProfile[] velocityProfiles = getVelocityProfiles2d(line.lineIndex);
        return (ByteBuffer)computeLineVelocityBuffer(line, false, false);
     }

    public FloatBuffer computeVelocity2dFloats(StsPreStackLine2d line)
    {
        if(!hasProfiles()) return null;
        //StsVelocityProfile[] velocityProfiles = getVelocityProfiles2d(line.lineIndex);
        return (FloatBuffer)computeLineVelocityBuffer(line, false, true);
    }

	 public float[] readRowPlaneFloatData(int nPlane)
	 {
		 return computeFloatVelocityPlane(nPlane);
	 }

	 public byte[] readRowPlaneByteData(int nPlane)
	 {
		 return computeByteVelocityPlane(nPlane);
	 }

	 public float[] computeFloatVelocityPlane(int nPlane)
	 {
		 int row = nPlane;
         int col = -1;
		 float[] velocities;
         int nLineSlices = lineSet.getNSlices(row);
         int nLineCols = lineSet.getNColsForRow(row);
         StsPreStackLineSet.checkTransparentTrace(nLineSlices);
		 int position = 0;
		 velocities = new float[nLineCols * nLineSlices];
		 for (col = 0; col < nLineCols; col++, position += nLineSlices)
			 if (!computeInterpolatedFloatVelocityProfile(row, col, velocities, position, nLineSlices)) return null;
		 return velocities;
	 }

	 public FloatBuffer computeFloatBufferVelocityPlane(int nPlane)
	 {
		 FloatBuffer buffer = null;
		 int row = nPlane;
         int col = -1;
         int nLineSlices = lineSet.getNSlices(row);
         int nLineCols = lineSet.getNColsForRow(row);
		 StsPreStackLineSet.checkTransparentTrace(nLineSlices);
		 buffer = FloatBuffer.allocate(nLineCols * nLineSlices);
		 for (col = 0; col < nLineCols; col++)
			 computeInterpolatedFloatBufferVelocityProfile(row, col, buffer, nSlices);
		 return buffer;
	 }

	 /*
	 private boolean computeInterpolatedFloatVelocityProfile(int row, int col, float[] velocities, int nLineSlices, int position)
	 {
		 float[] traceVelocities = computeInterpolatedFloatVelocityProfile(row, col, nLineSlices);
		 if (traceVelocities == null)
			 return false;
		 System.arraycopy(traceVelocities, 0, velocities, position, nLineSlices);
		 return true;
	 }
	 */

	public byte[] computeByteVelocityPlane(int nPlane)
	{
		int row = nPlane;
        int col = -1;
		byte[] velocities = null;
        int nLineSlices = lineSet.getNSlices(row);
        int nCols = lineSet.getNColsForRow(row);
        StsPreStackLineSet.checkTransparentTrace(nLineSlices);
		int position = 0;
		velocities = new byte[nCols * nLineSlices];
		for (col = 0; col < nCols; col++, position += nLineSlices)
			computeInterpolatedByteVelocityProfile(row, col, dataMin, dataMax, velocities, nLineSlices, position);
		return velocities;
	}

    // TODO a potentially expensive search since it's done every time we draw: reorganize profiles or make transient arrays for each row for 2D
    public boolean hasProfiles(StsPreStackLine2d line)
    {
        int nProfiles = velocityProfiles.getSize();
        if(nProfiles == 0) return false;
        int row = line.lineIndex;
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile)velocityProfiles.getElement(n);
            if(profile.row == row) return true;
        }
        return false;
    }

	 public boolean hasExportableData(int row)
	 {
		 if(!checkRunInterpolation()) return false;
		 return interpolation.rowHasDataPoint(row);
//         return	hasProfiles((StsPreStackLine2d)((StsPreStackLineSet2d)this.lineSet).lines[row]);
	 }
/*
     private Buffer computeLineVelocityByteBuffer(StsPreStackLine2d line, StsVelocityProfile[] velocityProfiles, boolean outputFloats)
     {
         ByteBuffer byteBuffer = null;
         FloatBuffer floatBuffer = null;
         byte[] byteVelocities = null;
         float[] floatVelocities = null;
         try
         {
             int nProfiles = velocityProfiles.length;
             if(nProfiles == 0) return null;
             int nGathers = line.nGathers;
             int nCroppedSlices = line.nCroppedSlices;

             StsPreStackLineSet.checkTransparentTrace(nCroppedSlices);

             floatVelocities = new float[nCroppedSlices];
             byteVelocities = new byte[nCroppedSlices];

             if(outputFloats)
                 floatBuffer = ByteBuffer.allocateDirect(nGathers * nCroppedSlices * 4).asFloatBuffer();
             else
                 byteBuffer = ByteBuffer.allocateDirect(nGathers * nCroppedSlices);

             if(nProfiles == 1)
             {
                 if(!computeSampleByteVelocityProfile(velocityProfiles[0], byteVelocities, floatVelocities))
                     return null;
                 for(int n = 0; n < nGathers; n++)
                     byteBuffer.put(byteVelocities);
             }
             else
             {
                 float[][] profileVelocities = new float[nProfiles][nCroppedSlices];
                 float[][] profileXYs = new float[nProfiles][2];
                 int[] profilePositions = new int[nProfiles];
                 for(int n = 0; n < nProfiles; n++)
                 {
                     int col =  velocityProfiles[n].col;
                     profilePositions[n] = col;
                     profileXYs[n][0] = line.cdpX[col];
                     profileXYs[n][1] = line.cdpY[col];
                     computedSampledFloatVelocityProfile(velocityProfiles[n], profileVelocities[n]);
                 }
                 int nextProfile = 0;
                 for(int n = 0; n < nGathers; n++)
                 {
                     if(nextProfile < nProfiles && n == profilePositions[nextProfile])
                     {
                         if(outputFloats)
                         {
                                floatBuffer.put(profileVelocities[nextProfile], 0, nCroppedSlices);
                         }
                         else
                         {
                             if(!computeSampleByteVelocityProfile(profileVelocities[nextProfile], byteVelocities))
                                 byteBuffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nCroppedSlices);
                             else
                                byteBuffer.put(byteVelocities);
                         }
                         nextProfile++;
                     }
                     else
                     {
                         float x = line.cdpX[n];
                         float y = line.cdpY[n];
                         float[] resultVelocities = interpolateVelocities(x, y, nProfiles, profileVelocities, profileXYs);
                         if(resultVelocities == null)
                         {
                            if(outputFloats)
                                floatBuffer.put(StsPreStackLineSet.floatTransparentTrace, 0, nCroppedSlices);
                             else
                                 byteBuffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nCroppedSlices);
                         }
                         else
                         {
                             if(outputFloats)
                             {
                                 floatBuffer.put(resultVelocities);
                             }
                             else
                             {
                                computeSampleByteVelocityProfile(resultVelocities, byteVelocities);
                                byteBuffer.put(byteVelocities, 0, nCroppedSlices);
                             }
                         }
                     }
                 }
             }
             byteBuffer.rewind();
             return byteBuffer;
         }
         catch(Exception e)
         {
             StsException.outputException("StsPreStackVelocityModel3d.computeInteprolatedByteVelocityProfile() failed.", e, StsException.WARNING);
             return null;
         }
     }

     private float[] interpolateVelocities(float x, float y, int nProfiles, float[][] profileVelocities, float[][] profileXYs)
     {
         int maxNPoints = 6;
         StsScatteredInterpolation interpolator = new StsScatteredInterpolation(maxNPoints);
         for(int n = 0; n < nProfiles; n++)
         {
             double dx = x - profileXYs[n][0];
             double dy = y - profileXYs[n][1];
             WeightedProfile weightedProfile = new WeightedProfile(profileVelocities[n], profileXYs[n][0], profileXYs[n][1], x, y);
             interpolator.addObject(weightedProfile);
         }
         interpolator.setRadialWeights();
         WeightedProfile resultProfile = new WeightedProfile(new float[nCroppedSlices]);
         resultProfile = (WeightedProfile)interpolator.interpolate(resultProfile);
         if(resultProfile == null) return null;
         return (float[])resultProfile.valueObject;
     }

     class WeightedProfile extends StsScatteredDataObject
     {

         WeightedProfile(float[] velocities, float x, float y, float x0, float y0)
         {
             super((x - x0)*(x - x0) + (y - y0)*(y - y0), velocities);
         }

         // constructor for a results profile
         WeightedProfile(float[] velocities)
         {
             super(0, velocities);
         }

         public StsScatteredDataObject interpolateValue(StsScatteredDataObject result)
         {
             WeightedProfile profileResult = (WeightedProfile)result;
             float[] resultVelocities = (float[])profileResult.valueObject;
             float[] velocities = (float[])valueObject;
             for(int n = 0; n < resultVelocities.length; n++)
                 resultVelocities[n] += velocities[n]*weight;
             return result;
         }
     }

     private boolean computeSampleByteVelocityProfile(StsVelocityProfile profile, byte[] byteVelocities, float[] floatVelocities)
     {
         try
         {
             if(!computedSampledFloatVelocityProfile(profile, floatVelocities)) return false;
             double scale = 254 / (dataMax - dataMin);
             for(int n = 0; n < nCroppedSlices; n++)
             {
                 byteVelocities[n] = StsMath.UnsignedIntToUnsignedByte((int)((floatVelocities[n] - dataMin) * scale));
             }
             return true;
         }
         catch(Exception e)
         {
             StsException.outputWarningException(this, "computeSampleByteVelocityProfile", e);
             return false;
         }
     }

     private boolean computeSampleByteVelocityProfile(float[] floatVelocities, byte[] byteVelocities)
     {
         try
         {
             double scale = 254 / (dataMax - dataMin);
             for(int n = 0; n < nCroppedSlices; n++)
             {
                 byteVelocities[n] = StsMath.UnsignedIntToUnsignedByte((int)((floatVelocities[n] - dataMin) * scale));
             }
             return true;
         }
         catch(Exception e)
         {
             StsException.outputWarningException(this, "computeSampleByteVelocityProfile", e);
             return false;
         }
     }

     private boolean computedSampledFloatVelocityProfile(StsVelocityProfile profile, float[] velocities)
     {
         StsPoint[] points = profile.getProfilePoints();
         if(points == null)return false;

         float[] vz0 = points[0].v;
         float[] vz1 = points[1].v;
         int nNextPoint = 2;
         int nPoints = points.length;
         float v1 = vz1[0];
         float v0 = v1;
         int n = 0;
         float z = zMin;
         // rms velocity is constant down to the first point and is
         // linearly extrapolated beyond the last point
         for(; n < nCroppedSlices; n++, z += zInc)
         {
             if(z > vz1[1])break;
             velocities[n] = v0;
         }
         for(; n < nCroppedSlices; n++, z += zInc)
         {
             while(z > vz1[1] && nNextPoint < nPoints)
             {
                 vz0 = vz1;
                 vz1 = points[nNextPoint++].v;
             }
             v0 = v1;
             float f = (z - vz0[1]) / (vz1[1] - vz0[1]);
             v1 = vz0[0] + f * (vz1[0] - vz0[0]);
             if(velocityType == StsParameters.SAMPLE_TYPE_VEL_RMS)
                 velocities[n] = v1;
             else if(velocityType == StsParameters.SAMPLE_TYPE_VEL_INTERVAL)
                 velocities[n] = (float)Math.sqrt((v1 * v1 * (z + zInc) - v0 * v0 * z) / zInc);
         }
         return true;
     }
*/

    public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null) return displayFields;
        velocityTypeBean = new StsComboBoxFieldBean(StsPreStackVelocityModel3d.class, "velocityTypeString", "Velocity type", StsParameters.VEL_STRINGS);
	    displayFields = new StsFieldBean[]
        {
//            new StsBooleanFieldBean(StsPreStackVelocityModel3d.class, "isVisible", "Visible"),
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

		 StsSeismicExportPanel volExportPanel = new StsSeismicExportPanel2d(currentModel, this, "Export Velocities as SEGY", true);
		 tabbedPane.add("Velocity PostStack3d", volExportPanel);
		 StsPreStackHandVelExportPanel handVelExportPanel = new StsPreStackHandVelExportPanel2d(currentModel, this);
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

    public Class[] getViewClasses() { return new Class[] {StsViewResidualSemblance2d.class, }; }

/*
    public boolean initialize(StsModel model)
    {
        if (!super.initialize(model)) return false;
        try
        {
            StsWindowFamily windowFamily = currentModel.getWindowFamily(currentModel.win3d);
            if (windowFamily != null)
            {
                windowFamily.adjustCursor(YDIR, currentY);
                windowFamily.adjustCursor(XDIR, currentX);
                return true;
            }
        }
        catch (Exception e)
        {
            e.printStackTrace();
            return false;
        }
        //int[] rowcol = lineSet.getRowColFromCoors(currentX, currentY);
        //lineSet.jumpToRowCol(rowcol, currentModel.win3d);
        return true;
    }
*/
}
