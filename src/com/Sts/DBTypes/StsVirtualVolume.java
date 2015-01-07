package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.text.*;
import java.util.*;

public abstract class StsVirtualVolume extends StsSeismicVolume implements StsTreeObjectI
{
    protected boolean dataRangeChanged = false;
    public StsObjectRefList volumes = null;

    transient public static final int UNDEFINED = 0;

    transient public static final int SEISMIC_MATH = 0;
    transient public static final int SEISMIC_XPLOT_MATH = 1;
    transient public static final int SEISMIC_BLEND = 2;
    transient public static final int SEISMIC_FILTER = 3;
    transient public static final int EP_SEISMIC_FILTER = 4;
    transient public static final int RGB_BLEND = 5;
    transient public static final int SENSOR_VOLUME = 6;
    transient public static final String[] TYPE_STRINGS = { "Math", "Blended", "Crossplot", "Filters", "EPFilter", "Color", "Sensor"};
    transient public int[] nRowCols = new int[2];

    public StsVirtualVolume()
    {
        System.out.println("VirtualVolume constructor called.");
    }
    public StsVirtualVolume(boolean persistent)
    {
        super(persistent);
	}
    public void setVolumes(StsObject[] volumeList)
    {
        volumes = StsObjectRefList.constructor(2, 2, "volumes", this);
        volumes.add(volumeList);
    }

    public void setVolume(StsObject volume)
    {
        volumes = StsObjectRefList.constructor(2, 2, "volumes", this);
        volumes.add(volume);
    }

    // if crossline+ direction is 90 degrees CCW from inline+, this is isColCCW; otherwise not
    // angle is from X+ direction to inline+ direction (0 to 360 degrees)
    public boolean initialize(StsModel model)
    {
        try
        {
            if(!super.initialize(model))
			{
				delete();
				return false;
			}
			rowFloatFilename = null;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsVirtualVolume.classInitialize() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }

    /** if a required virtualVolume is deleted (index will now be -1), delete this virtual volume from class. */
    public boolean isDisplayable()
    {
        Iterator iter = volumes.getIterator();
        while(iter.hasNext())
        {
            StsObject stsObject = (StsObject)iter.next();
            if(stsObject == null) return false;
            if(!stsObject.isPersistent()) return false;
        }
        return true;
    }

    /** copy routine copies StsSeismicVolume to StsVirtualVolume, but these members should be nulled. */
	public void clearNonRelevantMembers()
	{
		rowFloatFilename = null;
	}

    public boolean getReadoutEnabled()
    {
        return readoutEnabled;
    }

    public void setReadoutEnabled(boolean enabled)
    {
        readoutEnabled = enabled;
        return;
    }

    static public StsVirtualVolumeClass getVirtualVolumeClass()
    {
        return (StsVirtualVolumeClass) currentModel.getCreateStsClass(StsVirtualVolume.class);
    }

    static public StsVirtualVolumeClass getVirtualVolumeClass(StsObject object)
    {
        return (StsVirtualVolumeClass) currentModel.getCreateStsClass(object);
    }

    public boolean getIsPixelMode()
    {
        return getVirtualVolumeClass(this).getIsPixelMode();
    }

    // compute a rough histogram with a few planes
    public void computeHistogram()
    {
    	clearHistogram();
    	for(int i=0; i<nRows; i=i+10)
    	{
    		byte[] plane = readRowPlaneByteData(i);
    		if(plane == null) continue;
    		for(int j=0; j<plane.length; j++)
        		accumulateHistogram(plane[j]);
    	}
    	calculateHistogram();
    	setDataHistogram();
    }
    public boolean allocateVolumes(String mode)
    {
        return true;
    }
    public boolean allocateVolumes()
    {
        return true;
    }

/*
	public ByteBuffer readByteBufferPlane(int dirNo, float dirCoordinate)
	{
		int nVolumes = volumes.getSize();

	   // check that all volumes are available; if not, delete this volume
		for (int n = 0; n < nVolumes; n++)
		{
			StsObject volume = (StsObject)volumes.getElement(n);
			StsClass volumeClass = currentModel.getStsClass(volume.getClass());
			if(!volumeClass.contains(volume))
			{
				getStsClass().delete(this);
				return null;
			}
        }

		try
		{
			int nPlane = this.getCursorPlaneIndex(dirNo, dirCoordinate);
			if(nPlane == -1) return null;
			MappedByteBuffer[] planeBuffers = new MappedByteBuffer[nVolumes];
			for(int n = 0; n < nVolumes; n++)
			{
				StsSeismicVolume volume = (StsSeismicVolume)volumes.getElement(n);
				planeBuffers[n] = volume.filesMapBlocks[dirNo].getByteBufferPlane(nPlane, FileChannel.MapMode.READ_ONLY);
				planeBuffers[n].slice();
			}
			return processBytePlaneData(planeBuffers);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolume.readPlaneData() failed.", e, StsException.WARNING);
			return null;
		}
	}

	ByteBuffer processBytePlaneData(MappedByteBuffer[] planeBuffers)
	{
		return null;
	}
*/

   public ByteBuffer readByteBufferPlane(int dir, float dirCoordinate)
   {
	   byte[] byteData = readBytePlaneData(dir, dirCoordinate);
       if(byteData == null) return null;
       return ByteBuffer.wrap(byteData);
   }

    public StsVirtualVolumeClass getSeismicVolumeClass()
    {
        return (StsVirtualVolumeClass) currentModel.getCreateStsClass(getClass());
    }

    public byte[] readRowPlaneByteData(int nPlane)
    {
        return readBytePlaneData(YDIR, getYCoor(nPlane));
    }

    public byte[] readBytePlaneData(int dir, float dirCoordinate)
    {
        if(volumes == null) return null;
        int nVolumes = volumes.getSize();

       // check that all volumes are available; if not, delete this volume
        for (int n = 0; n < nVolumes; n++)
        {
            StsObject volume = (StsObject)volumes.getElement(n);
            StsClass volumeClass = currentModel.getStsClass(volume.getClass());
            if(!volumeClass.contains(volume))
            {
                getStsClass().delete(this);
                return null;
            }
        }

        byte[] planeData;

        if (volumes == null) return null;

        byte[][] planes = new byte[nVolumes][];

        try
        {
            // Calculate the requested plane
            for (int n = 0; n < nVolumes; n++)
            {
                StsVolumeDisplayable volume = (StsVolumeDisplayable) volumes.getElement(n);
                planes[n] = volume.readBytePlaneData(dir, dirCoordinate);
                if (planes[n] == null) return null;
            }
            nRowCols = getCursorDisplayNRowCols(dir);
            planeData = processPlaneData(planes);
//            dataCubes[dirNo].setPlaneData(nPlane, planeData);
            return planeData;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }

    public abstract byte[] processPlaneData(byte[][] planeValues);

    public boolean isAValueNull(byte value0, byte value1)
    {
        return value0 == -1 || value1 == -1;
    }

    public void initializeData()
    {
        dataMax = -Float.MAX_VALUE;
        dataMin = Float.MAX_VALUE;
        //clearHistogram();
        for(int i=getRowMin(); i<getRowMax(); i++)
        {
            readBytePlaneData(StsCursor3d.YDIR, getYCoor(i));
        }
        //calculateHistogram();
        //setDataHistogram();
    }

    public byte[] computeUnscaledByteValues(float[] values)
    {
        if(values == null) return null;
        int nValues = values.length;
        //boolean needHistogram = true;

        float oldDataMin = dataMin;
        float oldDataMax = dataMax;

        int nStart = 0;
        if(dataMin == Float.MAX_VALUE)
        {
            //clearHistogram();
            //needHistogram = true;
            for (nStart = 0; nStart < nValues; nStart++)
            {
                if (values[nStart] != StsParameters.nullValue)
                {
                    dataMin = values[nStart];
                    dataMax = values[nStart];
                    break;
                }
            }
        }
        else
            nStart = 0;

        for (int n = nStart + 1; n < nValues; n++)
        {
            if (values[n] == StsParameters.nullValue) continue;
            if (values[n] < dataMin)
                 dataMin = values[n];
             else if (values[n] > dataMax)
                dataMax = values[n];
        }
        if(dataMin != oldDataMin || dataMax != oldDataMax)
            resetDataRange(dataMin, dataMax, false);

        float scale = 254 / (dataMax - dataMin);

        byte[] planeData = new byte[nValues];
        for (int n = 0; n < nValues; n++)
        {
            if(values[n] == StsParameters.nullValue)
                planeData[n] = -1;
            else
            {
                float scaledFloat = (values[n] - dataMin) * scale;
                int scaledInt = Math.round(scaledFloat);
                planeData[n] = StsMath.unsignedIntToUnsignedByte(scaledInt);
                //if(needHistogram)
                    //accumulateHistogram(planeData[n]);
            }
        }
        /*
        if(needHistogram)
        {
            calculateHistogram();
            setDataHistogram();
        }
        */
        return planeData;
    }

    public void resetDataRange(float min, float max, boolean maintain)
    {
        float[] editRange;
        if(maintain)
        {
            editRange = colorscale.getEditRange();
            colorscale.setRange(min, max);
            colorscale.setEditRange(editRange[0], editRange[1]);
        }
        else
        {
            colorscale.setRange(min, max);
        }
        dataMin = min;
        dataMax = max;
        dataRangeChanged = false;
    }

    public String getSegyFileDate()
    {
        if (segyLastModified == 0)
        {
            File segyFile = new File(segyDirectory + segyFilename);
            if (segyFile != null)
            {
                segyLastModified = segyFile.lastModified();
            }
        }
        DateFormat dateFormat = DateFormat.getDateTimeInstance();
        return dateFormat.format(new Date(segyLastModified));
    }

    public String getDate()
    {
        return null;
    }

    public String getLabel()
    {
        return stemname;
    }

    /** For the 3 cursor directions, get the coordinates on each cursor plane
     *  translated to equivalent grid coordinates.  The XDIR plane is in directions y+,z+
     *  The YDIR plane is in directions x+,z+ and the ZDIR plane is in directions x+,y+
     */
    public float[][] getGridCoordinateRanges(float[][] axisRanges, int currentDirNo, boolean axesFlipped)
    {
        float[][] coordinateRanges = new float[2][2];

        if (currentDirNo == StsCursor3d.XDIR)
        {
            coordinateRanges[0][0] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[0][0]);
            coordinateRanges[0][1] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[0][1]);
            coordinateRanges[1][0] = axisRanges[1][0]; // zMax
            coordinateRanges[1][1] = axisRanges[1][1]; // zMin
        }
        else if (currentDirNo == StsCursor3d.YDIR)
        {
            coordinateRanges[0][0] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[0][0]);
            coordinateRanges[0][1] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[0][1]);
            coordinateRanges[1][0] = axisRanges[1][0]; // zMax
            coordinateRanges[1][1] = axisRanges[1][1]; // zMin
        }
        else if (currentDirNo == StsCursor3d.ZDIR)
        {
            if (!axesFlipped)
            {
                coordinateRanges[1][0] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[1][0]);
                coordinateRanges[1][1] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[1][1]);
                coordinateRanges[0][0] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[0][0]);
                coordinateRanges[0][1] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[0][1]);
            }
            else
            {
                coordinateRanges[1][0] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[1][0]);
                coordinateRanges[1][1] = getNumFromCoor(StsCursor3d.XDIR, axisRanges[1][1]);
                coordinateRanges[0][0] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[0][0]);
                coordinateRanges[0][1] = getNumFromCoor(StsCursor3d.YDIR, axisRanges[0][1]);
            }
        }
        return coordinateRanges;
    }

    public int getIntValue(int row, int col, int slice)
    {
        if (row < 0 || row >= nRows || col < 0 || col >= nCols || slice < 0 || slice >= nSlices)
        {
            StsException.systemError("StsSeismicVolume.getValue() failed for row: " +
                                     row + " col: " + col + " plane: " + slice);
            return 0;
        }
        int n = row * nCols + col;
        byte[] planeData = readBytePlaneData(StsCursor3d.ZDIR, getZCoor(slice));
        return (int) planeData[n];
    }

    public int getIntValue(float[] xyz)
    {
        int row = (int) getRowCoor(xyz[1]);
        int col = (int) getColCoor(xyz[0]);
        int plane = (int) getSliceCoor(xyz[2]);
        return getIntValue(row, col, plane);
    }

    public String getVolumeOneName()
    {
        return getVolumeName(volumes.getElement(0));
    }

    private String getVolumeName(Object volume)
    {
        return ( (StsObject) volume).getName();
    }

    public String getVolumeTwoName()
    {
        if (volumes.getSize() == 2)
        {
            return getVolumeName(volumes.getElement(1));
        }
        else
        {
            return "None";
        }
    }

    public String getVolumeThreeName()
    {
        if (volumes.getSize() == 3)
        {
            return getVolumeName(volumes.getElement(2));
        }
        else
        {
            return "None";
        }
    }

    public boolean getIsRegular()
    {return isRegular;
    }

    public void setIsRegular(boolean isRegular)
    {this.isRegular = isRegular;
    }

    public StsSpectrum getSpectrum()
    {return colorscale.getSpectrum();
    }
   
    // Currently this method is identical to superclass; need to understand if/when histogram is to be set.
    // Note also the currentColorscale business in superclass StsSeismicVolume which is now commented out because it's dysfunctional
    // or function not clear.  TJL 3/18/07

    public StsColorscale getColorscale()
    {
        setDataHistogram();
        return colorscale;
    }

    public void initializeColorscaleActionListener()
    {
        StsColorscale colorscale = getColorscale();
        colorscale.addActionListener(this);
    }
/*
    public void setColorscale(StsColorscale colorscale)
    {
        this.colorscale = colorscale;
        currentModel.win3dDisplayAll();
    }
*/
    /*
    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            seismicColorList.setColorListChanged(true);
            currentModel.clearDisplayTextured3dCursors(this);
            fireActionPerformed(e);
            currentModel.viewObjectRepaint(this, this);
        }
        else
        {
            fireActionPerformed(e);
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
    }

	public void fireActionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            currentModel.displayIfCursor3dObjectChanged(this);
        }
        else
        {
            fireActionPerformed(null);
        }
        return;
    }
    */
    /*
    public void initializeSpectrum()
    {
        try
        {
            StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
            StsSeismicVolumeClass virtualClass = (StsSeismicVolumeClass)currentModel.getStsClass(StsSeismicVolume.class);
            //
            // Determine Combined DataMin and DataMax
            //
            dataMin = Float.MAX_VALUE;
            dataMax = -Float.MAX_VALUE;
            colorscale = new StsColorscale("Virtual PostStack3d", spectrumClass.getSpectrum(virtualClass.getSeismicSpectrumName()), dataMin, dataMax);
            colorscale.setEditRange(dataMin, dataMax);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.initializeSpectrum() failed.",
                                         e, StsException.WARNING);
        }
    }

    public void setSpectrumDialog(StsSpectrumDialog spectrumDialog)
    {
        this.spectrumDialog = spectrumDialog;
    }
     */
    public boolean getIsVisibleOnCursor()
    {
        return isVisible && getVirtualVolumeClass(this).getIsVisibleOnCursor();
    }

    public StsColor getStsColor(float[] xyz)
    {
        int v = getIntValue(xyz);
        return colorscale.getSpectrum().getColor(v);
    }

    public StsColor getStsColor(int row, int col, float z)
    {
        int v = getIntValue(row, col, z);
        return colorscale.getSpectrum().getColor(v);
    }

    public StsColor getStsColor(float rowF, float colF, int plane)
    {
        int row = StsMath.roundOffInteger(rowF);
        int col = StsMath.roundOffInteger(colF);
        int v = getIntValue(row, col, plane);
        return colorscale.getSpectrum().getColor(getIntValue(row, col, plane));
    }

    public boolean getIsVisible()
    {return isVisible;
    }

    public byte[] getCurrentCursorPlaneData()
    {
        return null;
    }

    public abstract StsFieldBean[] getDisplayFields();
    public abstract StsFieldBean[] getPropertyFields();

    public Object[] getChildren()
    {return new Object[0];
    }

    abstract public StsObjectPanel getObjectPanel();

    public boolean anyDependencies()
    {
        return false;
    }

    public int getColorDisplayListNum(GL gl, boolean nullsFilled)
    {
        return 0;
    }


	void drawIsoSurfaces(GL gl)
	{
		drawByteIsoSurfaces(gl);
	}

    public void clearCache()
    {
        if(volumes == null) return;
        for(int i=0; i<volumes.getSize(); i++)
        {
            if(volumes.getElement(i) instanceof StsSeismicVolume)
                ((StsSeismicVolume) volumes.getElement(i)).clearCache();
            else if(volumes.getElement(i) instanceof StsCrossplot)
            {
                StsCrossplot cp = (StsCrossplot)volumes.getElement(i);
                for(int j=0; j<cp.getVolumes().length; j++)
                {
                    cp.getVolumes()[j].clearCache();
                }
            }
        }
//        dataCubeMemory.clearCache(dataCubes);
	}

    /** Virtual volumes don't use byteBlocks to store data. Return true indicating its ok. */
    public boolean setupRowByteBlocks() { return true; }
}
