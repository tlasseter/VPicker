package com.Sts.Types.PreStack;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;

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
 * All prestack lines, 2d and 3d, defined in both processing (StsPreStackSegyLine2d & StsPreStackSegyLine3d),
 * and in loading (StsPreStackLine2d & StsPreStackLine3d) are subclassed from this abstract class.
 */
abstract public class StsPreStackLine extends StsSeismicLine implements Comparable
{
    /** volume or set of 2d line set this line belongs to */
    public StsPreStackLineSet lineSet = null;
    /** If true, this line is an inline; otherwise its a crossline. Always true for 2d */
    public boolean isInline = true;
    /** maximum number of offsets in line */
    public int nOffsetsMax = 0;
    /** number of traces in line */
    public int nLineTraces;
    /** Index number of last trace in each gather */
    public int[] nLastTraceInGathers;
    /** Spacing between CDP points */
    //	public float cdpSpacing = 0.0f;
    /** Sample rate in seconds */
    public float sRate = 0.0f;
    /** min value of trace offset for entire line */
    public float traceOffsetMin = StsParameters.largeFloat;
    /** max value of trace offset for entire line */
    public float traceOffsetMax = -StsParameters.largeFloat;
    /** indicates data has been NMOed*/
	public boolean isNMOed = false;
    /** prestackSeismicVolume min gather index on line; if line is an inline, this is the crossline index of first gather */
    public int minGatherIndex;
    /** prestackSeismicVolume max gather index on line; if line is an inline, this is the crossline index of last gather */
    public int maxGatherIndex;
    /** current line number; also total number of lines when volume finished; used in forming 3d file name: group.format.stemname.line.nLine */
	public int nLine = 0;
    /** filename for output gather data */
    public String gatherFilename;
    /** sorted (min to max) offsets for traces in gather */
    transient public double[] traceOffsets;
	transient public StsPreStackLineFileBlocks fileBlocks = null;
//    transient StsRandomAccessFile gatherFile = null;
    transient public int nGatherTraces;
    
    public float cdpAvgInterval = 0.0f;  
	
	static public final float nullValue = StsParameters.nullValue;
    static public final byte nullByte = StsParameters.nullByte;
    static final double roundOff = StsParameters.roundOff;

    static final public boolean debug = false;

    static public final String gatherFormat = "gathers"; /** array of x values for a 2d line in rotated local coordinate system */

    public StsPreStackLine()
    {
    }

    public StsPreStackLine(boolean persistent)
    {
        super(persistent);
    }

	public StsPreStackLine(StsFile file, String stsDirectory)
	{
		super(false);
		segyDirectory = file.getDirectory();
		segyFilename = file.getFilename();
		File sfile = new File(segyDirectory + segyFilename);
		segyLastModified = sfile.lastModified();
		this.stsDirectory = stsDirectory;
	}

    protected StsPreStackLine(String directory, String filename, StsModel model, StsPreStackLineSet lineSet) throws FileNotFoundException, StsException
	{
		super(false);
        String pathname = directory + filename;
		if (! (new File(pathname)).exists())
		{
			throw new FileNotFoundException();
		}
		if(!StsParameterFile.initialReadObjectFields(pathname, this, StsPreStackLine.class, StsBoundingBox.class))
        {
//            throw new IOException();
        }
        setName(getStemname());
		stsDirectory = directory;
//		traceOffsets = getAttributeArray(stsDirectory, attributeFilename, "OFFSET");
		initializeScale();
        // a null lineSet indicates this line is temporary: is not part of a set and should not be added to db
        if(lineSet == null) return;
        this.lineSet = lineSet;
		addToModel();
		if (!initialize(model))
		{
			throw new FileNotFoundException(pathname);
		}
		isVisible = true;
	}

    protected StsPreStackLine(String directory, String filename, int lineIndex, StsModel model, StsPreStackLineSet lineSet) throws FileNotFoundException, StsException
    {
        this(directory, filename, model, lineSet);
        this.lineIndex = lineIndex;
    }

    protected StsPreStackLine(StsFile file, int lineIndex, StsModel model, StsPreStackLineSet lineSet) throws FileNotFoundException, StsException
    {
        this(file, model, lineSet);
        this.lineIndex = lineIndex;
    }



    public boolean checkFiles()
    {
        File file;
        file = new File(lineSet.stsDirectory + gatherFilename);
        if(!file.exists())
        {
            new StsMessage(currentModel.win3d, StsMessage.INFO, "Failed to find gather file " + stsDirectory + gatherFilename);
            return false;
        }
        file = new File(stsDirectory + attributeFilename);
        if(!file.exists())
        {
            new StsMessage(currentModel.win3d, StsMessage.INFO, "Failed to find attribute file " + stsDirectory + attributeFilename);
            return false;
        }
        return true;
    }

 /*
    protected StsPreStackLine(StsFile file, int setIndex, StsModel model, StsPreStackLineSet lineSet) throws FileNotFoundException, StsException
    {
        this(file.getDirectory(), file.getFilename(), model, lineSet);
    }
 */
    protected StsPreStackLine(StsFile file, StsModel model, StsPreStackLineSet lineSet) throws FileNotFoundException, StsException
    {
        this(file.getDirectory(), file.getFilename(), model, lineSet);
    }

    public boolean initialize(StsFile file, StsModel model)
    {
        try
        {
            String pathname = file.getDirectory() + file.getFilename();
            StsParameterFile.initialReadObjectFields(pathname, this, StsPreStackLine.class, StsBoundingBox.class);
            setName(getStemname());
            stsDirectory = file.getDirectory();
            if (!initialize(model))
            {
                return false;
            }
            isVisible = true;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsPreStackLine3d.loadFile() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean initialize(StsModel model)
    {
        try
        {
//            traceOffsets = getAttributeArray(stsDirectory, attributeFilename, "OFFSET");
            StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
			fileBlocks = new StsPreStackLineFileBlocks(nCols, nOffsetsMax, nLastTraceInGathers, nSlices, lineSet, gatherFilename, 4, memoryManager);
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError(this, "initialize(model)");
            return false;
        }
    }

    //TODO: convert offset array to a MappedDoubleBuffer and have it managed by StsBlocksMemoryManager
	public double[] getTraceOffsets()
	{
		try
		{
			if (traceOffsets == null)
				traceOffsets = getAttributeArray("OFFSET");
		}
		catch(FileNotFoundException e)
		{
			outputWarningException("getTraceOffsets()", " Couldn't find file: " + attributeFilename, e);
			return null;
		}
		return traceOffsets;
    }

    public void clearFileBlocks()
    {
        StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
        memoryManager.clearBlocks(fileBlocks);
    }

    public boolean initialize(StsModel model, String mode) throws FileNotFoundException
    {
        return true;
    }

	protected void initializeScale()
	{
		signedIntDataZero = Math.round(254 * dataMin / (dataMin - dataMax) - 127);
	}
/*
	public StsRandomAccessFile getGatherFile()
	{
		try
		{
			return new StsRandomAccessFile(stsDirectory + gatherFilename, "r");
		}
		catch (Exception e)
		{
			new StsMessage(currentModel.win3d, StsMessage.INFO, "Failed to open file " + stsDirectory + gatherFilename+"\n"+e);
			return null;
		}
	}
*/
    public int getMaxTracesPerGather()
    {
        int max = 0;
        for(int i=minGatherIndex; i<maxGatherIndex; i++)
        {
            if(max < getNTracesInGather(i))
                max = getNTracesInGather(i);
        }
        return max;
    }
/*
	protected boolean lineContainsGather(int lineIndex, int gatherIndex)
	{
		if(this.lineIndex != lineIndex)return false;
		if(this.minGatherIndex > gatherIndex)return false;
		if(this.maxGatherIndex < gatherIndex)return false;
		return true;
	}
*/

    public double[] getAttributeArray(String attributeName)  throws FileNotFoundException
    {
        return super.getAttributeArray(attributeName, nLineTraces);
    }

    public double[] getAttributeArray(String attributeName, int row, int col) throws FileNotFoundException
    {
        int start = 0;
        int idx = 0, num = 0;

        double[] attArray = getAttributeArray(attributeName, nLineTraces);
        if(attArray == null)
            return null;

        if(isInline)
            idx = col;
        else
            idx = row;

        start = getNFirstGatherTrace(idx);
        num = getNTracesInGather(idx);
        if(num == 0)
            return null;

        double[] newArray = new double[num];
        System.arraycopy(attArray, start, newArray, 0, num);
//        for(int i=0; i<num; i++)
//            newArray[i] = attArray[start+i];

        return newArray;
	}

    public double[][] getGatherAttributes(int row, int col)
    {
        double[][] allAttributes = new double[nAttributes][];
        try
        {
            for(int n = 0; n < nAttributes; n++)
                allAttributes[n] = getAttributeArray(attributeNames[n], nLineTraces);
            return allAttributes;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getGatherAttributes", "row: " + row + " col: " + col, e);
            return null;
        }
    }

    public int getNTraces()
	{
		return nLineTraces;
	}

    public int getNTracesInGather(int row, int col)
	{
        if(isInline)
            return getNTracesInGather(col);
        else
            return getNTracesInGather(row);
    }

    public int getNTracesInGather(int index)
	{
		int nFirstGatherTrace, nLastGatherTrace;

		int nGather = index - minGatherIndex;
		if (nGather < 0 || nGather >= nCols) return 0;
		if (nGather == 0)
			nFirstGatherTrace = 0;
		else
			nFirstGatherTrace = nLastTraceInGathers[nGather - 1] + 1;
		nLastGatherTrace = nLastTraceInGathers[nGather];
		return nLastGatherTrace - nFirstGatherTrace + 1;
	}

	public long getGatherFileOffset(int index)
	{
		int nFirstGatherTrace = getNFirstGatherTrace(index);
		if(nFirstGatherTrace == -1) return -1;
		return nSlices * nFirstGatherTrace;
	}

	public long getGatherFileOffset(int row, int col)
	{
		if(isInline)
			return getGatherFileOffset(col);
		else
			return getGatherFileOffset(row);
	}

	public int getGatherFileLength(int index)
	{
		int nFirstGatherTrace = getNFirstGatherTrace(index);
		if(nFirstGatherTrace == -1) return -1;
		int nLastGatherTrace = getNLastGatherTrace(index);
		if(nLastGatherTrace == -1) return -1;
		int length = nSlices * (nLastGatherTrace - nFirstGatherTrace + 1);
		return length;
	}

	public int getGatherFileLength(int row, int col)
	{
		if(isInline)
			return getGatherFileLength(col);
		else
			return getGatherFileLength(row);
	}

	public int getNFirstGatherTrace(int index)
	{
		int nGather = index - minGatherIndex;
		if (nGather < 0 || nGather >= nCols) return -1;
		if (nGather == 0)
			return 0;
		else
			return nLastTraceInGathers[nGather - 1] + 1;
	}

	public int getNLastGatherTrace(int index)
	{
		int nGather = index - minGatherIndex;
		if (nGather < 0 || nGather >= nCols) return -1;
		return nLastTraceInGathers[nGather];
	}
    public void setSegyFilename(String segyFilename)
    {
        this.segyFilename = segyFilename;
    }

    public String getSegyFilename()
    {
        return segyFilename;
    }

    public void setName(String name)
    {
        this.stemname = name;
    }

    public String getCreateHeaderFilename()
    {
        if(headerFilename == null)
            headerFilename = createFilename(headerFormat);
        return headerFilename;
    }

    public String getCreateGatherFilename()
    {
        if(gatherFilename == null)
            gatherFilename = createFilename(gatherFormat);
        return gatherFilename;
    }

    public String getCreateAttributeFilename()
    {
        if(attributeFilename == null)
            attributeFilename = createFilename(attributeFormat);
        return attributeFilename;
    }

	public boolean overwriteFiles()
	{
		File fi = new File(stsDirectory + getCreateGatherFilename());
		return fi.exists();
	}
/*
    public boolean writeHeaderFile()
    {
        try
        {
            StsParameterFile.writeObjectFields(stsDirectory + createHeaderFilename(), this, StsPreStackLine.class, StsBoundingBox.class);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException(getClass().getName() + ".writeHeaderFile() failed.", e, StsException.WARNING);
            return false;
        }
    }
*/
    public int compareTo(Object other)
	{
		StsPreStackLine otherLine = (StsPreStackLine)other;

		if(lineIndex < otherLine.lineIndex) return -1;
		if(lineIndex > otherLine.lineIndex) return 1;
		if(maxGatherIndex < otherLine.minGatherIndex) return -1;
		if(minGatherIndex > otherLine.maxGatherIndex) return 1;
		StsException.systemError("StsPreStackLine.compareTo() failed.  Lines overlap for lineIndex " + lineIndex + " firstLineGatherIndices " + minGatherIndex + "-" + maxGatherIndex +
			 " secondLineGatherIndices " + otherLine.minGatherIndex + "-" + otherLine.maxGatherIndex);
	    return 0;
	}

	public int compareTo(int gatherIndex)
	{
		if(this.maxGatherIndex < gatherIndex) return -1;
		if(this.minGatherIndex > gatherIndex) return 1;
		return 0;
	}
    public int compareTo(int lineIndex, int gatherIndex)
	{
		if(this.lineIndex < lineIndex) return -1;
		if(this.lineIndex > lineIndex) return 1;
		if(this.maxGatherIndex < gatherIndex) return -1;
		if(this.minGatherIndex > gatherIndex) return 1;
	    return 0;
	}

    public void drawTextureTileDepthSurface(StsPreStackVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
    }
public boolean setLineXYsFromCDPs()
    {
        try
		{
			double[] lineCdpX = getAttributeArray(StsSEGYFormat.CDP_X);
			double[] lineCdpY = getAttributeArray(StsSEGYFormat.CDP_Y);
			double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP);
			cdpX = new float[nCols];
			cdpY = new float[nCols];
			cdp = new int[nCols];

            // Check to see if CDPs are correct
            if((lineCdpX[0] == 0.0) && (lineCdpY[0] == 0.0))
                return false;

			for(int n = 0; n < nCols; n++)
			{
				int t = nLastTraceInGathers[n];
                float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(lineCdpX[t], lineCdpY[t]);
				cdpX[n] = xy[0];
				cdpY[n] = xy[1];
				cdp[n] = (int)lineCdp[t];
			}
            computeNormals();
            return true;
		}
		catch(FileNotFoundException fnfe)
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
			return false;
		}
	}public double[][] getLineXYCDPs()
    {
        try
		{
			double[] lineCdpX = getAttributeArray(StsSEGYFormat.CDP_X);
			double[] lineCdpY = getAttributeArray(StsSEGYFormat.CDP_Y);
			double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP);
            double[][] xyc = new double[nCols][3];

            for(int n = 0; n < nCols; n++)
			{
				int t = nLastTraceInGathers[n];
                xyc[n][0] = lineCdpX[t];
                xyc[n][1] = lineCdpY[t];
                xyc[n][2] = lineCdp[t];
			}
            return xyc;
		}
		catch(FileNotFoundException fnfe)
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
			return null;
		}
	}
	
	public void setCdpAvgInterval(float cdpAvgInterval) {this.cdpAvgInterval = cdpAvgInterval;}

    public float getCdpAvgInterval() {return cdpAvgInterval;}

    /**
     * returns CDP number from header data
     * 
     * @param column  index of gather in question (sequential counter from beginning of line)
     * @return cdp number or -1 if line is null or no CDP information is loaded
     */
    public int getCDP(int column)
    {
        int cdpNum = -1;
        if (cdp == null)
        {
            try
            {
                double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP);
                cdp = new int[nCols];
                for (int n = 0; n < nCols; n++)
                    cdp[n] = (int) lineCdp[nLastTraceInGathers[n]];
            } catch (Exception e)
            {
                e.printStackTrace();
            }
        }
        if (cdp != null && column < cdp.length) cdpNum = cdp[column];
        return cdpNum;
    }
    
}
