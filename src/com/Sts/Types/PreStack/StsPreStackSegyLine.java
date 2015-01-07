package com.Sts.Types.PreStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import sun.io.*;

import javax.media.opengl.*;
import java.awt.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public abstract class StsPreStackSegyLine extends StsPreStackLine
{
    transient protected int nVolume = 0;
	transient public StsSEGYFormat.TraceHeader firstTraceHeader, secondTraceHeader, lastTraceHeader, firstLineLastTraceHeader, lastLineFirstTraceHeader;
	/** number of traces per gather is the same for all gathers */
	transient public boolean isTracesPerGatherSame = false;
	/** number of CDPs per line is the same for all lines */
	transient public boolean isCDPsPerLineSame = true;
	/** Average number of traces per gather if not the same for each gather. */
	transient public int nTracesPerGather = 1;
	/** Average number of CDPs per line if if not the same for each line. */
	transient public int nCDPsPerLine = -1;
	/** current trace in gather; also total number of traces in gather when finished */
	transient public int nGatherTrace = 0;
	transient protected FileChannel randomAccessSegyChannel;
	transient protected MappedByteBuffer inputBuffer = null;
	transient protected StsMappedFloatBuffer gatherFloatBuffer = null;

    transient protected int nInputBytesPerBlock = 0;
	transient protected int nBlockTraces = 0;
	transient protected int nBlockSamples = 0;
	transient protected int outputRowMin = 0;
	transient protected int nBlock = 0;
	transient protected int gatherSampleCount = 0;
	transient protected boolean isLastBlock = false;
	/** number of blocks in memory waiting to be written out for xlines and slices */
	transient protected int nBlocksInMemory = 0;
	transient protected long inputPosition = 0;

	transient protected int nGatherTracesWritten = 0;
	transient protected long gatherFileStartPosition = 0; // start position in file of buffer for next write
	transient protected long gatherFileEndPosition = 0; // end position of buffer in file after data has been written

	/** attribute records for attributes written */
	transient public StsSEGYFormatRec[] attributeRecords;
	/** indexes for key attribute records */
	transient int ilineAttributeIndex = -1;
	transient int xlineAttributeIndex = -1;
	transient int cdpxAttributeIndex = -1;
	transient int cdpyAttributeIndex = -1;
    transient int shtxAttributeIndex = -1;
	transient int shtyAttributeIndex = -1;
    transient int recxAttributeIndex = -1;
	transient int recyAttributeIndex = -1;
	transient int offsetAttributeIndex = -1;

	transient StsSeismicBoundingBox overrideGeometryBox = null;
	transient int inlineDecimation = 1;
	transient int xlineDecimation = 1;

	transient float originRowNum = -1;
	transient float originColNum = -1;

	transient boolean isCropped = false;

	/** data scaling: 254/(dataMax - dataMin) */
	transient protected float scale = 1.0f;

	transient boolean scanHistogram = false;
	transient double scanPercentage = 0.1;

	transient static boolean isMinMemory = true;

	/** attributes for each trace; sequential for each attribute [nAttributes][nTraces] */
	transient protected double[][] traceAttributes;
	/** length of trace attributes array */
	transient protected int traceAttributesLength = 10000;
	/** number of traces written to attribute array */
	transient int nAttributeTrace = 0;

	transient byte[] blockGatherOrderBytes = null, blockOffsetOrderBytes = null;
	transient double[] ampsUncorrectedGather;
	transient double[] ampsCorrectedGather;
	transient double[] ampsSemblance;
	transient protected StsMemAllocPreStackProcess memoryAllocation = null;
	transient protected byte[] segyTraceBytes;
	transient protected byte[] traceBytes;
	transient protected byte[] paddedTraceBytes;

	transient boolean canceled = false;

	/** timer used in testing performance */
	static protected StsTimer timer = null;
	static protected StsTimer totalTimer = null;
	/** set to true if timing tests are desired */
	static protected boolean runTimer = false;

	/** wizard frame used for locating isVisible dialogs */
	transient Frame frame;

	static String[] segySuffixes = new String[] {".sgy", ".segy"};

//	static final boolean debug = false;
	static public final boolean debugLine = false;
	static public final boolean debugGather = false;


	/** computes volume geometry from trace CDP ilines, crossline, X, and Y values */
    public abstract boolean analyzeGeometry();
    protected abstract boolean areTracesInSameGather(StsSEGYFormat.TraceHeader firstTraceHeader, StsSEGYFormat.TraceHeader secondTraceHeader);
	public abstract boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry);
    public abstract boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry, boolean multiVolumeFile, int volumeNum);
    public abstract void renameGathersFile();

    public StsPreStackSegyLine()
	{
		super();
	}

	public StsPreStackSegyLine(boolean persistent)
	{
		super(persistent);
	}

	protected StsPreStackSegyLine(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws FileNotFoundException, StsException
	{
		super(file, stsDirectory);
		this.frame = frame;
		if(runTimer)timer = new StsTimer();
		try
		{
            segyData = StsSegyData.constructor(seismicWizard, file, this, segyFormat);
            if(segyData == null)
            {
                throw new StsException("StsPreStackSegyLine.constructor failed: couldn't construct StsSegyData.");
            }
//            segyData.randomAccessSegyFile = new RandomAccessFile(file.getPathname(), "r");

			if(segyFormat == null)
				this.setSegyFormat(StsSEGYFormat.constructor(currentModel, StsSEGYFormat.PRESTACK_RAW));
			else
				this.setSegyFormat(segyFormat);

//            segyData.readFileHeader();
//			analyzeBinaryHdr(false);
            setStemname();
   //         zDomain = segyFormat.getZDomainString();
        }
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume(file, frame) failed.", e, StsException.WARNING);
		}
	}

    public void setStemname()
    {
        if(segyData.isMultiVolumeFile())
            stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes) + nVolume;
        else
            stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes);
    }

    public void setStemname(boolean ignoreMultiVolume)
    {
       if(segyData.isMultiVolumeFile() && !ignoreMultiVolume)
            stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes) + nVolume;
        else
            stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes);
    }

    public void setNVolume(int volumeIdx, boolean ignoreMultiVolume)
    {
        nVolume = volumeIdx;
        setStemname(ignoreMultiVolume);
    }

	public void setOverrideGeometry(StsSeismicBoundingBox boundingBox)
	{
		this.overrideGeometryBox = boundingBox;
	}

	public void setIsNMOed(boolean value)
	{
		this.isNMOed = value;
	}

	public void setInlineDecimation(int value)
	{
		inlineDecimation = value;
	}

	public void setXlineDecimation(int value)
	{
		xlineDecimation = value;
	}

    public StsSegyData getSegyData() { return segyData; }

    public byte[] getTextHeader()
	{
		return segyData.getTextHeader();
	}

	public byte[] getBinaryHeader()
	{
		return segyData.getBinaryHeader();
	}

	public int analyzeTraces(double scanPercentage, StsProgressPanel progressPanel) // throws IOException
	{
		try
		{
			this.scanPercentage = scanPercentage;
			scanHistogram = true;
			initializeSegyIO();

			int result = segyIO.analyzeTraces(segyData.randomAccessSegyFile, (float)scanPercentage, progressPanel);
            dataMin = segyIO.dataMin;
            dataMax = segyIO.dataMax;
            return result;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackSegyLine.analyzeTraces() failed.", e, StsException.WARNING);
			return StsAnalyzeTraces.ANALYSIS_UNKNOWN_EXCEPTION;
		}
	}

	public int[] findGatherTraceRange(int nStartTrace, int nTracesPerGather)
	{
		try
		{
			if (nTracesPerGather < 50)
			{
				nTracesPerGather = 50;
			}
			// scan two times as many traces as avg gather forward and back to find start and end
			int startIdx = Math.max(nStartTrace - (2 * nTracesPerGather), 0);
			int endIdx = Math.min(nStartTrace + (2 * nTracesPerGather), getNTraces());
			if (endIdx < startIdx)
			{
				return null;
			}
			int nTracesRead = endIdx - startIdx + 1;
			int bytesToRead = getBytesPerTrace() * nTracesRead;
			byte[] tracesData = new byte[bytesToRead];

			long ioff = segyData.fileHeaderSize + (long)startIdx * getBytesPerTrace();
			segyData.randomAccessSegyFile.seek(ioff);
			segyData.randomAccessSegyFile.read(tracesData);

			int nFirstTrace = findFirstGatherTrace(tracesData, startIdx, nStartTrace);
			if (nFirstTrace == -1)
			{
			   return null;
			}
			int nLastTrace = findLastGatherTrace(tracesData, startIdx, nStartTrace, nTracesRead);
			if (nLastTrace == -1)
			{
			   return null;
			}
			return new int[] {nFirstTrace, nLastTrace};
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackSegyLine.getGatherData() failed.", e, StsException.WARNING);
			return null;
		}
	}

	private int findFirstGatherTrace(byte[] tracesData, int firstTrace, int searchPointTrace)
	{
		// If we're at the first trace then no need to look further
		if (searchPointTrace == 0)
		{
			return 0;
		}
		try
		{
			StsSEGYFormat.TraceHeader firstTraceHeader, traceHeader;
			int searchPointTraceOffset = getBytesPerTrace() * (searchPointTrace - firstTrace);
			firstTraceHeader = getSegyFormat().constructTraceHeader(tracesData, searchPointTraceOffset, searchPointTrace, getIsLittleEndian());

			int offset = searchPointTraceOffset;
			// look backwards from the search point trace
			for (int n = searchPointTrace; n > firstTrace; n--)
			{
				offset -= getBytesPerTrace();
				traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, n - 1, getIsLittleEndian());
				if (! areTracesInSameGather(firstTraceHeader, traceHeader))
				{
					return n;
				}
			}
			// if we get to the start then assume that that is the first trace in the gather
			return firstTrace;
		}
		catch(Exception e)
		{
			return -1;
		}
	}

	private int findLastGatherTrace(byte[] tracesData, int firstTrace, int searchPointTrace, int nTraces)
	{
		try
		{
			StsSEGYFormat.TraceHeader firstTraceHeader, traceHeader;
			int searchPointTraceOffset = getBytesPerTrace() * (searchPointTrace - firstTrace);
			firstTraceHeader = getSegyFormat().constructTraceHeader(tracesData, searchPointTraceOffset, searchPointTrace, getIsLittleEndian());

			int offset = searchPointTraceOffset;
			int lastTrace = firstTrace + nTraces - 1;
			// look forwards from the search point trace
			for (int n = searchPointTrace; n < lastTrace; n++)
			{
				offset += getBytesPerTrace();
				traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, n + 1, getIsLittleEndian());
				if (! areTracesInSameGather(firstTraceHeader, traceHeader))
				{
					return n;
				}
			}
			// if we get to the end then assume that that is the last trace in the gather
			return lastTrace;
		}
		catch(Exception e)
		{
			return -1;
		}
	}

	private int findLastGatherTrace(byte[] tracesData, int nStartTrace, int nTracesRead, int nFirstTrace,
									int nTracesPerGather)
	{
		StsSEGYFormat.TraceHeader firstTraceHeader, traceHeader;
		int offset = (nFirstTrace - nStartTrace) * getBytesPerTrace();
		firstTraceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, nFirstTrace, getIsLittleEndian());
		float firstIline = firstTraceHeader.iLine;
		float firstXline = firstTraceHeader.xLine;

		int nTraceInNextGather = nFirstTrace + nTracesPerGather;
		int nEndTrace = nStartTrace + nTracesRead - 1;
		if(nTraceInNextGather <= nEndTrace)
			offset = (nTraceInNextGather - nStartTrace) * getBytesPerTrace();
		else
		{
			nTraceInNextGather = nEndTrace;
			offset = (nTracesRead - 1) * getBytesPerTrace();
		}
		// starting from what we think is a trace in the next gather, search back for this gather or forward for next
		int nTrace = nTraceInNextGather;
		traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, nTrace, getIsLittleEndian());
		float iline = traceHeader.iLine;
		float xline = traceHeader.xLine;
		int inc;
		if(iline != firstIline || xline != firstXline)
		{
			// this trace is different from first, so we must be beyond last trace of gather: search backwards
			while(nTrace > nStartTrace)
			{
				nTrace--;
				offset -= getBytesPerTrace();
				traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, nTrace, getIsLittleEndian());
				iline = traceHeader.iLine;
				xline = traceHeader.xLine;
				if(iline == firstIline && xline == firstXline)
					return nTrace;
			}
		}
		else
		{
			// this trace is same iline-xline as first, so we are still in gather: search forward for last
			while(nTrace < nEndTrace)
			{
				nTrace++;
				offset += getBytesPerTrace();
				traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, nTrace, getIsLittleEndian());
				iline = traceHeader.iLine;
				xline = traceHeader.xLine;
				if(iline != firstIline || xline != firstXline)
					return nTrace - 1;
			}
		}
		return -1;
	}

	public double getBinaryHeaderValue(String recName)
	{
		return getSegyFormat().getBinaryHdrValue(recName, getBinaryHeader(), getIsLittleEndian());
	}

	public double getBinaryHeaderValue(StsSEGYFormatRec rec)
	{
		return getSegyFormat().getBinaryHdrValue(rec.getUserName(), getBinaryHeader(), getIsLittleEndian());
	}

	public byte[] getTraceHeaderBinary(int nTrace)
	{
		long offset = segyData.fileHeaderSize + (long)segyData.bytesPerTrace * (long)nTrace;
		byte[] hdrBuf = new byte[getTraceHeaderSize()];

		try
		{
			segyData.randomAccessSegyFile.seek(0);
			segyData.randomAccessSegyFile.seek(offset);
			segyData.randomAccessSegyFile.read(hdrBuf);
			return hdrBuf;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.getTraceHeaderBinary() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	public StsSEGYFormat.TraceHeader getTrace(int n)
	{
        return segyData.getTrace(n);
	}

	public boolean fileExists()
	{
		getCreateGatherFilename();
		File file = new File(stsDirectory + gatherFilename);
		return file.exists();
	}

	public String readTextHdr(String encoding)
	{
		try
		{
			if(encoding != null && encoding.equals("Cp500"))
			{
				ByteToCharCp500 converter = new ByteToCharCp500();
				return new String(converter.convertAll(getTextHeader()));
			}
			else
				return new String(getTextHeader(), 0, getSegyFormat().textHeaderSize, encoding);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.readTextHdr() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	/*
	 public boolean readWritePreStackLines(StsPreStackProcessPanel panel, StsSEGYFormatRec[] attributeRecords)
	 {
			return readWritePreStackLines(panel, attributeRecords, false);
	 }
	 */
	//public boolean readWritePreStackLines(StsPreStackProcessPanel panel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry)

	public boolean writeHeaderFile()
	{
		try
		{
            headerFilename = createFilename(headerFormat);
            StsParameterFile.writeObjectFields(stsDirectory + headerFilename, this, StsPreStackLine.class, StsBoundingBox.class);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException(StsToolkit.getSimpleClassname(this) + ".writeHeaderFile() failed.", e, StsException.WARNING);
			return false;
		}
	}

	protected void deletePartialFiles()
	{
		deleteExistingFile(stsDirectory, gatherFilename);
//		panel.setDescriptionAndLevel("\t" + gatherFilename + " deleted.\n");
		deleteExistingFile(stsDirectory, attributeFilename);
//		panel.setDescriptionAndLevel("\t" + attributeFilename + " deleted.\n");
		String headerFilename = createHeaderFilename();
		deleteExistingFile(stsDirectory, headerFilename);
//		panel.setDescriptionAndLevel("\t" + headerFilename + " deleted.\n");
	}

	private void attributeMissingMessage(String name)
	{
		new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + name + " attribute");
	}

	void initializeLine()
	{
		nLastTraceInGathers = new int[10000];
		nCols = 0;
		nOffsetsMax = 0;
		nLineTraces = 0;
		traceOffsetMin = StsParameters.largeFloat;
		traceOffsetMax = -StsParameters.largeFloat;
		if(this.overrideGeometryBox != null)
		{
			xOrigin = overrideGeometryBox.xOrigin;
			yOrigin = overrideGeometryBox.yOrigin;
			xInc = overrideGeometryBox.xInc;
			yInc = overrideGeometryBox.yInc;
			originRowNum = overrideGeometryBox.rowNumMin;
			rowNumInc = overrideGeometryBox.rowNumInc;
			originColNum = overrideGeometryBox.colNumMin;
			colNumInc = overrideGeometryBox.colNumInc;
			xMin = overrideGeometryBox.xMin;
			yMin = overrideGeometryBox.yMin;
			zMin = overrideGeometryBox.zMin;
			zInc = overrideGeometryBox.zInc;
			nSlices = overrideGeometryBox.nSlices;
			zMax = zMin + zInc*(nSlices-1);
/*
			if(inlineDecimation > 1)
			{
				xInc *= xlineDecimation;
				yInc *= inlineDecimation;
//				rowNumInc /= inlineDecimation;
//				colNumInc /= xlineDecimation;
			}
*/
		}
	}

	protected void setLastTraceInGather(int nTrace)
	{
		nLastTraceInGathers[nCols++] = nTrace;
	}

	protected void mapInputBuffer()
	{
		try
		{
//			if(debug)System.out.println("Remapping segy file channel. inputPosition: " + inputPosition + " nInputBytes: " + nInputBytesPerBlock);
			inputBuffer = null;
			inputBuffer = randomAccessSegyChannel.map(FileChannel.MapMode.READ_ONLY, inputPosition, nInputBytesPerBlock);
			inputPosition += nInputBytesPerBlock;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackSegyLine3d.mapInputBuffer() failed.", e, StsException.WARNING);
		}
	}

    public boolean openEmptyGatherFile()
	{
		try
		{
            gatherFilename = null;
            gatherFilename = createFilename(gatherFormat, "none");
			gatherFloatBuffer = StsMappedFloatBuffer.constructor(stsDirectory, gatherFilename, "rw");
			if(gatherFloatBuffer == null)return false;
            return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsPreStackSegyLine3d.openFiles() failed.");
			return false;
		}
	}

    public boolean openFiles()
	{
		try
		{
            gatherFilename = null;
            getCreateGatherFilename();
			deleteExistingFile(stsDirectory, gatherFilename);
			gatherFloatBuffer = StsMappedFloatBuffer.constructor(stsDirectory, gatherFilename, "rw");
			if(gatherFloatBuffer == null)return false;

            attributeFilename = null;
            getCreateAttributeFilename();
            deleteExistingFile(stsDirectory, attributeFilename);

            headerFilename = null;
            getCreateHeaderFilename();
            deleteExistingFile(stsDirectory, headerFilename);

            return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsPreStackSegyLine3d.openFiles() failed.");
			return false;
		}
	}

    protected boolean writeLineFiles()
    {
        closeGatherFile();
        trimGathersArray();
        writeAttributesFile();
        writeHeaderFile();
        return true;
    }

    protected void writeAttributesFile()
	{
		StsMappedDoubleBuffer attributeBuffer = null;

		try
		{
            attributeFilename = createFilename(attributeFormat);
			deleteExistingFile(stsDirectory, attributeFilename);
            attributeBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "rw");
			if(attributeBuffer == null)return;
			if(!attributeBuffer.map(0, nAttributes * nLineTraces)) return;
			for(int n = 0; n < nAttributes; n++)
				attributeBuffer.put(traceAttributes[n], 0, nLineTraces);
			attributeBuffer.close();
            nAttributeTrace = 0;
		}
		catch(Exception e)
		{
			StsException.systemError("StsPreStackSegyLine3d.writeAttributesFile() failed.");
		}
		finally
		{
			if(attributeBuffer != null) attributeBuffer.close();
		}
	}

	public void closeGatherFile()
	{
        if(gatherFloatBuffer != null)
        {
            gatherFloatBuffer.clear0();
            gatherFloatBuffer.clean();
            gatherFloatBuffer.close();
        }
	}

	protected void increaseTraceAttributesArray()
	{
		int length = traceAttributesLength;
		int newLength = length + 10000;
		double[][] newTraceAttributes = new double[nAttributes][newLength];
		for(int n = 0; n < nAttributes; n++)
			System.arraycopy(traceAttributes[n], 0, newTraceAttributes[n], 0, length);
		traceAttributes = newTraceAttributes;
		traceAttributesLength = newLength;
	}

	protected void clearInputBuffer()
	{
		if(runTimer)timer.start();

		if(inputBuffer != null)
		{
			inputBuffer.clear();
			StsToolkit.clean(inputBuffer);
			inputBuffer = null;
		}

		if(runTimer)timer.stopPrint("Time to clear input buffer for block " + nBlock + ":");
	}

    protected void clearOutputBuffer()
	{
		if(runTimer)timer.start();

		if(gatherFloatBuffer != null)
		{
			gatherFloatBuffer.clear();
			inputBuffer = null;
		}

		if(runTimer)timer.stopPrint("Time to clear input buffer for block " + nBlock + ":");
	}

    public StsSEGYFormat getSEGYFormat()
	{
		return getSegyFormat();
	}

	public void setSEGYFormat(StsSEGYFormat segyFormat)
	{
		this.setSegyFormat(segyFormat);
	}

    public int getSampleFormat()
	{
		return segyData.getSampleFormat();
	}

	public void setNSamples(int nSlices)
	{
		this.nSlices = nSlices;
	}

	public int getNSamples()
	{
		return nSlices;
	}

	public void setTMin(float tMin)
	{
		this.zMin = tMin;
	}

	public float getTMin()
	{
		return zMin;
	}

	public void setTMax(float tMax)
	{
		this.zMax = tMax;
	}

	public float getTMax()
	{
		return zMax;
	}

	public void setTInc(float tInc)
	{
		this.zInc = tInc;
	}

	public float getTInc()
	{
		return zInc;
	}

	public int getTotalTraces()
	{
		return getNTraces();
	}
	public int getNTracesPerGather()
	{
		return this.nTracesPerGather;
	}
    public int getNumVolumes()
    {
        return segyData.getNumVolumes();
    }
	public float getVerticalScalar()
	{
		return verticalScalar;
	}

	public float getHorizontalScalar()
	{
		return horizontalScalar;
	}

	public void setVerticalScalar(float vertScale)
	{
		verticalScalar = vertScale;
	}

	public void setHorizontalScalar(float horzScale)
	{
		horizontalScalar = horzScale;
	}

	public boolean isDone()
	{
		return false;
	}

	public  void drawTextureTileDepthSurface(StsPreStackVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
		System.out.println("drawTextureTileDepthSurfaceInTime called ");
	}

    protected void trimGathersArray()
    {
        int[] newArray = new int[nCols];
        System.arraycopy(nLastTraceInGathers, 0, newArray, 0, nCols);
        nLastTraceInGathers = newArray;
    }

    public void cancelProcess() { canceled = true; }
}
