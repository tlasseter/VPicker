package com.Sts.Types;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 *
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

/** reads a SEGY file in blocks into memory, converts samples to scaled bytes and outputs
 *  them into 3 byte cubes without headers ordered in line, crossline, and slice planes.
 *  Input sample data can be in any format supported by SEGY.
 */
public class StsSegyVolume extends StsSeismicBoundingBox implements Serializable
{
    /** Is the clipped data to be set to null or the maximum and minimum */
//	transient protected boolean clipToNull = false;
	/** null byte is always the 255th index */
	transient static byte NULLBYTE = (byte)255;
	/** list of start and end inline traces */
	transient public SeismicLine[] lines;
	/** Total samples in each of 255 steps */
	transient int dataCnt[] = new int[255];
	transient int ttlHistogramSamples = 0;

//	transient protected StsSEGYFormat segyFormat = null;

//	transient byte[] textHeader;
//	transient byte[] binaryHeader;
//	transient byte[] traceHeader;

	/** random access file for reading input segy file */
//	transient RandomAccessFile randomAccessSegyFile;
	transient FileChannel randomAccessSegyChannel;
	transient MappedByteBuffer inputBuffer = null;
	transient FileChannel inlineChannel;
//	transient MappedByteBuffer inlineBuffer = null;
	transient FileChannel inlineFloatChannel;
	transient MappedByteBuffer inlineFloatByteBuffer = null;
	transient FloatBuffer inlineFloatBuffer = null;
//	transient int fileHeaderSize = 0;
//	transient int bytesPerTrace;
//	transient long segyFileSize;

//	transient int nTotalTraces = 0;
	transient int nInputBytesPerBlock = 0;
	transient int nBlockTraces = 0;
	transient int nBlocksPerWrite = 0;
	transient int outputRowMin = 0;
	transient int nOutputSamplesPerWrite = 0;
	transient int nBlock = 0;
	transient boolean isLastBlock = false;
	/** number of blocks in memory waiting to be written out for xlines and slices */
	transient int nBlocksInMemory = 0;
	transient long inputPosition = 0;
	transient long outputPosition = 0;

	transient boolean scanHistogram = false;
	/** data scaling: 254/(dataMax - dataMin) */
	transient float scale = 1.0f;
	/** useful traces in cube */
	transient StsSEGYFormat.TraceHeader firstTraceHeader, lastTraceHeader;
	transient StsSEGYFormat.TraceHeader endFirstRowTraceHeader, endFirstColTraceHeader;

	/** wizard frame used for locating displayed dialogs */
	transient Frame frame;

    StsProgressPanel progressPanel;

    /** timer used in testing performance */
	transient StsTimer timer = null;
	transient StsTimer totalTimer = null;
	/** set to true if timing tests are desired */
	transient boolean runTimer = false;

	/** crop limits relative to row, col and slice limits */
	transient public StsCroppedBoundingBox cropBox;

	//transient StsSegyVolumeProcessPanel panel;

//	transient int textHeaderSize = 3200;
//	transient int binaryHeaderSize = 400;
//	transient int traceHeaderSize = 240;
//	transient int bytesPerSample = 4;

	// These four need to be persisted in database; i.e., moved up to StsSeismicBoundingBox
//    transient float verticalScalar = 1.0f;
//	transient float horizontalScalar = 1.0f;
//    transient String verticalUnits;
//	transient String horizontalUnits;

	/** 1=IBM, 2=short, 3=int, 5=IEEE, 8=byte; @see StsSEGYFormat */
//	transient int sampleFormat = 1;
//	transient float sampleSpacing = 0;
//	transient boolean getIsLittleEndian = false;

	transient static boolean isMinMemory = true;
	transient boolean displayMessages = true;

	transient RandomAccessFile rowFile = null, colFile = null, sliceFile = null;
	transient RandomAccessFile floatRowFile = null;
	transient RandomAccessFile attributeFile = null;
	transient byte[] blockColBytes = null, blockSliceBytes = null;
	transient StsMemAllocVolumeProcess memoryAllocation = null;

    transient StsSeismicBoundingBox overrideGeometryBox = null;
    transient int inlineDecimation = 1;
    transient int xlineDecimation = 1;

    transient boolean interrupt = false;

//    transient boolean progressBar = false;
	transient int progressInc = 0;
	transient StsProgressPanel pd = null;
	transient JDialog pdialog = null;
	transient double scanPercentage = 0.1;

	transient static final int INLINE = 0;
	transient static final int XLINE = 1;
	transient static final int TD = 2;

	transient byte[] segyTraceBytes;
	transient byte[] traceBytes;
	transient byte[] paddedTraceBytes;
	transient float[] nativePaddedTraceFloats;

    /** attributes for each trace; sequential for each attribute [nAttributes][nTraces] */
    transient double[][] traceAttributes;
    /** number of traces written to attribute array */
    transient int nAttributeTrace = 0;
    /** Attribute records for extraction */
    StsSEGYFormatRec[] attributeRecords = null;

	static public float defaultScanPercentage = 0.1f;

    static public final String group = StsSeismicBoundingBox.group3d;

    /** acceptable suffixes for segy files (upper or lower case ok) */
	static String[] segySuffixes = new String[] {".sgy", ".segy"};
	static final boolean debug = false;
	static final boolean memoryDebug = false;


	static final float nullValue = StsParameters.nullValue;

	public transient StsAnalyzeTraces analyzeTraces = null;

    public String getGroupname()
    {
        return group3d;
    }

    public StsSegyVolume()
	{}

	/**
	 *
	 * @param file segy input file
	 * @param stsDirectory directory where S2S seismic byte volumes are written
	 * @param frame frame used for centering this dialog
	 * @param segyFormat format information for reading segy file
	 */
	StsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws StsException
	{
		super(false);

		try
		{
			if(segyFormat == null)
			{
				throw new StsException("StsSegyVolume.constructTraceAnalyzer() failed. segyFormat cannot be null.");
			}
			if(runTimer) timer = new StsTimer();
            segyData = StsSegyData.constructor(seismicWizard, file, this, segyFormat);
            if(segyData == null)
            {
                throw new StsException("StsSegyVolume.constructor failed: couldn't construct StsSegyData.");
            }
            segyDirectory = file.getDirectory();
			segyFilename = file.getFilename();
			File sfile = new File(segyDirectory + segyFilename);
			segyLastModified = sfile.lastModified();
			this.stsDirectory = stsDirectory;
			this.frame = frame;
            stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes);
            setName(stemname);
//			getIsLittleEndian = segyFormat.getIsLittleEndian();
            this.zDomain = segyFormat.getZDomainString();
//            segyData.randomAccessSegyFile = new RandomAccessFile(file.getPathname(), "r");

//			readFileHeader();
//			analyzeBinaryHdr();
			zMin = 0.0f;
            setFilenames(stemname);
        }
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume(file, frame) failed.", e, StsException.WARNING);
		}
	}


    /*
		static public StsSegyVolume constructTraceAnalyzer(StsFile file, String stsDirectory, Frame frame)
		{
			try
			{
				return new StsSegyVolume(file, stsDirectory, frame, null);
			}
			catch (Exception e)
			{
				StsException.outputException("StsSegyVolume.constructTraceAnalyzer(file, frame) failed.", e, StsException.WARNING);
				return null;
			}
		}

		static public StsSegyVolume constructTraceAnalyzer(String segyDirectory, String stsDirectory, String filename, Frame frame)
		{
			return constructTraceAnalyzer(StsFile.constructTraceAnalyzer(segyDirectory, filename), stsDirectory, frame);
		}
	 */
	static public StsSegyVolume constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsSegyVolume(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.constructTraceAnalyzer(file, frame) failed.", e, StsException.WARNING);
			return null;
		}
	}
/*
	public void setParametersFromFile(String filename)
	{
		try
		{
			StsParameterFile.readObjectFields(filename, this, StsSegyVolume.class, StsBoundingBox.class);
			analyzeBinaryHdr(false);
			analyzeGeometry(false);
			System.out.println("Calling analyzeTraces(0) from StsSegyVolume for volume:" + this.getName());
			if(analyzeTraces(0) != StsAnalyzeTraces.ANALYSIS_OK)
				StsMessageFiles.errorMessage("StsSegyVolume:setParametersFromFile() failed, type appears to be IEEE"); ;
		}
		catch(Exception e)
		{
			StsMessageFiles.errorMessage("StsSegyVolume:setParametersFromFile() failed");
		}
	}
*/
	public String readTextHdr(String encoding)
	{
        return segyData.readTextHdr(encoding);
	}

	/*
		public String readTextHdr(String encoding)
		{
			try
			{
				return new String(textHeader, 0, segyFormat.textHeaderSize, encoding);
			}
			catch (Exception e)
			{
				StsException.outputException("StsSegyVolume.readTextHdr() failed.",
					e, StsException.WARNING);
				return null;
			}
		}
	 */
	public byte[] getTraceHeaderBinary(int nTrace)
	{
        return segyData.getTraceHeaderBinary(nTrace);
	}

    public String getGroup()
    {
        return group;
    }

    public StsSegyData getSegyData() { return segyData; }

	public void setNSamples(int nSamples)
	{
		nSlices = nSamples;
	}

	public int getNSamples()
	{
		return nSlices;
	}

    // Used to get the exact statistics of irregular volume. Only run if the user processes the volume.
	//
	// constructs a set of line descriptors used to analyze an irregular volume
	//

    public boolean analyzeIrregularGeometry()
	{
		StsSEGYFormat.TraceHeader firstRowTraceHeader, lastRowTraceHeader, prevFirstTrace;

		try
		{
			boolean geometryChanged = false;

			if(nRows < 0)
				return false;
			lines = new SeismicLine[nRows];

			// first line has been determined in analyzeGeometry
			firstRowTraceHeader = firstTraceHeader;
			lastRowTraceHeader = endFirstRowTraceHeader;
			lines[0] = new SeismicLine(firstRowTraceHeader, lastRowTraceHeader);
			int nRowTraces = nCols;
			prevFirstTrace = firstRowTraceHeader;
            progressPanel.initialize(nRows);
            progressPanel.appendLine("Analyzing irregular geometry for volume " + getName());
			for(int n = 1; n < nRows; n++)
			{
				int nFirstRowTrace = lastRowTraceHeader.nTrace + 1;
//                System.out.println("nRows=" + nRows + " nTraces=" + nFirstRowTrace);
				firstRowTraceHeader = getTrace(nFirstRowTrace);

				// fill in any blank lines
				int traceRow = getTraceRowIndex(firstRowTraceHeader);
				while(traceRow > n)
				{
					lines[n++] = new SeismicLine(null, null);
				}

				int nLastRowTrace = lastRowTraceHeader.nTrace + nRowTraces;
				lastRowTraceHeader = getTrace(nLastRowTrace);

				if(lastRowTraceHeader.iLine != firstRowTraceHeader.iLine)
					lastRowTraceHeader = searchBackwardForEndOfLine(nLastRowTrace, firstRowTraceHeader.iLine);
				else
					lastRowTraceHeader = searchForwardForEndOfLine(nLastRowTrace, firstRowTraceHeader.iLine);

				if(lastRowTraceHeader == null)
					return false;

				if(lineGeometryChanged(firstRowTraceHeader, lastRowTraceHeader))
					geometryChanged = true;
				if(xLineNumIncChanged(prevFirstTrace, firstRowTraceHeader))
					geometryChanged = true;
				lines[n] = new SeismicLine(firstRowTraceHeader, lastRowTraceHeader);
				prevFirstTrace = firstRowTraceHeader;
				nRowTraces = lastRowTraceHeader.nTrace - firstRowTraceHeader.nTrace + 1;
                progressPanel.setDescription("Analyzing irregular volume row " + n + " of " + (nRows-1));
                progressPanel.setValue(n+1);
			}
//            if(geometryChanged) return analyzeAngle();

			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.analyzeIrregularGeometry() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	private boolean lineGeometryChanged(StsSEGYFormat.TraceHeader firstRowTraceHeader,
										StsSEGYFormat.TraceHeader lastRowTraceHeader)
	{
		boolean changed = false;
		boolean findInc = false;

		float xLineNumMin = firstRowTraceHeader.xLine;
		float xLineNumMax = lastRowTraceHeader.xLine;

		// number of traces based on trace index
		int nLineIndexTraces = lastRowTraceHeader.nTrace - firstRowTraceHeader.nTrace + 1;

		// check if limits have changed
		if(nLineIndexTraces > 1 && (xLineNumMin != this.colNumMin || xLineNumMax != this.colNumMax))
		{
			changed = true;

			if(xLineNumMax > xLineNumMin)
			{
				this.colNumMin = Math.min(xLineNumMin, this.colNumMin);
				this.colNumMax = Math.max(xLineNumMax, this.colNumMax);
			}
			else
			{
				this.colNumMin = Math.max(xLineNumMin, this.colNumMin);
				this.colNumMax = Math.min(xLineNumMax, this.colNumMax);
			}

			// Check number of dead traces; if any found, recompute xLineNumInc.
			// If number of traces based on xLine numbering is not integral, we have a different
			//     xLineNumInc and there must be dead traces.
			//
			// If actual number of traces is less than number based on xline numbers, we have that many dead traces.
			// If actual number of traces is greater than number based on xLine numbers, then we have intermediate
			// traces which means there must be dead traces.

			// number of traces based on xline number
			float nLineNumTracesF = (colNumMax - colNumMin) / colNumInc + 1.0f;
			int nLineNumTraces = Math.round(nLineNumTracesF);

			// check if we have a fractional trace spacing at either end
			float fraction = Math.abs(nLineNumTracesF - nLineNumTraces);
			findInc = fraction > 0.01f;

			// check estimated number of dead traces and find increment if necessary
			if(!findInc)
			{
				int nEstDeadTraces = nLineNumTraces - nLineIndexTraces;
				findInc = nEstDeadTraces < 0 || (nEstDeadTraces > 0 && nEstDeadTraces > 0.1f * nLineIndexTraces);
			}
		}

		// check if xLine increment (xLineNumInc) has possibly changed.
		// If number of traces computed from xline number range (nLineNumTraces) is different than
		// number of traces computed from index range (nLineIndexTraces), we have dead traces.
		// Search the line to correct a possible change in xLineNumInc
		if(findInc)
		{
			StsSEGYFormat.TraceHeader prevRowTraceHeader, nextRowTraceHeader;

			int nFirstRowTrace = firstRowTraceHeader.nTrace;
			int nLastRowTrace = lastRowTraceHeader.nTrace;
			nextRowTraceHeader = firstRowTraceHeader;
			for(int nTrace = nFirstRowTrace; nTrace < nLastRowTrace; nTrace++)
			{
				prevRowTraceHeader = nextRowTraceHeader;
				nextRowTraceHeader = getTrace(nTrace + 1);
				float newInc = nextRowTraceHeader.xLine - prevRowTraceHeader.xLine;
				if(Math.abs(newInc) < Math.abs(colNumInc) && newInc != 0)
				{
					float ratio = colNumInc / newInc;
					float error = Math.abs(ratio - (int)ratio);
					if(error > 0.01)
					{
						StsMessageFiles.errorMessage("Crossline numbering inconsistent between trace[" +
							prevRowTraceHeader.iLine +
							"," + prevRowTraceHeader.xLine + "]" + " and trace[" +
							nextRowTraceHeader.iLine + "," +
							nextRowTraceHeader.xLine + "]. Ignoring it.");
					}
					else
					{
						colNumInc = newInc;
						changed = true;
					}
				}
			}
		}
		if(changed)
            nCols = Math.max(nCols, Math.round((colNumMax - colNumMin) / colNumInc) + 1);
		return changed;
	}

	private boolean xLineNumIncChanged(StsSEGYFormat.TraceHeader prevRowTraceHeader,
									   StsSEGYFormat.TraceHeader nextRowTraceHeader)
	{
		boolean changed = false;

		float newInc = nextRowTraceHeader.iLine - prevRowTraceHeader.iLine;
		if(Math.abs(newInc) < Math.abs(rowNumInc) && newInc != 0)
		{
			float ratio = rowNumInc / newInc;
			float error = Math.abs(ratio - (int)ratio);
			if(error < 0.01)
			{
				StsMessageFiles.errorMessage("Line numbering inconsistent between trace[" + prevRowTraceHeader.iLine +
											 "," + prevRowTraceHeader.xLine + "]" + " and trace[" +
											 nextRowTraceHeader.iLine + "," +
											 nextRowTraceHeader.xLine + "]. Ignoring it.");
			}
			else
			{
				rowNumInc = newInc;
				nRows = Math.round((rowNumMax - rowNumMin) / rowNumInc) + 1;
				changed = true;
			}
		}
		return changed;
	}

	public void setIsRegular(boolean value)
	{
		isRegular = value;
	}

	public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages)
	{
		displayMessages = messages;
        this.progressPanel = progressPanel;
        return analyzeGeometry();
	}

    private boolean analyzeGeometry()
	{
		try
		{
            originSet = false;
 /*
            if(useKnownPoints)
			{
				return analyzeGeometryWithKnownPoints(segyData.getSegyFormat());
			}
 */
			StsSEGYFormat.TraceHeader secondTrace, trace;
			int nTrace;

			isRegular = true;
        /*
            bytesPerSample = StsSegyIO.getBytesPerSample(sampleFormat);
			traceHeaderSize = segyFormat.traceHeaderSize;
			traceHeader = new byte[traceHeaderSize];
			bytesPerTrace = nCroppedSlices * bytesPerSample + traceHeaderSize;
			segyFileSize = randomAccessSegyFile.length();
			nTotalTraces = (int)((segyFileSize - fileHeaderSize) / bytesPerTrace);
	    */
		//	if(nTotalTraces == 0) return false;
//            float nTotalTracesFloat = (float)((float)(segyFileSize - fileHeaderSize) / (float)bytesPerTrace);
//            if(nTotalTracesFloat != nTotalTraces)
//            {
//                StsException.outputException("Unexpected file size, please verify header sizes.", null, StsException.WARNING);
//                return false;
//            }
			if(nSlices == 0) return false;
//            zMin = segyFormat.startZ * verticalScalar;
			zMax = zMin + (nSlices - 1) * zInc;

			long ioff = segyData.getFileHeaderSize();
            int nTotalTraces = segyData.getNTotalTraces();
			// Get first trace to set line and xLine number minimum
			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null)
				return false;

			rowNumMin = firstTraceHeader.iLine;
			colNumMin = firstTraceHeader.xLine;

			// Get last trace to compute line and xLine number maximums
			lastTraceHeader = getTrace(nTotalTraces - 1);
			if(lastTraceHeader == null)
				return false;

			rowNumMax = lastTraceHeader.iLine;
			colNumMax = lastTraceHeader.xLine;

			// Sanity check the last trace to determine if format and/or line numbering is correct
			// If not, don't bother with analysis
// SAJ - This is a problem with a file I'm working and I don't understand the logic --- Hopefully not needed.
//            if((Math.abs(rowNumMax - rowNumMin) > segyFileSize/(nCroppedSlices+traceHeaderSize)) ||
//               (Math.abs(colNumMax - colNumMin) > segyFileSize/(nCroppedSlices+traceHeaderSize)))
//            {
//                new StsMessage(frame, StsMessage.WARNING, "Line numbers not in specified locations and/or format size wrong.");
//                return false;
//            }

			if(rowNumMax == rowNumMin)
            {
                progressPanel.appendLine("First and last trace have same inline number: " + rowNumMin);
                nCols = 24;
                nRows = nTotalTraces/24;
                return false;
            }
            if(colNumMax == colNumMin)
			{
               progressPanel.appendLine("First and last trace have same xline number: " + colNumMin);
               nCols = 24;
               nRows = nTotalTraces/24;
               return false;
			}
			// search a few traces to compute xLine number increment
			if(!estimateXLineNumInc(progressPanel)) return false;

			nCols = Math.round((colNumMax - colNumMin) / colNumInc) + 1;
			if(nCols <= 1)
			{
				isRegular = false;
				nCols = 2;
			}

			// Get what we think is end of first line (if not on same line, set isRegular to false)
			endFirstRowTraceHeader = getTrace(nCols - 1);
			if(endFirstRowTraceHeader == null)
				return false;
			if(endFirstRowTraceHeader.iLine != firstTraceHeader.iLine) // we must be on one of the next iLines
			{
				isRegular = false;
				endFirstRowTraceHeader = searchBackwardForEndOfLine(nCols - 1, firstTraceHeader.iLine);
				if(endFirstRowTraceHeader == null)
					return false;
			}
			else
			{
				// check if next trace after endFirstRowTraceHeader is on next line; if on same line,
				// search forward for end of line
				trace = getTrace(nCols);
				if(trace.iLine == firstTraceHeader.iLine)
				{
					isRegular = false;
					endFirstRowTraceHeader = searchForwardForEndOfLine(nCols, firstTraceHeader.iLine);
				}
				nCols = Math.max(nCols, endFirstRowTraceHeader.nTrace + 1);
			}

			// get first trace on next line
			trace = getTrace(endFirstRowTraceHeader.nTrace + 1);
			rowNumInc = trace.iLine - firstTraceHeader.iLine;
			nRows = Math.round((rowNumMax - rowNumMin) / rowNumInc) + 1;

			// check file size to make sure its really regular
			if(isRegular && nRows * nCols != nTotalTraces)
				isRegular = false;

			// if irregular, reestimate xLine max using the first line in addition to lastTrace above
			if(!isRegular)
			{
				if(colNumInc > 0)
					colNumMax = Math.max(colNumMax, endFirstRowTraceHeader.xLine);
				else
					colNumMax = Math.min(colNumMax, endFirstRowTraceHeader.xLine);
			}

			// get trace at start of last line which is the end of the first xLine for
			// a regular cube
			nTrace = nTotalTraces - nCols;
			endFirstColTraceHeader = getTrace(nTrace);

			// compute rotation angle and bin spacings in line and xline directions
//			if(!analyzeAngle()) return false;

//            if(!isXLineCCW)
//                flipLineNumOrder();

			xOrigin = firstTraceHeader.x * horizontalScalar;
			yOrigin = firstTraceHeader.y * horizontalScalar;

            if(!isRegular) analyzeIrregularGeometry();
//			initializeAngle();
//			initializeRange();
//            initializeCropBoxRange();
//            progressPanel.appendLine("Completed volume " + getName());
//            progressPanel.setDescription("Completed volume " + getName());
//            progressPanel.finished();
            return true;
		}
		// again just use (Exception e) and you get everything
		catch(Exception e)
		{
			if(displayMessages)
				StsMessageFiles.errorMessage("Couldn't analyze geometry: select a different number of bytes per sample.");
			return false;
		}
	}

	public boolean analyzeAngle(StsProgressPanel progressPanel, boolean messages)
	{
		displayMessages = messages;
        this.progressPanel = progressPanel;
        return analyzeAngle();
	}

    private boolean analyzeAngle()
	{
 //       if(useKnownPoints)
 //           return analyzeGeometryWithKnownPoints(segyData.getSegyFormat());
        if(super.analyzeAngle(firstTraceHeader, endFirstRowTraceHeader, lastTraceHeader))
			return true;
		if(super.analyzeAngle(firstTraceHeader, endFirstColTraceHeader, lastTraceHeader))
			return true;
		return false;
    }  

    public void progressPanelAppendTraceHeaderDescriptions(StsProgressPanel panel)
    {
        panel.appendLine("Trace header info for volume: " + getName());
        if(firstTraceHeader != null) panel.appendLine("    First trace Header: " + firstTraceHeader.toString());
        if(endFirstRowTraceHeader != null) panel.appendLine("    First inline last trace Header: " + endFirstRowTraceHeader.toString());
        if(endFirstColTraceHeader != null) panel.appendLine("    First xline last trace Header: " + endFirstColTraceHeader.toString());
        if(lastTraceHeader != null) panel.appendLine("    Last trace Header: " + lastTraceHeader.toString());
    }

    public void setCropBox(StsCroppedBoundingBox cropBox)
    {
        this.cropBox = cropBox;
    }
/*
    private void initializeCropBoxRange()
    {
        cropBox.rowMin = 0;
		cropBox.rowMax = (nRows - 1);
		cropBox.rowInc = 1;
		cropBox.colMin = 0;
		cropBox.colMax = (nCols - 1);
		cropBox.colInc = 1;
		cropBox.sliceMin = 0;
		cropBox.sliceMax = (nCroppedSlices - 1);
		cropBox.sliceInc = 1;
	}
*/
    /** Searching from nStartTrace to nEndTrace, find xLineNumInc.
	 *  If there are dead traces, they may be skipped, so find at least
	 *  10 sequential trace pairs that agree on this number.
	 */
	private boolean estimateXLineNumInc(StsProgressPanel progressPanel)
	{
		StsSEGYFormat.TraceHeader prevTrace, nextTrace;
		float inc = StsParameters.largeFloat;

        int nTotalTraces = segyData.getNTotalTraces();
        int nTraces = nTotalTraces / 10;
		int nEndTrace = Math.min(nTotalTraces - 1, nTraces);
        int minAgreeablePairs = Math.min(10, nEndTrace);
        int nAgreeableTracePairs = 0;
		nextTrace = getTrace(0);
		for(int n = 1; n <= nEndTrace; n++)
		{
			prevTrace = nextTrace;
			nextTrace = getTrace(n);
			if(prevTrace.iLine == nextTrace.iLine) // traces are sequential on same line
			{
				float newInc = nextTrace.xLine - prevTrace.xLine;
                if(newInc <= 0)
                {
                    if(progressPanel != null)
                    {
                       progressPanel.appendLine("Adjacent traces on line are not in increasing crossline order:");
                       progressPanel.appendLine("    " + prevTrace.toString());
                       progressPanel.appendLine("    " + nextTrace.toString());
                       return false;
                    }
                }

                if(Math.abs(newInc) < Math.abs(inc))
				{
					inc = newInc;
					nAgreeableTracePairs = 1;
				}
				else if(Math.abs(newInc) == Math.abs(inc))
				{
					nAgreeableTracePairs++;
                    if(nAgreeableTracePairs >= minAgreeablePairs)
					{
						colNumInc = inc;
						return true;
					}
				}
			}
		}
        if(progressPanel != null) progressPanel.appendLine("Unable to determine crossline increment. Check crossline numbering.");
        return false;
	}

	StsSEGYFormat.TraceHeader getTrace(int n)
	{
        return segyData.getTrace(n);
	}

	private StsSEGYFormat.TraceHeader searchForwardForEndOfLine(int nTrace, float iLine)
	{
		try
		{
            int nTotalTraces = segyData.getNTotalTraces();
            if(nTrace >= nTotalTraces)
				return searchBackwardForEndOfLine(nTotalTraces, iLine);

			StsSEGYFormat.TraceHeader startTrace = getTrace(nTrace);
			if(startTrace == null)
				return null;
			StsSEGYFormat.TraceHeader prevTrace = startTrace;
			while(nTrace <= nTotalTraces - 1)
			{
				nTrace++;
				StsSEGYFormat.TraceHeader trace = getTrace(nTrace);
				if(trace == null)
					return null;
				if(trace.iLine != iLine)
					return prevTrace;
				if(nTrace == nTotalTraces)
					return prevTrace;
				prevTrace = trace;
			}
			return null;
		}
		catch(Exception e)
		{
			StsException.systemError("StsSegyVolume.searchForEndOfLine() failed for trace number : " + nTrace);
			return null;
		}
	}

	private StsSEGYFormat.TraceHeader searchBackwardForEndOfLine(int nTrace, float iLine)
	{
		try
		{
			if(nTrace <= 0)
				return searchForwardForEndOfLine(0, iLine);

			while(nTrace > 0)
			{
				nTrace--;
				StsSEGYFormat.TraceHeader trace = getTrace(nTrace);
				if(trace.iLine == iLine)
					return trace;
			}
			return null;
		}
		catch(Exception e)
		{
			StsException.systemError("StsSegyVolume.searchForwardForLine() failed for iLine: " + iLine);
			return null;
		}
	}

	public void resetDataRange()
	{
        segyData.resetDataRange();
	}

	// JKF 04JUNE2006
/*
    public boolean analyzeTraces(StsProcessStepQueue processStepQueue, double scanPercentage) throws IOException
	{
		try
		{
			this.processStepQueue = processStepQueue;
			this.scanPercentage = scanPercentage;
			scanHistogram = true;
			analyzeTraces();
			return true;
		}
		catch(Exception e)
		{
			return false;
		}
	}

	public boolean analyzeTraces(StsProgressPanel panel, double scanPercentage) throws IOException
	{
		try
		{
			this.pd = panel;
			this.scanPercentage = scanPercentage;
			scanHistogram = true;
			analyzeTraces();
			return true;
		}
		catch(Exception e)
		{
			return false;
		}
	}
*/
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
            if( result == StsSegyIO.ANALYSIS_OK)
                calculateHistogram(segyIO.histogramSamples, segyIO.nHistogramSamples);
            return result;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackSegyLine.analyzeTraces() failed.", e, StsException.WARNING);
			return StsAnalyzeTraces.ANALYSIS_UNKNOWN_EXCEPTION;
		}
	}

/*
	private void runAnalyzeTraces()
	{
        Runnable runnable = new Runnable()
		{
			public void run()
			{
				analyzeTraces();
			}
		};
		Thread thread = new Thread(runnable);
		thread.start();
	}
*/
/*
    private boolean analyzeTraces()
	{
		try
		{
			boolean fmtOkay = true;
//			progressValue = 0.0;
			interrupt = false;
			initializeSegyIO();

			if(segyFormat != null)
			{
				resetDataRange();
			}

			//
			// If irregular, must determine actual extent prior to crop screen
			//
			if(!isRegular) // && !scanHistogram)
			{
				progressDescription = "Analyzing irregular geometry for " + getName();
				if(pd != null)
					pd.appendLine("Analyzing irregular geometry for " + getName());
				analyzeIrregularGeometry();
			}
			else
				progressDescription = "";

            int ntraces = nRows * nCols;
			int nTracesToScan = Math.max(((int)((float)nTotalTraces * (scanPercentage / 100.0f))), 10); // Default to ~1% of traces

			// Set to run entire file if scan percentage is greater than 50%, will be faster
			if(scanPercentage > 50.0f && ntraces > 10)
				nTracesToScan = ntraces;

            // Set to run entire volume if small enough  (<10 Mbytes).
            if(randomAccessSegyFile.length() < 10000000)
                nTracesToScan = ntraces;

            analyzeTraces = StsAnalyzeTraces.constructTraceAnalyzer(randomAccessSegyFile,
				fileHeaderSize,
				nTracesToScan,
				nTotalTraces,
				nCroppedSlices,
				bytesPerSample,
				traceHeaderSize,
				getIsLittleEndian,
				sampleFormat,
				segyFormat.getTraceRecFromUserName("B2SCALCO"),
				segyFormat.getTraceRecFromUserName("TRC_TYPE"));

			analyzeTraces.analyze();
			fmtOkay = true;
	//		fmtOkay = analyzeTraces.getFmtOkay();
			if(fmtOkay)
			{
				dataMin = analyzeTraces.getDataMin();
				dataMax = analyzeTraces.getDataMax();

				dataAvg = (float)analyzeTraces.getAverage();
				if(scanHistogram)
				{
					dataCnt = analyzeTraces.getHistogramData();
					ttlHistogramSamples = analyzeTraces.getTtlHistogramSamples();
					calculateHistogram();
				}
			}
			return fmtOkay;
		}
		catch(Exception ex)
		{
			progressDescription = ex.getMessage();
			return false;
		}
	}
*/
	public long getOutputFileSize()
	{
		long traceDerived = (long)nRows * (long)nCols * (long)nSlices;
		long fileDerived = (segyData.segyFileSize - (long)segyData.binaryHeaderSize - (long)segyData.textHeaderSize -
							((long)nRows * (long)nCols * (long)segyData.traceHeaderSize)) / segyData.bytesPerSample;
		if(traceDerived < fileDerived)
			return fileDerived;
		else
			return traceDerived;
	}
/*
	public void setFilenames()
	{
		if(groupname != null)
			setFilenames(groupname + "." + stemname);
		else
			setFilenames(stemname);
	}
*/
	public void setFilenames(String stemname)
	{
		createVolumeFilenames(stemname);
		getCreateAttributeFilename(stemname);
		createFloatRowVolumeFilename(stemname);
	}

	private void resetSeismicBoundingBoxParametersToCrop()
	{
		nRows = cropBox.nRows;
		nCols = cropBox.nCols;
		nSlices = cropBox.nSlices;

		// if the first xline is not CCW from inline, the origin needs to be
		// at the other end of the first xline so xline is then CCW from inline
//        if(!isXLineCCW)
//            flipLineNumOrder();

		// origin is temporarily at location of first trace
		// move it to the lower left corner of the cube
/*
		int originalOriginRow, originalOriginCol;
		if(useKnownPoints)
		{
			originalOriginRow = Math.round((knownPoints[0].inline - rowNumMin) / rowNumInc);
			originalOriginCol = Math.round((knownPoints[0].crossline - colNumMin) / colNumInc);
		}
		else
		{
			originalOriginRow = Math.round((firstTraceHeader.iLine - rowNumMin) / rowNumInc);
			originalOriginCol = Math.round((firstTraceHeader.xLine - colNumMin) / colNumInc);
		}

		int originRow = cropBox.rowMin;
		int originCol = cropBox.colMin;

		// adjustments from original origin to new origin
		float dX = (originCol - originalOriginCol) * xInc;
		float dY = (originRow - originalOriginRow) * yInc;

		double[] origin = this.getAbsoluteXY(dX, dY);
		xOrigin = origin[0];
		yOrigin = origin[1];
*/
        xMax = xMin + cropBox.colMax*xInc;        
        xMin += cropBox.colMin*xInc;
		xInc = cropBox.colInc * xInc;

		yMax = yMin + cropBox.rowMax * yInc;
        yMin += cropBox.rowMin*yInc;
        yInc = cropBox.rowInc * yInc;

        zMax = zMin + cropBox.sliceMax*zInc;
        zMin += (cropBox.sliceMin * zInc);
		zInc = cropBox.sliceInc * zInc;

		rowNumMin = (cropBox.rowMin * rowNumInc) + rowNumMin;
		rowNumMax = ((cropBox.rowMax - cropBox.rowMin) * rowNumInc) + rowNumMin;
		rowNumInc *= cropBox.rowInc;
		colNumMin = (cropBox.colMin * colNumInc) + colNumMin;
		colNumMax = ((cropBox.colMax - cropBox.colMin) * colNumInc) + colNumMin;
		colNumInc *= cropBox.colInc;

//		if(!isXLineCCW) flipRowNumOrder();
	}

	/** Flip inLine numbering if xLine is CW from line as
	 *  origin has been moved to lower-left from upper-left.
	 */
	private void flipLineNumOrder()
	{
		float temp = rowNumMin;
		rowNumMin = rowNumMax;
		rowNumMax = temp;
		rowNumInc = -rowNumInc;
		temp = yMin;
		yMin = yMax;
		yMax = temp;
		yInc = -yInc;
	}

	public boolean readWriteVolume(StsModel model, StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecs)
	{
        segyIO.initializeDataRange(dataMin, dataMax);
		nAttributes = attributeRecs.length;
		attributeRecords = attributeRecs;
        new StsPostStackSegyVolumeConstructor(model, this, progressPanel);
        return true;
    }

	public int getTraceColIndex(StsSEGYFormat.TraceHeader traceHeader)
	{
		return Math.round((traceHeader.xLine - colNumMin) / colNumInc);
	}

	public int getTraceColIndex(byte[] segyTraceBytes)
	{
		float xLine = segyData.getXLine(segyTraceBytes);
		return Math.round((xLine - colNumMin) / colNumInc);
	}

	private int getTraceRowIndex(StsSEGYFormat.TraceHeader traceHeader)
	{
		return Math.round((traceHeader.iLine - rowNumMin) / rowNumInc);
	}

	public void cancel()
	{
		interrupt = true;
	}

	public void close()
	{
		try
		{
			memoryAllocation.freeMemory();
		}
		catch(Exception e)
		{}

		if(segyData.randomAccessSegyFile != null)
			try
			{
				segyData.randomAccessSegyFile.close();
			}
			catch(Exception e)
			{}
		if(rowFile != null)
			try
			{
				rowFile.close();
			}
			catch(Exception e)
			{}
		if(colFile != null)
			try
			{
				colFile.close();
			}
			catch(Exception e)
			{}
		if(sliceFile != null)
			try
			{
				sliceFile.close();
			}
			catch(Exception e)
			{}
		if(attributeFile != null)
			try
			{
				attributeFile.close();
			}
			catch(Exception e)
			{}
		if(floatRowFile != null)
			try
			{
				floatRowFile.close();
			}
			catch(Exception e)
			{}

	}

    public String getLabel()
	{
		return stemname;
	}

	public void setIsMinMemory(boolean isMinMemory)
	{
		this.isMinMemory = isMinMemory;
	}

	public boolean isMinMemory()
	{
		return isMinMemory;
	}

	public float getTInc()
	{
		return zInc;
	}

	public void resetDataMinMax()
	{
		resetDataMinMaxByFormat(segyData.getSampleFormat());
	}

	public void setDataMin()
	{
		setDataMinByFormat(segyData.getSampleFormat());
	}

	public void setDataMax()
	{
		setDataMaxByFormat(segyData.getSampleFormat());
	}

	// Set geom from user specified min/max/inc values:
	public void setILineGeom(int iLineNumMin, int iLineNumMax, int iLineNumInc) throws StsException
	{
		setRowNumMin((float)iLineNumMin);
		setRowNumMax((float)iLineNumMax);
		setRowNumInc((float)iLineNumInc);
		float nRowsF = (rowNumMax - rowNumMin) / rowNumInc;
		setNRows(Math.round(nRowsF));
	}

	public void setXLineGeom(int xLineNumMin, int xLineNumMax, int xLineNumInc) throws StsException
	{
		setColNumMin((float)xLineNumMin);
		setColNumMax((float)xLineNumMax);
		setColNumInc((float)xLineNumInc);
		float nColsF = (colNumMax - colNumMin) / colNumInc;
		setNCols(Math.round(nColsF));
	}

	public class SeismicLine
	{
		public StsSEGYFormat.TraceHeader firstTraceHeader, lastTraceHeader;

		public SeismicLine(StsSEGYFormat.TraceHeader firstTraceHeader, StsSEGYFormat.TraceHeader lastTraceHeader)
		{
			this.firstTraceHeader = firstTraceHeader;
			this.lastTraceHeader = lastTraceHeader;
		}

		public String toString()
		{
			return " SeismicLine from " + firstTraceHeader.toString() + " to " + lastTraceHeader.toString();
		}
	}

	public void initializeAttributes()
	{
		int nAttributes = attributeRecords.length;
		attributeNames = new String[nAttributes];
		for(int n = 0; n < nAttributes; n++)
		{
			StsSEGYFormatRec record = attributeRecords[n];
			String userName = record.getUserName();
			attributeNames[n] = userName;
		}
	}

	private void outputAttributeBlock(int firstRow, int lastRow)
	{
		StsMappedDoubleBuffer attributeBuffer = null;
		try
		{
            getCreateAttributeFilename();
            attributeBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "rw");
			if(attributeBuffer == null)
				return;

			// Write Attributes to file
//			System.out.println("Writing attributes to file.");
			int offset = 0;
			long samplesPerSlice = (long)cropBox.nSamplesPerSlice;
			long startPosition = (long)firstRow * (long)cropBox.nCols;
			int tracesInBlock = cropBox.nBlockSliceSamples;
//			long tracesInBlock = ( (long) lastRow - (long) firstRow + 1) * (long) cropBox.nCols;

			for(int j = 0; j < attributeRecords.length; j++)
			{
				long position = (long)j * samplesPerSlice + startPosition;
				attributeBuffer.map(position, tracesInBlock);
				attributeFile.seek(position);
				attributeBuffer.put(traceAttributes[j], 0, tracesInBlock);
				attributeBuffer.clean();
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.outputBlock.process() failed.", e, StsException.WARNING);
			return;
		}
		finally
		{
			if(attributeBuffer != null)
				attributeBuffer.close();
		}
	}

	private void outputAttributeBlockX(int firstRow, int lastRow)
	{
		// Write Attributes to file...................SAJ
		System.out.println("Writing attributes to file.");
		int offset = 0;
		long bytesPerSlice = (long)cropBox.nSamplesPerSlice * (long)8;
		long startPosition = (long)firstRow * (long)cropBox.nCols * (long)8;
		long samplesInBlock = ((long)lastRow - (long)firstRow + 1) * (long)cropBox.nCols;
		try
		{
			for(int j = 0; j < attributeRecords.length; j++)
			{
				long position = (long)j * bytesPerSlice + startPosition;
				attributeFile.seek(position);
				for(int i = 0; i < samplesInBlock; i++)
					attributeFile.writeDouble(traceAttributes[j][i]);
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.outputBlock.process() failed.", e, StsException.WARNING);
			return;
		}
	}

	public void setVerticalUnits(byte vUnits)
	{
		StsParameters.getDistanceUnitString(vUnits);
	}

    public void setOverrideGeometry(StsSeismicBoundingBox boundingBox)
	{
        xInc = boundingBox.xInc;
        yInc = boundingBox.yInc;
        angle = boundingBox.angle;
        
        this.overrideGeometryBox = boundingBox;
	}

	public void setInlineDecimation(int value)
	{
		inlineDecimation = value;
	}

	public void setXlineDecimation(int value)
	{
		xlineDecimation = value;
	}
/*
	class CropBox
	{
		// describes the cropped box in terms of row,col,slice and increments
		int rowMin = 0;
		int rowMax = 0;
		int rowInc = 0;
		int nRows = 0;
		int colMin = 0;
		int colMax = 0;
		int colInc = 0;
		int nCols = 0;
		int sliceMin = 0;
		int sliceMax = 0;
		int sliceInc = 0;
		int nCroppedSlices = 0;

		// number of samples in a row plane of the output cube
		int nSamplesPerRow = 0;
		// number of samples in a col plane of the output cube
		int nSamplesPerCol = 0;
		// number of samples in a slice plane of the output cube
		int nSamplesPerSlice = 0;

		// row number of the first and last rows in the current block being processed
		int blockRowMin = 0;
		int blockRowMax = 0;

		// number of rows in this block
		int nBlockRows = 0;

		// convenience values for the block
		int nBlockColSamples;
		int nBlockSliceSamples;
		int nBlockSamples;

		int outputRowMin = 0;
		int outputRowMax = -1;

		CropBox()
		{
		}

		void initialize(int nUncroppedRows, int nUncroppedCols, int nUncroppedSlices)
		{
			// Make sure range is divisible by increment
			while(((rowMax - rowMin) % rowInc) != 0)
				rowMax--;
			nRows = 1 + Math.round((rowMax - rowMin) / rowInc);

			while(((colMax - colMin) % colInc) != 0)
				colMax--;
			nCols = 1 + Math.round((colMax - colMin) / colInc);

			while(((sliceMax - sliceMin) % sliceInc) != 0)
				sliceMax--;
			nCroppedSlices = 1 + Math.round((sliceMax - sliceMin) / sliceInc);

			// Determine if the volume is to be cropped
			if((rowMin > 0) || (rowMax < (nUncroppedRows - 1)) || (rowInc != 1))
				isCropped = true;
			else if((colMin > 0) || (colMax < (nUncroppedCols - 1)) || (colInc != 1))
				isCropped = true;
			else if((sliceMin > 0) || (sliceMax < (nUncroppedSlices - 1)) || (sliceInc != 1))
				isCropped = true;

			nSamplesPerRow = nCols * nCroppedSlices;
			nSamplesPerCol = nRows * nCroppedSlices;
			nSamplesPerSlice = nRows * nCols;
		}

		void adjustCropBlockRows(int uncropBlockRowMin, int uncropBlockRowMax)
		{
			if(isCropped)
			{
				blockRowMin = StsMath.intervalRoundUp(uncropBlockRowMin, rowInc);
				blockRowMax = StsMath.intervalRoundDown(uncropBlockRowMax, rowInc);
			}
			else
			{
				blockRowMin = uncropBlockRowMin;
				blockRowMax = uncropBlockRowMax;
			}

			nBlockRows = Math.max(0, (blockRowMax - blockRowMin) / rowInc + 1);
			nBlockSamples = nBlockRows * nSamplesPerRow;
			nBlockColSamples = nBlockRows * nCroppedSlices;
			nBlockSliceSamples = nBlockRows * nCols;

			outputRowMin = outputRowMax + 1;
			outputRowMax = outputRowMin + nBlockRows - 1;

			if(traceAttributes == null || traceAttributes[0].length != nBlockSliceSamples)
				traceAttributes = new double[nAttributes][cropBox.nBlockSliceSamples];
		}
	}
*/
}