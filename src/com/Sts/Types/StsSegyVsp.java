package com.Sts.Types;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
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
public class StsSegyVsp extends StsSeismicBoundingBox implements Serializable
{
    /** Is the clipped data to be set to null or the maximum and minimum */
//	transient protected boolean clipToNull = false;

    transient boolean multiComponent = false;
    transient int traceId = -1;
    transient String componentString = "";
//	transient protected StsSEGYFormat segyFormat = null;
    transient int actualTraces = 0;

//	transient byte[] textHeader;
//	transient byte[] binaryHeader//	transient byte[] traceHeader;

	/** random access file for reading input segy file */
//	transient RandomAccessFile randomAccessSegyFile;
	transient FileChannel randomAccessSegyChannel;
	transient MappedByteBuffer inputBuffer = null;
//	transient FileChannel inlineChannel;
//	transient MappedByteBuffer inlineBuffer = null;
//	transient FileChannel inlineFloatChannel;
//	transient MappedByteBuffer inlineFloatByteBuffer = null;
//	transient FloatBuffer inlineFloatBuffer = null;
//	transient int fileHeaderSize = 0;
//	transient int bytesPerTrace;

    transient boolean useKnownPoints = false;

	transient int nTraces = 0;
	transient long inputPosition = 0;
	transient long outputPosition = 0;

	transient boolean scanHistogram = false;
	/** data scaling: 254/(dataMax - dataMin) */
	transient float scale = 1.0f;
	/** useful traces in cube */
	transient StsSEGYFormat.TraceHeader firstTraceHeader, lastTraceHeader;
	transient StsSEGYFormat.TraceHeader endFirstRowTraceHeader, endFirstColTraceHeader;

	transient boolean displayMessages = true;

	/** wizard frame used for locating displayed dialogs */
	transient Frame frame;

	/** timer used in testing performance */
	transient StsTimer timer = null;
	transient StsTimer totalTimer = null;
	/** set to true if timing tests are desired */
	transient boolean runTimer = false;

    /** list of trace types detected dirung scan */
    transient int[] traceTypes = null;

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

	transient StsBinaryFile rowFile = null;
	transient StsBinaryFile floatRowFile = null;
	transient RandomAccessFile attributeFile = null;
	transient byte[] blockColBytes = null, blockSliceBytes = null;
//	transient StsMemAllocVolumeProcess memoryAllocation = null;

    transient StsSEGYFormatRec[] selectedAttributes = null;

	transient boolean interrupt = false;

	transient boolean progressBar = false;
	transient int progressInc = 0;
	transient StsProgressPanel pd = null;
	transient JDialog pdialog = null;
	transient double scanPercentage = 10.0;

	transient static final int INLINE = 0;
	transient static final int XLINE = 1;
	transient static final int TD = 2;

	transient byte[] segyTraceBytes;
	transient byte[] traceBytes;
	transient byte[] paddedTraceBytes;
	transient float[] nativePaddedTraceFloats;

	/** temp float array used for output trace */
	transient float[] outputTraceFloats;
    /** temp byte array used for output trace */
    transient byte[] outputTraceBytes;
    /** attributes for each trace; sequential for each attribute [nAttributes][nTraces] */
	transient double[][] traceAttributes;
	/** number of traces written to attribute array */
	transient int nAttributeTrace = 0;
	/** Attribute records for extraction */
	StsSEGYFormatRec[] attributeRecords = null;

	StsProgressPanel panel;

	private transient StsAnalyzeTraces analyzeTraces = null;

	static public float defaultScanPercentage = 10.0f;

	public final static byte nullByte = StsParameters.nullByte;

	/** acceptable suffixes for segy files (upper or lower case ok) */
	static String[] segySuffixes = new String[] { ".sgy", ".segy" };

	static final boolean debug = false;
	static final boolean memoryDebug = false;

	static final float nullValue = StsParameters.nullValue;

    public String getGroupname() { return groupVsp; }

    public StsSegyVsp() {}

	/**
	 *
	 * @param file segy input file
	 * @param stsDirectory directory where S2S seismic byte volumes are written
	 * @param frame frame used for centering this dialog
	 * @param segyFormat format information for reading segy file
	 */
	StsSegyVsp(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws StsException
	{
		super(false);

		try
		{
			if(segyFormat == null)
			{
				throw new StsException("StsSegyVsp.constructor failed: segyFormat cannot be null.");
			}
			if(runTimer) timer = new StsTimer();

            segyData = StsSegyData.constructor(seismicWizard, file, this, segyFormat);
            if(segyData == null)
            {
                throw new StsException("StsSegyVsp.constructor failed: couldn't construct StsSegyData.");
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
//			segyData.analyzeBinaryHdr();
			zMin = 0.0f;
            setFilenames();
            StsSEGYFormatRec[] records = segyFormat.getAllTraceRecords();
            selectedAttributes = new StsSEGYFormatRec[segyFormat.getNumberOfRequiredRecords()];
            int cnt = 0;
            for (int i = 0; i < records.length; i++)
            {
                if (records[i].required || records[i].userRequired)
                    selectedAttributes[cnt++] = records[i];
            }
        }
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp(file, frame) failed.", e, StsException.WARNING);
		}
	}
/*
	static public StsSegyVsp constructTraceAnalyzer(StsFile file, String stsDirectory, Frame frame)
	{
		try
		{
			return new StsSegyVsp(file, stsDirectory, frame, null);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.constructTraceAnalyzer(file, frame) failed.", e, StsException.WARNING);
			return null;
		}
	}

	static public StsSegyVsp constructTraceAnalyzer(String segyDirectory, String stsDirectory, String filename, Frame frame)
	{
		return constructTraceAnalyzer(StsFile.constructTraceAnalyzer(segyDirectory, filename), stsDirectory, frame);
	}
*/
	static public StsSegyVsp constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsSegyVsp(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.constructTraceAnalyzer(file, frame) failed.", e, StsException.WARNING);
			return null;
		}
	}

    public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages)
	{
		displayMessages = messages;
		return analyzeGrid();
    }

	public boolean analyzeGrid()
	{
		try
		{
			nTraces = segyData.getNTotalTraces();

            if(nTraces == 0) return false;
			if(nSlices == 0) return false;
			zMax = zMin + (nSlices - 1)*zInc;

			long ioff = segyData.fileHeaderSize;

			// Get first trace to set line and xLine number minimum
			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null) return false;

			// Get last trace to compute line and xLine number maximums
			lastTraceHeader = getTrace(nTraces-1);
			if(lastTraceHeader == null) return false;

			rowNumMin = 1;
			rowNumMax = 1;
			rowNumInc = 1;
			colNumMin = 1;
			colNumMax = nTraces;
			colNumInc = 1;

			nRows = 1;
			nCols = nTraces;
			angle = 0.0f;
			xInc = 1.0f;
			yInc = 1.0f;
			xOrigin = firstTraceHeader.x * horizontalScalar;
			yOrigin = firstTraceHeader.y * horizontalScalar;

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

	StsSEGYFormat.TraceHeader getTrace(int n)
	{
        return segyData.getTrace(n);
	}
/*
	public boolean computeScanHistogram()
	{
		scanHistogram = true;
		try
		{
			if(!analyzeTraces(0))
			{
				StsMessageFiles.errorMessage("StsSegyVsp:computeScanHistogram() failed, type appears to be IEEE");
				return false;
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsSegyVsp.computeScanHistogram failed to build Histogram");
			return false;
		}
	}
*/
	public void resetDataRange()
	{
        segyData.resetDataRange();
	}

	// JKF 03JULY2006

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
/*
	private boolean analyzeTraces()
	{
		boolean fmtOkay = true;
		interrupt = false;

		if (segyFormat != null)
		{
			resetDataRange();
		}

		int nTracesToScan = Math.max(((int)((float)nTraces * (scanPercentage / 100.0f))), 10); // Default to ~1% of traces

		// Set to run entire file if scan percentage is greater than 50%, will be faster
		if (scanPercentage > 50.0f && nTraces > 10)
			nTracesToScan = nTraces;

		initializeSegyIO();

		analyzeTraces = StsAnalyzeTraces.constructor(randomAccessSegyFile,
			fileHeaderSize,
			nTracesToScan,
			nTraces,
			nCroppedSlices,
			bytesPerSample,
			traceHeaderSize,
			getIsLittleEndian,
			sampleFormat,
			segyFormat.getTraceRecFromUserName("B2SCALCO"),
            segyFormat.getTraceRecFromUserName("TRC_TYPE"),
            getName());

		analyzeTraces.analyze();
		fmtOkay = true;
//		fmtOkay = analyzeTraces.getFmtOkay();
		if (fmtOkay)
		{

			dataMin = analyzeTraces.getDataMin();
			dataMax = analyzeTraces.getDataMax();
            traceTypes = analyzeTraces.getTraceTypes();

			dataAvg = (float)analyzeTraces.getAverage();
			if (scanHistogram)
			{
				dataCnt = analyzeTraces.getHistogramData();
				ttlHistogramSamples = analyzeTraces.getTtlHistogramSamples();
				calculateHistogram();
			}
		}
		return fmtOkay;
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

    public void setFilenames()
	{
	    setFilenames(getGroupname() + "." + stemname);
	}

    public void setFilenames(String stemname)
	{
        getCreateAttributeFilename();
        rowCubeFilename = getVolumeFilename(inline, stemname);
		createFloatRowVolumeFilename(stemname);
	}

    public int[] getTraceTypes()
    {
        return traceTypes;
    }

	public boolean readWriteVolume(StsProgressPanel panel, StsSEGYFormatRec[] attributeRecs, boolean multiComponent)
	{
        boolean success = false;
        this.panel = panel;
        this.multiComponent = multiComponent;

        if((multiComponent) && (traceTypes.length > 0))
        {
            int tempCols = nCols;
            for(int i=0; i<traceTypes.length; i++)
            {
                nCols = tempCols;
                actualTraces = 0;
                panel.setDescription("Procesing " + traceTypes[i] + " component\n");
                componentString = "_T" + Integer.toString(traceTypes[i]);
                traceId = traceTypes[i];
                success = readWriteVolume(panel, attributeRecs);
                panel.setDescription("     Found " + actualTraces + " traces of this component type out of " + nTraces + " traces.\n");
            }
        }
        else
            success = readWriteVolume(panel, attributeRecs);

        if(segyData.randomAccessSegyFile != null)
        {
            try { segyData.randomAccessSegyFile.close(); }
            catch (Exception e)  {}
        }
        componentString = "";
        return success;
	}

    public boolean readWriteVolume(StsProgressPanel panel, StsSEGYFormatRec[] attributeRecs)
    {
		long startTime = 0, stopTime = 0;
        this.panel = panel;
		nAttributes = attributeRecs.length;
		attributeRecords = attributeRecs;
        interrupt = false;

        setFilenames();
		initializeSegyIO();
        segyIO.initializeDataRange(dataMin, dataMax);
        try
		{
            segyData.randomAccessSegyFile.seek(0);
			long nTotalSamples = nTraces*nSlices;
			if(interrupt)
                return false;

            panel.setDescription("Reading SEGY VSP: " + segyFilename + "\n");
			deleteExistingFiles();

			if(rowCubeFilename != null)
			{
                StsFile file = StsFile.constructor(stsDirectory, rowCubeFilename);
                if (file == null) return false;
                rowFile = new StsBinaryFile(file);
                if(rowFile == null)return false;
                rowFile.openWrite();
            }

			if(attributeFilename != null)
			{
				attributeFile = getRandomAccessFile(stsDirectory + attributeFilename);
				if(attributeFile == null)return false;
			}
			if(rowFloatFilename != null)
			{
				StsFile file = StsFile.constructor(stsDirectory, rowFloatFilename);
				if (file == null) return false;
				floatRowFile = new StsBinaryFile(file);
                if (floatRowFile == null) return false;
                floatRowFile.openWrite();
			}

			if(runTimer)
			{
				startTime = System.currentTimeMillis();
//                totalTimer = new StsTimer();
//                totalTimer.start();
			}

			randomAccessSegyChannel = segyData.randomAccessSegyFile.getChannel();

			outputTraceFloats = new float[nSlices];
            outputTraceBytes = new byte[nSlices];
            traceAttributes = new double[nAttributes][nTraces];
			initializeAttributes();
            panel.setDescription("     Processing " + nTotalSamples + " samples.\n");

			inputPosition = (long)segyData.fileHeaderSize;
			outputPosition = 0;

//			scale = 254 / (dataMax - dataMin);
			segyTraceBytes = new byte[segyData.bytesPerTrace];
			traceBytes = new byte[nSlices];
			int nSamples = nTraces*nSlices;
			int bytesPerInputTrace = nSlices*segyData.bytesPerSample + segyData.traceHeaderSize;
			int nBytesPerOutputSample = getIsDataFloat() ? 4 : 1;

			clearMappedBuffers();

			if(runTimer) timer.start();

			inputBuffer = randomAccessSegyChannel.map(FileChannel.MapMode.READ_ONLY, inputPosition, nTraces*bytesPerInputTrace);
/*
			if(rowFloatFilename != null)
			{
				inlineFloatByteBuffer = inlineFloatChannel.map(FileChannel.MapMode.READ_WRITE, 0, nSamples*4);
				inlineFloatBuffer = inlineFloatByteBuffer.asFloatBuffer();
			}
*/
			if(runTimer) timer.stopPrint("   mapping of input and inline files: ");

            if(analyzeTraces != null)
                analyzeTraces.resetProgressBar();

			if(!processBlock(inputBuffer, rowFile, floatRowFile, outputTraceFloats, outputTraceBytes, nTraces))
//			if(!doProcessBlock(inputBuffer, inlineBuffer, inlineFloatBuffer, nTraces))
               return false;

			if(interrupt)
               return false;

			outputAttributes(actualTraces);

			if(runTimer)
			{
				stopTime = System.currentTimeMillis();
				long elapsedTime = (stopTime - startTime)/1000;
				System.out.println("Total time to process seismic data: " + elapsedTime + " secs");
//                totalTimer.stop("Total time to process seismic data: ");
			}
            getCreateHeaderFilename();
            panel.setDescription("     Writing volume header file: " + headerFilename + componentString + "\n");
			calculateHistogram(segyIO.histogramSamples, segyIO.nHistogramSamples);

            nCols = actualTraces;
            panel.appendLine("Completed processing of " + getName());

			if(!writeHeaderFile(componentString)) return false;
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVsp.readWriteVolume() failed.", e, StsException.WARNING);
			String message = new String("Failed to load volume " + stemname + ".\n" + "Error: " + e.getMessage());
            panel.appendLine("Processing of " + getName() + " failed.");
            new StsMessage(frame, StsMessage.WARNING, message);
			return false;
		}
		finally
		{
			if(memoryDebug) checkMemory("Before freeing memory.");
			clearMappedBuffers();
			close();
			if(memoryDebug) checkMemory("After freeing memory.");
		}
	}

	private void clearMappedBuffers()
	{
		if(runTimer) timer.start();

		if(inputBuffer != null)
		{
			inputBuffer.clear();
			inputBuffer = null;
		}
	/*
		if(inlineBuffer != null)
		{
			inlineBuffer.force();
			inlineBuffer.clear();
			inlineBuffer = null;
		}
	*/
	/*
		if(inlineFloatByteBuffer != null)
        {
            inlineFloatByteBuffer.force();
			inlineFloatByteBuffer.clear();
			inlineFloatByteBuffer = null;
		}
	*/

//        Runtime.getRuntime().gc();
		System.gc();
//        Runtime.getRuntime().runFinalization();

		if(runTimer) timer.stopPrint("Time to clear mapped buffers:");
	}

	private void checkMemory(String message)
	{
		long maxMemory = Runtime.getRuntime().maxMemory();
		long totalMemory = Runtime.getRuntime().totalMemory();
		long freeMemory = Runtime.getRuntime().freeMemory();
		long usedMemory = totalMemory - freeMemory;
		long availMemory = maxMemory - usedMemory;
		StsMessageFiles.infoMessage("SEGY Memory Debug for volume " + segyFilename + " " + message + " Max memory: " + maxMemory + " avail memory: " + availMemory);
	}

	public double getBinaryHeaderValue(StsSEGYFormatRec rec)
	{
        return segyData.getBinaryHeaderValue(rec);
	}

	private boolean processBlock(MappedByteBuffer inputBuffer, StsBinaryFile rowFile, StsBinaryFile floatRowFile, float[] outputTraceFloats, byte[] outputTraceBytes, int nBlockTraces)
	{
		int nTrace = -1;

		try
		{
			if(runTimer) timer.start();
			inputBuffer.position(0);
            panel.initialize(nBlockTraces);
            panel.setDescription("Processing " + getName());
            for(nTrace = 0; nTrace < nBlockTraces; nTrace++)
			{
                panel.appendLine("Processing trace " + nTrace + " of " + nBlockTraces);
				inputBuffer.get(segyTraceBytes);
				if(extractHeaderAttributes(segyTraceBytes, nTrace, actualTraces))
				{
					processTrace(segyTraceBytes, nTrace, outputTraceFloats);
					floatRowFile.setFloatValues(outputTraceFloats, false);
                    segyIO.processTrace(outputTraceFloats, nSlices, nTrace, outputTraceBytes);
                    rowFile.setByteValues(outputTraceBytes);
                    actualTraces++;
				}
                panel.setValue(nTrace+1);
			}
            panel.finished();
            panel.setDescription("Process completed.");
            if(runTimer) timer.stopPrint("process " + nTraces + " traces:");
			return true;
		}
		catch (Exception e)
		{
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            StsException.outputException("StsSegyVsp.processTraces() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

	final private boolean processTrace(byte[] segyTraceBytes, int nTrace, float[] outputSamples)
	{
		int n = 0, pos = 0, s = 0;

		try
		{
			segyIO.processTrace(segyTraceBytes, nTrace, outputSamples);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.processTrace() failed at  trace " + nTrace + " sample " + s + " pos " + pos, e, StsException.WARNING);
			return false;
		}
	}
/*
	private boolean doProcessBlock(MappedByteBuffer inputBuffer, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer, int nBlockTraces)
	{
		int nTrace = -1;

		try
		{
			if(runTimer) timer.start();
			inputBuffer.position(0);
			for(nTrace = 0; nTrace < nBlockTraces; nTrace++)
			{
				inputBuffer.get(segyTraceBytes);
				if(extractHeaderAttributes(segyTraceBytes, nTrace, actualTraces))
                {
                    processTrace(segyTraceBytes, inlineBuffer, nativeBuffer, nTrace);
                    actualTraces++;
                }
			}
			if(runTimer) timer.stopPrint("process " + nTraces + " traces:");

			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.processTraces() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

	final private boolean processTrace(byte[] segyTraceBytes, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer,int nTrace)
	{
		int n = 0, pos = 0, s = 0;

		try
		{
			int offset = traceHeaderSize;
			//Process Samples
			for(s = 0; s < nCroppedSlices; s++)
			{
				segyIO.processSampleBytes(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				offset += bytesPerSample;
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.processTrace() failed at  trace " + nTrace + " sample " + s + " pos " + pos, e, StsException.WARNING);
			return false;
		}
	}
*/
	private boolean extractHeaderAttributes(byte[] segyTraceBytes, int traceNum, int actualNum)
	{
        return segyData.extractHeaderAttributes(segyTraceBytes, traceNum, attributeRecords, traceAttributes);
	}
/*
	final private void processSampleBytes(byte[] segyTraceBytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		//		int sampleFormat = getSampleFormat();
		switch(sampleFormat)
		{
			case StsSEGYFormat.IBMFLT:
				sampleProcessIBMfloatSample(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			case StsSEGYFormat.IEEEFLT:
				sampleProcessIEEEfloatSample(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			case StsSEGYFormat.INT2:
				sampleProcessInt2ToBytes(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			case StsSEGYFormat.INT4:
				sampleProcessInt4ToBytes(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			case StsSEGYFormat.BYTE:
				sampleProcessByteToBytes(segyTraceBytes, offset, inlineBuffer);
				break;
			case StsSEGYFormat.FLOAT8:
				sampleProcessFloatByteToBytes(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			case StsSEGYFormat.FLOAT16:
				sampleProcessFloatInt2ToBytes(segyTraceBytes, offset, inlineBuffer, nativeBuffer);
				break;
			default:
				new StsMessage(null, StsMessage.WARNING, "Sample processing failed: unsupported sample format " + sampleFormat);
				break;
		}
	}

	final private void sampleProcessIBMfloatSample(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		float sampleFloat =  StsMath.convertIBMFloatBytes(bytes, offset, getIsLittleEndian);
		if(nativeBuffer != null) nativeBuffer.put(sampleFloat);
		checkSetFloatToByte(inlineBuffer, sampleFloat);
	}

	final private void checkSetFloatToByte(MappedByteBuffer inlineBuffer, float sampleFloat)
	{
		if(((sampleFloat > dataMax) || (sampleFloat < dataMin)) && (clipToNull))
			inlineBuffer.put(nullByte);
		else
		{
			accumulateHistogram(sampleFloat);
			inlineBuffer.put(StsMath.unsignedIntToSignedByte254(Math.round((sampleFloat - dataMin) * scale)));
		}
	}

	final private void sampleProcessIEEEfloatSample(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		float sampleFloat = Float.intBitsToFloat(StsMath.convertIntBytes(bytes, offset, getIsLittleEndian));
		if(nativeBuffer != null)
			nativeBuffer.put(sampleFloat);
		checkSetFloatToByte(inlineBuffer, sampleFloat);
	}

	final private void sampleProcessInt2ToBytes(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		int sampleInt = StsMath.convertInt2Bytes(bytes, offset, getIsLittleEndian);
		if(nativeBuffer != null)
			nativeBuffer.put((float)sampleInt);
		checkSetFloatToByte(inlineBuffer, (float)sampleInt);
	}

	final private void sampleProcessInt4ToBytes(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		System.out.println("Do we need to add support for a 4 byte integer?????? Why?");
	}

	final private void sampleProcessByteToBytes(byte[] bytes, int offset, MappedByteBuffer inlineBuffer)
	{
		checkSetFloatToByte(inlineBuffer, (float)bytes[offset]);
	}

	final private void sampleProcessFloatInt2ToBytes(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		int scaledInt = (int)StsMath.convertInt2Bytes(bytes, offset, getIsLittleEndian);
		int dataScale = this.segyFormat.getDataScale(bytes, getIsLittleEndian);
		float scaledFloat = (float)scaledInt / (float)(Math.pow(2.0f, dataScale));
		if(nativeBuffer != null)
			nativeBuffer.put(scaledFloat);
		checkSetFloatToByte(inlineBuffer, scaledFloat);
	}

	final private void sampleProcessFloatByteToBytes(byte[] bytes, int offset, MappedByteBuffer inlineBuffer, FloatBuffer nativeBuffer)
	{
		int dataScale = this.segyFormat.getDataScale(bytes, getIsLittleEndian);
		float scaledFloat = (float)bytes[offset] / (float)(Math.pow(2.0f, dataScale));
		if(nativeBuffer != null)
			nativeBuffer.put(scaledFloat);
		checkSetFloatToByte(inlineBuffer, scaledFloat);
	}
*/
	/** converts an unsigned int to a signedByte value between/including 0 to 254 */
/* moved to StsMath for more general use.
	final public byte unsignedIntToSignedByte254(int i)
	{
		if(i >= 255) i = 254;
		if(i < 0) i = 0;
		dataCnt[i] += 1;
		ttlHistogramSamples++;
		return (byte)i;
	}
*/
	public boolean hasOutputFiles()
	{
		setFilenames();
        File file;
        file = new File(stsDirectory + attributeFilename);
		if(file.exists()) return true;
		file = new File(stsDirectory + rowCubeFilename);
		if (file.exists()) return true;
		return false;
	}

	public void cancel() { interrupt = true; }

	public void close()
	{
    /*
        try
		{
			memoryAllocation.freeMemory();
		}
		catch(Exception e) {}
     */
		if(rowFile != null)       try { rowFile.close();       } catch(Exception e) {}
		if(attributeFile != null) try { attributeFile.close(); } catch(Exception e) {}
		if(floatRowFile != null)  try { floatRowFile.close();  } catch(Exception e) {}
	}

	private RandomAccessFile getRandomAccessFile(String pathname)
	{
		try
		{
			File file = new File(pathname);
		/*
			if(file.exists())
			{
				boolean deleteOK = file.delete();
				if(!deleteOK)
				{
					new StsMessage(this.frame, StsMessage.WARNING, "Failed to delete file " + pathname + ".\n" +
						"Exit wizard and check if another application has it open.");
					return null;
				}
			}
		*/
			return new RandomAccessFile(pathname, "rw");
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVsp.getRandomAccessFile() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

	private float testCubeConversion(float sampleFloat)
	{
		if(sampleFloat == 0.0f) return 0.0f;
		if(sampleFloat < 100.0f) return 1.0f;
		else if(sampleFloat < 1000.0f) return 2.0f;
		else return 3.0f;
	}

	// if crossline+ direction is 90 degrees CCW from inline+, this is isXLineCCW; otherwise not
	// angle is from X+ direction to inline+ direction (0 to 360 degrees)
	private void initializeAngle()
	{
		angle = angle % 360.0f;
		if (angle < 0.0f) angle += 360.0f;
		setAngle();
	}

/*    private void accumulateHistogram(float bindex)
	{
		addToHistogram(bindex);
	}

	private void accumulateHistogram(short bindex)
	{
		addToHistogram((float)bindex);
	}

	private void accumulateHistogram(int bindex)
	{
		System.out.println("Do we need to add support for a 4 byte integer?????? Why?");
	}

	private void accumulateHistogram(byte bindex)
	{
		addToHistogram((float)bindex);
	}
*/
/*
	private void accumulateHistogram(float bindex)
	{
		byte bsamp = 0;
		float scale = 254 / (dataMax - dataMin);
		float scaledFloat = (bindex - dataMin)*scale;
		int scaledInt = Math.round(scaledFloat);
		bsamp = StsMath.unsignedIntToSignedByte254(scaledInt);

		int index = StsMath.signedByteToUnsignedInt(bsamp);
		if(index > 254) index = 254;
		if(index < 0) index = 0;
		dataCnt[index] += 1;
		ttlHistogramSamples++;
	}
*/

	public String getLabel() { return stemname; }

	public void setIsMinMemory(boolean isMinMemory) {this.isMinMemory = isMinMemory;}
	public boolean isMinMemory(){return isMinMemory;}

	public int getNTraces() { return segyData.nTotalTraces; };

	public float getTInc()
	{
		return zInc;
	}
	public static void main(String[] args)
	{
		String segyDirectory = "c:/CoreLabsData";
		String stsDirectory = "c:/CoreLabsData";
		String name = "small_seis517-540.sgy";
		try
		{
			StsFile file = StsFile.constructor(segyDirectory, name);
			StsModel model = new StsModel();
			StsSEGYFormat segyFormat = StsSEGYFormat.constructor(model, StsSEGYFormat.VSP);
			StsSegyVsp segyVolume = StsSegyVsp.constructor(null, file, stsDirectory, null, segyFormat);
			segyVolume.readWriteVolume(null,null);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVsp.main() failed.",
				e, StsException.WARNING);
		}
	}

   public void initSelectedAttributes(int count)
   {
       if(selectedAttributes != null)
           selectedAttributes = null;
       selectedAttributes = new StsSEGYFormatRec[count];
   }
   public void setSelectedAttribute(int index, StsSEGYFormatRec record)
   {
       selectedAttributes[index] = record;
   }
   public StsSEGYFormatRec[] getSelectedAttributes() { return selectedAttributes; }

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

	void initializeAttributes()
	{
		int nAttributes = attributeRecords.length;
		attributeNames = new String[nAttributes];
		for (int n = 0; n < nAttributes; n++)
		{
			StsSEGYFormatRec record = attributeRecords[n];
			String userName = record.getUserName();
			attributeNames[n] = userName;
		}
	}

	private void outputAttributes(int numTraces)
	{
		StsMappedDoubleBuffer attributeBuffer = null;
		try
		{
            getCreateAttributeFilename();
            attributeBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "rw");
			if (attributeBuffer == null)
               return;
			attributeBuffer.map(0, numTraces*nAttributes);
			for (int j = 0; j < attributeRecords.length; j++)
               attributeBuffer.put(traceAttributes[j], 0, numTraces);
			attributeBuffer.clean();
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVsp.outputBlock.process() failed.", e, StsException.WARNING);
			return;
		}
		finally
		{
			if(attributeBuffer != null) attributeBuffer.close();
		}
	}

	public void setVerticalUnits(byte vUnits)
	{
		StsParameters.getDistanceUnitString(vUnits);
	}

    public boolean isAValidFileSize()
    {
        return segyData.isAValidFileSize();
    }

	public void setZDomain(byte zDomain)
	{
		this.zDomain = StsParameters.TD_ALL_STRINGS[zDomain];
	}

    public String getTypeAsString()
    {
        if(segyData.segyFormat != null)
            return StsVspClass.vspTypeStrings[segyData.segyFormat.getType()];
        else
            return StsVspClass.vspTypeStrings[0];
    }


}
