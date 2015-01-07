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

/** reads a SEGY file in blocks into memory, converts samples to scaled bytes and outputs
 *  them into: seis2d.txt.name, seis2d.attributes.name, seis2d.bytes.name, seis2d.floats.name
 *  Input sample data can be in any format supported by SEGY.
 */
public class StsSegyLine2d extends StsSeismicLine2d implements Serializable
{
	/** Is the clipped data to be set to null or the maximum and minimum */
//	transient protected boolean clipToNull = false;
	/** Total samples in each of 255 steps */
//	transient int dataCnt[] = new int[255];
//	transient int ttlHistogramSamples = 0;
	transient int traceId = -1;
//	transient protected StsSEGYFormat segyFormat = null;
	transient int actualTraces = 0;

//	transient byte[] textHeader;
//	transient byte[] binaryHeader;
//	transient byte[] traceHeader;

	/** random access file for reading input segy file */
//	transient RandomAccessFile randomAccessSegyFile;
	transient FileChannel randomAccessSegyChannel;
	transient MappedByteBuffer inputBuffer = null;
//	transient int fileHeaderSize = 0;
//	transient int bytesPerTrace;
//	transient long segyFileSize;
	transient int nTotalTraces = 0;
	transient long inputPosition = 0;

	transient boolean scanHistogram = true;
    /** data scaling: 254/(dataMax - dataMin) */
	transient float scale = 1.0f;
    /** data scaling constant: -dataMin*scale */
    transient float scaleConstant = 0.0f;
    /** useful traces in cube */
	transient StsSEGYFormat.TraceHeader firstTraceHeader, lastTraceHeader;
	transient StsSEGYFormat.TraceHeader endFirstRowTraceHeader, endFirstColTraceHeader;

	transient boolean displayMessages = true;

	/** wizard frame used for locating displayed dialogs */
	transient Frame frame;

    transient StsProgressPanel progressPanel;

	/** timer used in testing performance */
	transient StsTimer timer = null;
	transient StsTimer totalTimer = null;
	/** set to true if timing tests are desired */
	transient boolean runTimer = false;

	transient static boolean isMinMemory = true;

	transient StsBinaryFile rowFile = null;
	transient StsBinaryFile floatRowFile = null;
	transient StsBinaryFile attributeFile = null;

	transient StsSEGYFormatRec[] selectedAttributes = null;

	transient boolean interrupt = false;

//    transient boolean progressBar = false;
	transient int progressInc = 0;
	transient StsProgressPanel pd = null;
	transient JDialog pdialog = null;

	transient static final int INLINE = 0;
	transient static final int XLINE = 1;
	transient static final int TD = 2;

    /** byte array for input trace */
    transient byte[] segyTraceBytes;
	/** temp float array used for output trace */
	transient float[] outputFloats;
    /** temp byte array used for output trace */
    transient byte[] outputBytes;
    /** attributes for each trace; sequential for each attribute [nAttributes][nTraces] */
	transient double[][] traceAttributes;
	/** number of traces written to attribute array */
	transient int nAttributeTrace = 0;
	/** Attribute records for extraction */
	transient StsSEGYFormatRec[] attributeRecords = null;

    transient int cdpxAttributeIndex = -1;
    transient int cdpyAttributeIndex = -1;
    transient int cdpAttributeIndex = -1;

	static public float defaultScanPercentage = 10.0f;

	public final static byte nullByte = StsParameters.nullByte;

    transient public StsObject assocLine = null;

	/** acceptable suffixes for segy files (upper or lower case ok) */
	static String[] segySuffixes = new String[] {".sgy", ".segy"};

    static final boolean debug = false;
	static final boolean memoryDebug = false;

	static final float nullValue = StsParameters.nullValue;

	public StsSegyLine2d()
	{}

	public StsSegyLine2d(boolean persist)
	{
		super(persist);
	}
	/**
	 *
	 * @param file segy input file
	 * @param stsDirectory directory where S2S seismic byte volumes are written
	 * @param frame frame used for centering this dialog
	 * @param segyFormat format information for reading segy file
	 */
	StsSegyLine2d(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws StsException
	{
		super(false);

		try
		{
			if(segyFormat == null)
			{
				throw new StsException("StsSegyLine2d.constructor failed: segyFormat cannot be null.");
			}
			if(runTimer)timer = new StsTimer();
            segyData = StsSegyData.constructor(seismicWizard, file, this, segyFormat);
            if(segyData == null)
            {
                throw new StsException("StsSegyLine2d.constructor failed: couldn't construct StsSegyData.");
            }
			segyDirectory = file.getDirectory();
			segyFilename = file.getFilename();
			File sfile = new File(segyDirectory + segyFilename);
			segyLastModified = sfile.lastModified();
			this.stsDirectory = stsDirectory;
			this.frame = frame;
            this.zDomain = segyFormat.getZDomainString();
//            segyData.randomAccessSegyFile = new RandomAccessFile(file.getPathname(), "r");
            
            zMin = 0.0f;

//            readFileHeader();
//			segyData.analyzeBinaryHdr();

			stemname = StsStringUtils.trimSuffix(segyFilename, segySuffixes);
			setName(stemname);

			StsSEGYFormatRec[] records = segyFormat.getAllTraceRecords();
			selectedAttributes = new StsSEGYFormatRec[segyFormat.getNumberOfRequiredRecords()];
			int cnt = 0;
			for(int i = 0; i < records.length; i++)
			{
				if(records[i].required || records[i].userRequired)
					selectedAttributes[cnt++] = records[i];
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d(file, frame) failed.", e, StsException.WARNING);
		}
	}

	static public StsSegyLine2d constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsSegyLine2d(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d.constructTraceAnalyzer(file, frame) failed.", e, StsException.WARNING);
			return null;
		}
	}

	public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages)
	{
		displayMessages = messages;
        this.progressPanel = progressPanel;
        return analyzeGeometry();
	}

	public boolean analyzeGeometry()
	{
		try
		{
        /*
            traceHeaderSize = segyFormat.traceHeaderSize;
			traceHeader = new byte[traceHeaderSize];
			bytesPerTrace = nCroppedSlices * bytesPerSample + traceHeaderSize;
			segyFileSize = randomAccessSegyFile.length();
	    */
            nTotalTraces = segyData.getNTotalTraces();
//			nTraces = (int)((segyFileSize - fileHeaderSize) / bytesPerTrace);
			if(nTotalTraces == 0)return false;
			if(nSlices == 0)return false;
			zMax = zMin + (nSlices - 1) * zInc;

			long ioff = segyData.getFileHeaderSize();

			// Get first trace to set line and xLine number minimum
			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null)return false;

			// Get last trace to compute line and xLine number maximums
			lastTraceHeader = getTrace(nTotalTraces - 1);
			if(lastTraceHeader == null)return false;

			rowNumMin = 1;
			rowNumMax = 1;
			rowNumInc = 1;
			colNumMin = 1;
			colNumMax = nTotalTraces;
			colNumInc = 1;

			nRows = 1;
			nCols = nTotalTraces;
			angle = 0.0f;
			xInc = 1.0f;
			yInc = 1.0f;
			xOrigin = firstTraceHeader.x * horizontalScalar;
			yOrigin = firstTraceHeader.y * horizontalScalar;

			return true;
		}
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
				StsMessageFiles.errorMessage("StsSegyLine2d:computeScanHistogram() failed, type appears to be IEEE");
				return false;
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsSegyLine2d.computeScanHistogram failed to build Histogram");
			return false;
		}
	}
*/
	public void resetDataRange()
	{
        segyData.resetDataRange();
	}

    public boolean isAValidFileSize()
    {
        return segyData.isAValidFileSize();
    }

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

    public int getTotalTraces()
	{
		return segyData.nTotalTraces;
	}

    public void setFilenames()
	{
		setIsDataFloat(true);
        rowFloatFilename = createFloatDataFilename();
        rowCubeFilename = createByteDataFilename();
        attributeFilename = getCreateAttributeFilename();
    }

	public boolean readWriteLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecs)
	{
        this.progressPanel = progressPanel;
        boolean success = false;

		success = readWrite(attributeRecs);

		if(segyData.randomAccessSegyFile != null)
		{
			try
			{
				segyData.randomAccessSegyFile.close();
			}
			catch(Exception e)
			{}
		}
		return success;
	}

    public boolean readWrite(StsSEGYFormatRec[] attributeRecs)
	{
		long startTime = 0, stopTime = 0;
        initializeSegyIO();
        segyIO.initializeDataRange(dataMin, dataMax);

		nAttributes = attributeRecs.length;
		attributeRecords = attributeRecs;
		interrupt = false;
//		cdpX = new double[65000];
//		cdpY = new double[65000];
//		cdp = new int[65000];

		setFilenames();
		try
		{
			segyData.randomAccessSegyFile.seek(0);
			long nTotalSamples = nTotalTraces * nSlices;
			if(interrupt)
				return false;

			progressPanel.appendLine("Reading 2D SEGY : " + segyFilename);
			deleteExistingFiles();
			/*
			   if(rowCubeFilename != null)
			   {
			 rowFile = getRandomAccessFile(stsDirectory + rowCubeFilename);
			 if(rowFile == null)return false;
//				inlineChannel = rowFile.getChannel();
			   }
			 */
            if(rowFloatFilename != null)
			{
                String floatFilename = createFloatDataFilename();
                StsFile file = StsFile.constructor(stsDirectory, rowFloatFilename);
				if(file == null)return false;
				floatRowFile = new StsBinaryFile(file);
				floatRowFile.openWrite();
				if(floatRowFile == null)return false;
			}
			if(rowCubeFilename != null)
			{
				StsFile file = StsFile.constructor(stsDirectory, rowCubeFilename);
				if(file == null)return false;
				rowFile = new StsBinaryFile(file);
				rowFile.openWrite();
				if(rowFile == null)return false;
			}
			if(runTimer)
			{
				startTime = System.currentTimeMillis();
//                totalTimer = new StsTimer();
//                totalTimer.start();
			}

			randomAccessSegyChannel = segyData.randomAccessSegyFile.getChannel();

			outputFloats = new float[nSlices];
            outputBytes = new byte[nSlices];
            traceAttributes = new double[nAttributes][nTotalTraces];
			if(!initializeAttributes())return false;
			progressPanel.appendLine("     Processing " + nTotalSamples + " samples.");

            inputPosition = (long)segyData.fileHeaderSize;

			scale = 254 / (dataMax - dataMin);
            scaleConstant = -dataMin*scale;
            segyTraceBytes = new byte[segyData.bytesPerTrace];
			int bytesPerInputTrace = nSlices * segyData.bytesPerSample + segyData.traceHeaderSize;

			clearMappedBuffers();

			if(runTimer)timer.start();

			inputBuffer = randomAccessSegyChannel.map(FileChannel.MapMode.READ_ONLY, inputPosition, nTotalTraces * bytesPerInputTrace);

//			System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " cropBox.nBlockSamples: " + nSamples);
//			inlineBuffer = inlineChannel.map(FileChannel.MapMode.READ_WRITE, 0, nSamples);


			if(runTimer)timer.stopPrint("   mapping of input  files: ");


			if(!processBlock(inputBuffer, floatRowFile, rowFile, outputFloats, nTotalTraces))
				return false;

            if(interrupt)
				return false;

//			calculateHistogram(segyIO.histogramSamples, segyIO.nHistogramSamples);
			outputAttributes(actualTraces);

			if(runTimer)
			{
				stopTime = System.currentTimeMillis();
				long elapsedTime = (stopTime - startTime) / 1000;
				System.out.println("Total time to process seismic data: " + elapsedTime + " secs");
//                totalTimer.stop("Total time to process seismic data: ");
			}
            getCreateHeaderFilename();
            progressPanel.appendLine("     Writing 2d line header file: " + headerFilename);

			nCols = actualTraces;
            adjustXYRange();

			if(!writeHeaderFile()) return false;
			progressPanel.appendLine("\tLine Header File Written");
            progressPanel.setDescription("Completed line " + getName());
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d.readWriteLines() failed.", e, StsException.WARNING);
			String message = new String("Failed to load  " + stemname + ".\n" + "Error: " + e.getMessage());
			new StsMessage(frame, StsMessage.WARNING, message);
			return false;
		}
		finally
		{
			if(memoryDebug)checkMemory("Before freeing memory.");
			clearMappedBuffers();
			close();
			if(memoryDebug)checkMemory("After freeing memory.");
		}
	}

    public boolean writeHeaderFile(String extension)
	{
		try
		{
			String origStemname = stemname;
			setName(stemname + extension);
			writeHeaderFile();
			stemname = origStemname;
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException(StsToolkit.getSimpleClassname(this) + ".writeHeaderFile() failed.", e, StsException.WARNING);
			return false;
		}
	}

    public boolean writeHeaderFile()
	{
		try
		{
            getCreateHeaderFilename();
            StsParameterFile.writeObjectFields(stsDirectory + headerFilename, this, StsSegyLine2d.class, StsMainObject.class);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException(StsToolkit.getSimpleClassname(this) + ".writeHeaderFile() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private void clearMappedBuffers()
	{
		if(runTimer)timer.start();

		if(inputBuffer != null)
		{
			inputBuffer.clear();
			inputBuffer = null;
		}

//        Runtime.getRuntime().gc();
		System.gc();
//        Runtime.getRuntime().runFinalization();

		if(runTimer)timer.stopPrint("Time to clear mapped buffers:");
	}

	private void checkMemory(String message)
	{
		long maxMemory = Runtime.getRuntime().maxMemory();
		long totalMemory = Runtime.getRuntime().totalMemory();
		long freeMemory = Runtime.getRuntime().freeMemory();
		long usedMemory = totalMemory - freeMemory;
		long availMemory = maxMemory - usedMemory;
		StsMessageFiles.infoMessage("SEGY Memory Debug for line " + segyFilename + " " + message + " Max memory: " + maxMemory + " avail memory: " + availMemory);
	}

	private void attributeMissingMessage(String name)
	{
		new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + name + " attribute");
	}

	public double getBinaryHeaderValue(StsSEGYFormatRec rec)
	{
		return segyData.getBinaryHeaderValue(rec);
	}

    static public int getAssocCdpIndex( int[] cdpArray, int searchValue)
    {
        for( int i=0; i<cdpArray.length; i++)
        {
            if( StsMath.sameAs(cdpArray[i], searchValue))
                return i;
        }
        return -1;
    }

    private void getRotatedCdpXY( int[] rotatedCdp, float[] rotatedCdpX, float[] rotatedCdpY)
    {
        StsSeismicLine seismicLine = (StsSeismicLine)assocLine;
        System.arraycopy( seismicLine.cdp,  0, rotatedCdp,  0, nTotalTraces);
        System.arraycopy( seismicLine.cdpX, 0, rotatedCdpX, 0, nTotalTraces);
        System.arraycopy( seismicLine.cdpY, 0, rotatedCdpY, 0, nTotalTraces);
    }

    private boolean initializeCdpXY()
    {
        StsSeismicBoundingBox seismicBoundingBox = (StsSeismicBoundingBox)assocLine;
        float cdpInc = cdpInterval;
        if( cdpInc <= 0) return false;

        float bbXmin = seismicBoundingBox.getXMin();
        float bbXmax = seismicBoundingBox.getXMax();
        float bbYmin = seismicBoundingBox.getYMin();
        float bbYmax = seismicBoundingBox.getYMax();

        float dx = bbXmax - bbXmin;
        float dy = bbYmax - bbYmin;
        float length = (float)Math.sqrt(dx*dx + dy*dy);

        nTotalTraces = (int)(length/cdpInc) + 1;
        if( nTotalTraces < 1) return false;

        cdp  = new int[nTotalTraces];
        cdpX = new float[nTotalTraces];
        cdpY = new float[nTotalTraces];

        float xinc = 0.0f;
        float yinc = 0.0f;
        if( nTotalTraces > 1)
        {
            xinc = dx/(nTotalTraces - 1);
            yinc = dy/(nTotalTraces - 1);
        }
        for( int i=0; i<nTotalTraces; i++)
        {
            cdp[i] = i + 1;
            cdpX[i] = (float)i*xinc;
            cdpY[i] = (float)i*yinc;
        }

        setAssocLine( this);
        setAssocLineName( "User-defined line");
        return true;
    }

    private int checkAssocCdp( int assocCdpIndex, StsSeismicLine seismicLine, int traceNo, boolean force, int previous)
    {
        int cdpIndex = assocCdpIndex;

        if( cdpIndex >= 0) return cdpIndex;
        if( previous < 0) return cdpIndex;

        if( force)
            cdpIndex = previous;
        else
        {
            if( StsYesNoDialog.questionValue(this.frame, "Cdp "+ seismicLine.cdp[traceNo] + " missing for line " + seismicLine.getName() +
                                         ".\n\n Yes - Use previous Cdp for XY \n No - Cancel Line Load"))
                cdpIndex = previous;
        }
        return cdpIndex;
    }

    private boolean getUnrotatedCdpXY( int[] rotatedCdp, float[] rotatedCdpX, float[] rotatedCdpY,
                                       double[] unRotatedCdp, double[] unRotatedCdpX, double[] unRotatedCdpY,
                                       boolean force)
    {
        int previous = -1;
        StsSeismicLine seismicLine = (StsSeismicLine)assocLine;
        for( int i=0; i<nTotalTraces; i++)
        {
            int assocCdpIndex = getAssocCdpIndex( rotatedCdp, seismicLine.cdp[i]);
            assocCdpIndex = checkAssocCdp( assocCdpIndex, seismicLine, i, force, previous);
            if( assocCdpIndex < 0 ) return false;

            float[] unrotatedXY = currentModel.getProject().getUnrotatedRelativeXYFromRotatedXY(rotatedCdpX[assocCdpIndex], rotatedCdpY[assocCdpIndex]);
            unRotatedCdp[i]  = (double)seismicLine.cdp[i];
            unRotatedCdpX[i] = seismicLine.xOrigin + (double)unrotatedXY[0];
            unRotatedCdpY[i] = seismicLine.yOrigin + (double)unrotatedXY[1];
            previous =  assocCdpIndex;
        }
        return true;
    }

    public boolean checkGeomAttributes()
    {
        if( assocLine == null) return false;
        
        if( assocLine.getName() == "User-defined line")
            return initializeCdpXY();
        else
            return overrideGeomAttributes( true, false);
    }

    private boolean overrideGeomAttributes( boolean check, boolean force)
    {
        if( !(assocLine instanceof StsSeismicLine)) return false;
        
        int[]   rotatedCdp  = new int[nTotalTraces];
        float[] rotatedCdpX = new float[nTotalTraces];
        float[] rotatedCdpY = new float[nTotalTraces];
        getRotatedCdpXY( rotatedCdp, rotatedCdpX, rotatedCdpY);

        double[] unRotatedCdp  = new double[nTotalTraces];
        double[] unRotatedCdpX = new double[nTotalTraces];
        double[] unRotatedCdpY = new double[nTotalTraces];
        if( !getUnrotatedCdpXY( rotatedCdp, rotatedCdpX, rotatedCdpY, unRotatedCdp, unRotatedCdpX, unRotatedCdpY, force))
            return false;

        if( check) return true;

        System.arraycopy( unRotatedCdp,  0, traceAttributes[cdpAttributeIndex],  0, nTotalTraces);
        System.arraycopy( unRotatedCdpX, 0, traceAttributes[cdpxAttributeIndex], 0, nTotalTraces);
        System.arraycopy( unRotatedCdpY, 0, traceAttributes[cdpyAttributeIndex], 0, nTotalTraces);
        return true;
    }

	private boolean processBlock(MappedByteBuffer inputBuffer, StsBinaryFile floatRowFile, StsBinaryFile rowFile,
                                 float[] outputFloats, int nBlockTraces)
	{
		int nTrace = -1;

		try
		{
			if(runTimer)timer.start();
			inputBuffer.position(0);
        /*
            int bytesPerInputTrace = nCroppedSlices * bytesPerSample + traceHeaderSize;
            byte[] byteDump = new byte[nBlockTraces*bytesPerInputTrace];
            inputBuffer.get(byteDump);
            int nFloatValues = nBlockTraces*bytesPerInputTrace/4;
            float[] floats = new float[nFloatValues];
            StsMath.convertIBMFloatBytes(byteDump, floats, nFloatValues, getIsLittleEndian);
            inputBuffer.position(0);
        */
            progressPanel.setDescriptionAndLevel("Processing line " + getName(), StsProgressBar.INFO);
            progressPanel.setMaximum(nBlockTraces);
            int tracesPerIncrement = nBlockTraces/100;
            int incrementCounter = 0;
			for(nTrace = 0; nTrace < nBlockTraces; nTrace++)
			{
				inputBuffer.get(segyTraceBytes);
				if(extractHeaderAttributes(segyTraceBytes, nTrace, actualTraces))
				{
					processTrace(segyTraceBytes, nTrace, outputFloats);
                    normalizeAndScaleFloatsToBytes(outputFloats, outputBytes);
                    for (int kk=0; kk < outputFloats.length; kk++)
					   segyIO.accumulateHistogramSamples(outputFloats[kk]);
					floatRowFile.setFloatValues(outputFloats, false);
                    rowFile.setByteValues(outputBytes, false);
                    actualTraces++;
				}
                incrementCounter++;
                if(incrementCounter == tracesPerIncrement)
                {
                    progressPanel.setValue(nTrace);
                    incrementCounter = 0;
                }
            }

            if( assocLine != null)
            {
                if( !overrideGeomAttributes(false, true))
                    return false;
            }

            progressPanel.finished();
            if(runTimer)timer.stopPrint("process " + nTotalTraces + " traces:");
			return true;
		}
		catch(Exception e)
		{
            progressPanel.finished();
            progressPanel.setDescriptionAndLevel("Exception thrown: " + e.getMessage(), StsProgressBar.ERROR);
            StsException.outputException("StsSegyLine2d.processTraces() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

    private void adjustXYRange()
    {
        double[] cdpX = traceAttributes[cdpxAttributeIndex];
        double[] cdpY = traceAttributes[cdpyAttributeIndex];
        double[] xRange = StsMath.arrayMinMax(cdpX);
        double[] yRange = StsMath.arrayMinMax(cdpY);

        xOrigin = xRange[0];
        yOrigin = yRange[0];
        xMin = 0.0f;
        yMin = 0.0f;
        xMax = (float)(xRange[1] - xRange[0]);
        yMax = (float)(yRange[1] - yRange[0]);

        // leave incs zero; otherwise project gets confused on how to set them reasonably
//        xInc = (xMax - xMin)/(nCols - 1); // approx inc as spacing may vary
//        yInc = (yMax - yMin)/(nCols - 1); // approx inc as spacing may vary

        rowNumMin = 0;
		colNumMin = 0;
		rowNumMax = 0;
		colNumMax = actualTraces - 1;
		rowNumInc = 1;
		colNumInc = 1;
    }

    private void scaleTrace(float[] outputFloats, byte[] outputBytes)
    {
        for(int n = 0; n < nSlices; n++)
            outputBytes[n] = StsMath.unsignedIntToUnsignedByte((int)(outputFloats[n]*scale + scaleConstant));
    }

    final private boolean processTrace(byte[] segyTraceBytes, int nTrace, float[] outputSamples)
	{
		int n = 0, pos = 0, s = 0;

		try
		{
			segyIO.processTrace(segyTraceBytes, nTrace, outputSamples);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d.processTrace() failed at  trace " + nTrace + " sample " + s + " pos " + pos, e, StsException.WARNING);
			return false;
		}
	}

	private boolean extractHeaderAttributes(byte[] segyTraceBytes, int traceNum, int actualNum)
	{
        return segyData.extractHeaderAttributes(segyTraceBytes, actualNum, attributeRecords, traceAttributes);

	}


	public boolean overwriteFiles()
	{
		setFilenames();
//		String fullStemname = volumeFilePrefix+stemname+".nativePlanes";
		//System.out.println("check for exist "+fullStemname);
		File fi = getHeaderFile();
		return fi.exists();
	}

	public void cancel()
	{
		interrupt = true;
	}

	public void close()
	{
		if(attributeFile != null)try
		{
			attributeFile.close();
		}
		catch(Exception e)
		{}
		if(floatRowFile != null)
		{
			try
			{
				floatRowFile.close();
			}
			catch(Exception e)
			{}
		}
        if(this.rowFile != null)
        {
			try
			{
				rowFile.close();
			}
			catch(Exception e)
			{}
        }
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
			StsException.outputException("StsSegyLine2d.getRandomAccessFile() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	private float testCubeConversion(float sampleFloat)
	{
		if(sampleFloat == 0.0f)return 0.0f;
		if(sampleFloat < 100.0f)return 1.0f;
		else if(sampleFloat < 1000.0f)return 2.0f;
		else return 3.0f;
	}

	// if crossline+ direction is 90 degrees CCW from inline+, this is isXLineCCW; otherwise not
	// angle is from X+ direction to inline+ direction (0 to 360 degrees)
	private void initializeAngle()
	{
		angle = angle % 360.0f;
		if(angle < 0.0f)angle += 360.0f;
		setAngle();
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

    public int getNTraces()
	{
		return nTotalTraces;
	}

	public float getTInc()
	{
		return zInc;
	}


    public boolean hasOutputFiles()
    {
        File file = getHeaderFile(stemname);
        return file.exists();
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
			StsSEGYFormat segyFormat = StsSEGYFormat.constructor(model, StsSEGYFormat.POSTSTACK2D);
			StsSegyLine2d segyVolume = StsSegyLine2d.constructor(null, file, stsDirectory, null, segyFormat);
//			segyVolume.readWriteLines(null, null);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d.main() failed.",
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

	public StsSEGYFormatRec[] getSelectedAttributes()
	{
		return selectedAttributes;
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

	boolean initializeAttributes()
	{
		int nAttributes = attributeRecords.length;
		attributeNames = new String[nAttributes];
		for(int n = 0; n < nAttributes; n++)
		{
			StsSEGYFormatRec record = attributeRecords[n];
			String userName = record.getUserName();
			attributeNames[n] = userName;
			if(userName.equals(StsSEGYFormat.CDP))
				cdpAttributeIndex = n;
//		if(userName.equals("ILINE_NO"))
//			ilineAttributeIndex = n;
//		else if(userName.equals("XLINE_NO"))
//			xlineAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_X))
				cdpxAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_Y))
				cdpyAttributeIndex = n;
//		else if(userName.equals("OFFSET"))
//			offsetAttributeIndex = n;
		}
		boolean attributesOk = true;

		if(cdpAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace CDP number");
		}
		if(cdpxAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace CDP X value");
		}
		if(cdpyAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace CDP Y value");
		}

		return attributesOk;

	}

	private void outputAttributes(int numTraces)
	{
		StsMappedDoubleBuffer attributeBuffer = null;
		try
		{
            getCreateAttributeFilename();
            attributeBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "rw");
			if(attributeBuffer == null)
				return;
			if(!attributeBuffer.map(0, numTraces * nAttributes)) return;
			for(int j = 0; j < attributeRecords.length; j++)
				attributeBuffer.put(traceAttributes[j], 0, numTraces);
			attributeBuffer.clean();
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyLine2d.outputBlock.process() failed.", e, StsException.WARNING);
			return;
		}
		finally
		{
			if(attributeBuffer != null)attributeBuffer.close();
		}
	}

	public void setVerticalUnits(byte vUnits)
	{
		StsParameters.getDistanceUnitString(vUnits);
	}

	public void setZDomain(byte zDomain)
	{
		this.zDomain = StsParameters.TD_ALL_STRINGS[zDomain];
	}

    public StsObject getAssocLine() { return assocLine;}
    public void setAssocLine( StsObject assocLine) { this.assocLine = assocLine;}
}