package com.Sts.Types;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.IO.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import sun.io.*;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 28, 2007
 * Time: 7:52:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSegyData
{
    public StsSeismicWizard seismicWizard;
    public StsSeismicBoundingBox dataSet;
    public StsSEGYFormat segyFormat;

    private byte[] textHeader;
	public byte[] binaryHeader;
	public byte[] traceHeader;
	public int textHeaderSize = StsSEGYFormat.defaultTextHeaderSize;
	public int binaryHeaderSize = StsSEGYFormat.defaultBinaryHeaderSize;
	public int traceHeaderSize = StsSEGYFormat.defaultTraceHeaderSize;
	public int bytesPerSample = 4;
    public int fileHeaderSize = 0;
	public int bytesPerTrace;
    public int nSamples;
	public int headerSampleFormat = StsSEGYFormat.defaultSampleFormat;
    public int overrideSampleFormat = StsSEGYFormat.NONE;
 //   public boolean isSampleFormatOverride = false;
	public int nTotalTraces = 0;
    private float sampleSpacing = 0.0f;
//    public float dataMin;
//    public float dataMax;
    public boolean override = false;
    public int overrideNSamples = 0;
    public float overrideSampleSpacing = 0.0f;
    public float startZ = 0.0f;
    public boolean isLittleEndian = false;
    public String textHeaderFormatString = StsSEGYFormat.EBCDIC;
    public long segyFileSize;

    public float verticalScalar = 1.0f;
	public float horizontalScalar = 1.0f;
    public int numVolumes = 1;
    /** Indicates that segy volume contains multiple volumes (multi-component data for example) */
    public boolean isMultiVolumeFile = false;

    public RandomAccessFile randomAccessSegyFile;

    public final String[] sampleFormatStrings = StsSEGYFormat.sampleFormatStrings;

    public StsSegyData()
    {
    }

    static public StsSegyData constructor(StsSeismicWizard seismicWizard, StsFile file, StsSeismicBoundingBox data, StsSEGYFormat segyFormat)
    {
        StsSegyData segyData = new StsSegyData();
        if(!segyData.initialize(seismicWizard, file, data, segyFormat))
            return null;
        else
            return segyData;
    }

    public boolean initialize(StsSeismicWizard seismicWizard, StsFile file, StsSeismicBoundingBox data, StsSEGYFormat segyFormat)
    {
		try
		{
            this.seismicWizard = seismicWizard;
            this.dataSet = data;
            randomAccessSegyFile = new RandomAccessFile(file.getPathname(), "r");
            segyFileSize = randomAccessSegyFile.length();
            setSegyFormat(segyFormat);
            return true;
        }
		catch(Exception e)
		{
			StsException.outputWarningException(this, "initialize", e);
            return true;
        }
    }

    public void initializeFileHeader()
    {
		try
		{
            readFileHeader();
            initializeBinaryHdr();
        }
		catch(Exception e)
		{
			StsException.outputWarningException(this, "initializeFileHeader", e);
        }
    }

    public void readFileHeader() throws IOException
	{
        int textHeaderSize = segyFormat.textHeaderSize;
		textHeader = new byte[textHeaderSize];
		randomAccessSegyFile.seek(0);
		randomAccessSegyFile.read(textHeader);

		binaryHeaderSize = segyFormat.binaryHeaderSize;
		binaryHeader = new byte[binaryHeaderSize];
		randomAccessSegyFile.seek(textHeaderSize);
		randomAccessSegyFile.read(binaryHeader);
		fileHeaderSize = (textHeaderSize + binaryHeaderSize); // + # extended file text header sizes
    }

    public boolean initializeBinaryHdr() //throws IOException
	{

        return analyzeBinaryHdr(null);
    }

    //TODO add progressPanel as argument and output lines to it with more info on failures (return false)
    public boolean analyzeBinaryHdr(StsProgressPanel progressPanel) //throws IOException
	{
//        isLittleEndian = segyFormat.getIsLittleEndian();
//****wrw this overrides menu selection and may be bogus
//        sampleFormat = segyFormat.getSampleFormat(isLittleEndian, binaryHeader);
//        setBytesPerSampleFromFormat(sampleFormat);
        //For NMO, need sample interval in sec, so convert from ms.
		if(!override)
		{
            sampleSpacing = segyFormat.getSampleSpacing(isLittleEndian, binaryHeader);
            /*  This causes problems as switching endianness makes it positive in some cases, but it's still garbarge.
             *  So the user overrides the garbage, but doesn't know that he needs to reset the Endianness.  If he doesn't reset,
             *  everything else is garbage.  So I'm commenting it out.  TJL 4/8/08
             */
            /*
            if(sampleSpacing <= 0.0)
            {
            	// Determine if the problem is endianness.
            	isLittleEndian = !isLittleEndian;
            	sampleSpacing = segyFormat.getSampleSpacing(isLittleEndian, binaryHeader);
            	if(sampleSpacing <= 0.0) 
            	{
                    if(progressPanel != null) progressPanel.appendLine("    Volume " + dataSet.getName() + " sampleSpacing in error: " + sampleSpacing);
                    isLittleEndian = !isLittleEndian;  // Not an endianess problem so flip back to user value.
            		return false;      
            	}
            	// Set new endianness in segyformat since it is used on trace definition page.
                // Nope.  Each file could have a different littleEndianness.  segyFormat is the initial and saved value. TJL 3/30/08
                // Each segyData carries its own copy of isLittleEndian which is used in all subsequent processing of this file.
                // segyFormat.setIsLittleEndian(isLittleEndian);
            }
            */
			nSamples = segyFormat.getNSamples(isLittleEndian, binaryHeader);
            overrideNSamples = nSamples;
            overrideSampleSpacing = sampleSpacing;
        }
        else
        {
            nSamples = overrideNSamples;
            sampleSpacing = overrideSampleSpacing;
           	if(sampleSpacing <= 0.0)
            {
                if(progressPanel != null) progressPanel.appendLine("    Volume " + dataSet.getName() + " sampleSpacing in error: " + sampleSpacing);
                return false;
            }
        }

        dataSet.zInc = verticalScalar * sampleSpacing;
        if(nSamples <= 0) return false;
        dataSet.nSlices = nSamples;
        // Not sure why this is being set to 0...messing with pre-stack wizard. JKF
		//zMin = 0.0f;
		//Only setting to zero if it has never been set
		if (dataSet.zMin >= StsParameters.largeFloat)
			dataSet.zMin = 0.0f;
		dataSet.zMax = dataSet.zMin + (dataSet.nSlices - 1) * dataSet.zInc;

 //       isLittleEndian = segyFormat.getIsLittleEndian();
        headerSampleFormat = (int)segyFormat.getSampleFormat(isLittleEndian, binaryHeader);
        setHeaderSampleFormat(headerSampleFormat);
        bytesPerTrace = nSamples * bytesPerSample + traceHeaderSize;
        nTotalTraces = (int)((segyFileSize - fileHeaderSize) / bytesPerTrace);
        isMultiVolumeFile = segyFormat.isMultiVolume(isLittleEndian, binaryHeader);
        if(isMultiVolumeFile) numVolumes = (int)segyFormat.getBinaryHdrValue("NUMVOLS", binaryHeader, isLittleEndian);
		// When the user enters the file format screen the sample format is set to the value in the header. If they manually set
        // it to something else, the userOverrideFormat is set to true and the header value is ignored.

//        overrideSampleFormat = segyFormat.overrideSampleFormat;
        // sample format assumed to be a 4-byte float
//		dataMin = Float.MAX_VALUE;
//		dataMax = -dataMin;
//		bytesPerSample = 4;
        dataSet.zDomain = segyFormat.getZDomainString();
        return true;
	}
/*
    public void setIsSampleFormatOverride(boolean b)
    {
        if(isSampleFormatOverride == b) return;
        isSampleFormatOverride = b;
        if(!isSampleFormatOverride) headerSampleFormat = (int)segyFormat.getSampleFormat(isLittleEndian, binaryHeader);
    }
*/
    public StsSEGYFormat getSegyFormat()
    {
        return segyFormat;
    }

    public void setSegyFormat(StsSEGYFormat segyFormat)
    {
        if(!copySegyFormat(segyFormat)) return;
        initializeFromSegyFormat();
    }

    private boolean copySegyFormat(StsSEGYFormat segyFormat)
    {

        this.segyFormat = (StsSEGYFormat)StsToolkit.copy(segyFormat);
        return this.segyFormat != null;
    }

    private void initializeFromSegyFormat()
    {
        isLittleEndian = segyFormat.getIsLittleEndian();
        textHeaderFormatString = segyFormat.getTextHeaderFormatString();
        traceHeaderSize = segyFormat.getTraceHeaderSize();
        traceHeader = new byte[traceHeaderSize];
        dataSet.resetStatus();
        initializeFileHeader();
    }

    public byte[] getTextHeader()
    {
        return textHeader;
    }

    public boolean getIsLittleEndian()
    {
        return isLittleEndian;
    }

    public void setIsLittleEndian(boolean isLittleEndian)
    {
        if(this.isLittleEndian == isLittleEndian) return;
        this.isLittleEndian = isLittleEndian;
        dataSet.resetStatus();
    }

    public void setOverrideHeader(boolean b)
    {
        override = b;
        updateTraceGeometry();
        dataSet.resetStatus();
    }

    public boolean getOverrideHeader() { return override; }

    public void setOverrideNSamples(int nSamples)
    {
        overrideNSamples = nSamples;
//        dataSet.setNSlices(nSamples);
        updateTraceGeometry();
        dataSet.resetStatus();
    }

    public int getOverrideNSamples() { return overrideNSamples; }

    public void setOverrideSampleSpacing(float sampleSpacing)
    {
        overrideSampleSpacing = sampleSpacing;
        updateTraceGeometry();
        dataSet.resetStatus();
    }

    public float getOverrideSampleSpacing() { return overrideSampleSpacing; }

    public void setTextHeaderFormatString(String string)
    {
        textHeaderFormatString = string;
        segyFormat.setTextHeaderFormatString(string);
    }
    public String getTextHeaderFormatString() { return textHeaderFormatString; }

    public void setStartZ(float z)
    {
        startZ = z;
        updateTraceGeometry();
    }

    private void updateTraceGeometry()
    {
        int nSlices = getNSamples();
        dataSet.setNSlices(nSlices);
        float zInc = getSampleSpacing();
        dataSet.setZInc(zInc);
        float zMin = getStartZ();
        dataSet.setZMin(zMin);
        float zMax = zMin + (nSlices-1)*zInc;
        dataSet.setZMax(zMax);
    }

    public float getStartZ() { return startZ; }

    public int getNSamples()
    {
        if(override)
            return overrideNSamples;
        else
            return nSamples;
    }

    public void setNSamples(int nSamples)
    {
        if(override)
            overrideNSamples = nSamples;
        else
            this.nSamples = nSamples;
    }

    public float getSampleSpacing()
    {
        if(override)
            return overrideSampleSpacing;
        else
            return sampleSpacing;
    }

    public void setSampleSpacing(float spacing)
    {
        if(override)
            overrideSampleSpacing = spacing;
        else
            sampleSpacing = spacing;
    }

    public String readTextHdr(String encoding)
	{
		try
		{
			if(encoding != null && encoding.equals("Cp500"))
			{
				ByteToCharCp500 converter = new ByteToCharCp500();
				return new String(converter.convertAll(textHeader));
			}
			else
				return new String(getTextHeader(), 0, textHeaderSize, encoding);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyData.readTextHdr() failed.", e, StsException.WARNING);
			return null;
		}
	}

    public void setTextHeader(byte[] textHeader)
    {
        this.textHeader = textHeader;
    }

    public byte[] getBinaryHeader()
    {
        return binaryHeader;
    }

    public void setBinaryHeader(byte[] binaryHeader)
    {
        this.binaryHeader = binaryHeader;
    }

    public byte[] getTraceHeaderBinary(int nTrace)
    {
        long offset = fileHeaderSize + (long)bytesPerTrace * (long)nTrace;
        byte[] hdrBuf = new byte[getTraceHeaderSize()];

        try
        {
            randomAccessSegyFile.seek(0);
            randomAccessSegyFile.seek(offset);
            randomAccessSegyFile.read(hdrBuf);
            return hdrBuf;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSegyData.getTraceHeaderBinary() failed.", e, StsException.WARNING);
            return null;
        }
    }

	public StsSEGYFormat.TraceHeader getTrace(int n)
	{
		int nSampleBytesPerTrace = nSamples * bytesPerSample;
		byte[] data = new byte[nSampleBytesPerTrace];

		try
		{
			if(n < 0)
				n = 0;
			else if(n >= nTotalTraces)
				n = nTotalTraces - 1;

			long ioff = (long)fileHeaderSize + (long)n * bytesPerTrace;
			randomAccessSegyFile.seek(ioff);
			randomAccessSegyFile.read(traceHeader);
//            randomAccessSegyFile.read(data);
			return segyFormat.constructTraceHeader(traceHeader, n, isLittleEndian);
		}
		catch(Exception e)
		{
			StsException.systemError("StsSegyVolume.getTrace() failed. Couldn't find trace " +
									 n + " of " + nTotalTraces + " total traces.");
			return null;
		}
	}

    public void setBytesPerSampleFromFormat()
    {
        int sampleFormat = getSampleFormat();

        switch(sampleFormat)
        {
            case StsSEGYFormat.BYTE: // format type 8: 8 bit integer
            case StsSEGYFormat.FLOAT8: // format type 6: LGC - Float8 (int8 with scalar applied)
                bytesPerSample = 1;
                return;
            case StsSEGYFormat.INT2: // format type 3: 2 byte integer
            case StsSEGYFormat.FLOAT16: // format type 7: LGC - Float8 (int16 with scalar applied)
                bytesPerSample = 2;
                return;
            case StsSEGYFormat.INT4: // format type 2: 4 byte integer
            case StsSEGYFormat.IBMFLT: // format type 1: 4 byte IBM floating point
            case StsSEGYFormat.IEEEFLT: // format type 5: 4 byte IEEE Float
            case StsSEGYFormat.FIXED: // format type 4: 4 byte fixed point w/gain code (obsolete)
                bytesPerSample = 4;
                return;
            default:
                bytesPerSample = 4;
//                StsException.systemError(this, "setBytesPerSampleFromFormat", "Illegal sample format: " + sampleFormat);
        }
    }

    public void resetDataRange()
    {

        switch(getSampleFormat())
        {
            case StsSEGYFormat.BYTE: // format type 8: 8 bit integer
                bytesPerSample = 1;
                dataSet.dataMin = Byte.MAX_VALUE;
                dataSet.dataMax = Byte.MIN_VALUE;
                break;
            case StsSEGYFormat.INT2: // format type 3: 2 byte integer
                bytesPerSample = 2;
                dataSet.dataMin = Short.MAX_VALUE;
                dataSet.dataMax = Short.MIN_VALUE;
                break;
            case StsSEGYFormat.INT4: // format type 2: 4 byte integer
                bytesPerSample = 4;
                dataSet.dataMin = Integer.MAX_VALUE;
                dataSet.dataMax = Integer.MIN_VALUE;
                break;
            case StsSEGYFormat.IBMFLT: // format type 1: 4 byte IBM floating point
                // FALL_THROUGH
            case StsSEGYFormat.IEEEFLT: // format type 5: 4 byte IEEE Float
                dataSet.dataMin = Float.MAX_VALUE;
                dataSet.dataMax = -dataSet.dataMin;
                bytesPerSample = 4;
                break;
            case StsSEGYFormat.FLOAT8: // format type 6: LGC - Float8 (int8 with scalar applied)
                dataSet.dataMin = Float.MAX_VALUE;
                dataSet.dataMax = Float.MIN_VALUE;
                bytesPerSample = 1;
                break;
            case StsSEGYFormat.FLOAT16: // format type 7: LGC - Float8 (int16 with scalar applied)
                dataSet.dataMin = Float.MAX_VALUE;
                dataSet.dataMax = Float.MIN_VALUE;
                bytesPerSample = 2;
                break;
            case StsSEGYFormat.FIXED: // format type 4: 4 byte fixed point w/gain code (obsolete)
                // FALL_THROUGH
        }
    }

	public double getBinaryHeaderValue(StsSEGYFormatRec rec)
	{
        return segyFormat.getBinaryHdrValue(rec.getUserName(), binaryHeader, isLittleEndian);
	}

	public boolean extractHeaderAttributes(byte[] segyTraceBytes, int traceNum, StsSEGYFormatRec[] attributeRecords, double[][] traceAttributes)
	{
        // Cleaned up.  scale adjustments (set to 1.0 if equal to 0.0 for example is handled in getHdrValue
//		StsSEGYFormat.TraceHeader header = null;
//		double val = 0.0f;
//		int xyScale = 1;
//        int edScale = 1;
		if(traceNum < 0) return false;
//		StsSEGYFormat.TraceHeader header = segyFormat.constructTraceHeader(traceHeader, 0, isLittleEndian);
		int xyScale = segyFormat.getCoordinateScale(segyTraceBytes, isLittleEndian);
//		if(xyScale == 0) xyScale = 1;
        int edScale = segyFormat.getElevationScale(segyTraceBytes, isLittleEndian);
//        if(edScale == 0) edScale = 1;
//		int tid = segyFormat.getTraceId(segyTraceBytes, getIsLittleEndian);
//		if((tid != traceId) && (traceId != -1)) return false;
		for(int i = 0; i < attributeRecords.length; i++)
		{
			traceAttributes[i][traceNum] = attributeRecords[i].getHdrValue(segyTraceBytes, 0, isLittleEndian, xyScale, edScale);
//            System.out.println("Attribute=" + attributeRecords[i].getName() + " Value=" + val);
//			traceAttributes[i][traceNum] = val;
		}
        return true;
    }

    public boolean extractHeaderAttributes(byte[] segyTraceBytes, StsSEGYFormatRec[] attributeRecords, double[] traceAttributes)
	{
		int xyScale = segyFormat.getCoordinateScale(segyTraceBytes, isLittleEndian);
        int edScale = segyFormat.getElevationScale(segyTraceBytes, isLittleEndian);
		for(int i = 0; i < attributeRecords.length; i++)
			traceAttributes[i] = attributeRecords[i].getHdrValue(segyTraceBytes, 0, isLittleEndian, xyScale, edScale);
        return true;
    }

    public boolean isAValidFileSize()
    {
        try
        {

            bytesPerSample = StsSegyIO.getBytesPerSample(getSampleFormat());
            long fileSize = randomAccessSegyFile.length();
            long allTracesSize = fileSize - (segyFormat.getBinaryHeaderSize() + segyFormat.getTextHeaderSize());
            long traceSize = 0;
            int traceHeaderSize = segyFormat.getTraceHeaderSize();
            int samples = 0;
            if (override)
            {
                samples = overrideNSamples;
            }
            else
            {
                samples = nSamples;
            }
            traceSize = (long)((samples * bytesPerSample) + traceHeaderSize);
            if (allTracesSize % traceSize != 0)
                return false;
        }
        catch (IOException ex)
        {
            return false;
        }
        return true;
    }

    public float getXLine(byte[] segyTraceBytes)
    {
        return segyFormat.getXLine(segyTraceBytes, isLittleEndian);
    }

    public String getName()
    {
        return dataSet.getName();
    }

    public byte[] getTraceHeader()
    {
        return traceHeader;
    }

    public void setTraceHeader(byte[] traceHeader)
    {
        this.traceHeader = traceHeader;
    }

    public int getTextHeaderSize()
    {
        return textHeaderSize;
    }

    public void setTextHeaderSize(int textHeaderSize)
    {
        if(this.textHeaderSize == textHeaderSize) return;
        this.textHeaderSize = textHeaderSize;
        segyFormat.setTextHeaderSize(textHeaderSize);
        dataSet.resetStatus();
        seismicWizard.analyzeHeaders();
    }

    public int getBinaryHeaderSize()
    {
        return binaryHeaderSize;
    }

    public void setBinaryHeaderSize(int binaryHeaderSize)
    {
        if(this.binaryHeaderSize == binaryHeaderSize) return;
        this.binaryHeaderSize = binaryHeaderSize;
        segyFormat.setBinaryHeaderSize(binaryHeaderSize);
        dataSet.resetStatus();
    }

    public int getTraceHeaderSize()
    {
        return traceHeaderSize;
    }

    public void setTraceHeaderSize(int traceHeaderSize)
    {
        if(this.traceHeaderSize == traceHeaderSize) return;
        this.traceHeaderSize = traceHeaderSize;
        traceHeader = new byte[traceHeaderSize];
        segyFormat.setTraceHeaderSize(traceHeaderSize);
        dataSet.resetStatus();
    }

    public int getBytesPerSample()
    {
        return bytesPerSample;
    }

    public void setBytesPerSample(int bytesPerSample)
    {
        this.bytesPerSample = bytesPerSample;
    }

    public int getFileHeaderSize()
    {
        return fileHeaderSize;
    }

    public void setFileHeaderSize(int fileHeaderSize)
    {
        this.fileHeaderSize = fileHeaderSize;
    }

    public int getBytesPerTrace()
    {
        return bytesPerTrace;
    }

    public void setBytesPerTrace(int bytesPerTrace)
    {
        this.bytesPerTrace = bytesPerTrace;
    }

    public int getNSamplesPerTrace()
    {
        if(override)
            return overrideNSamples;
        else
            return nSamples;
    }

    public void setNSamplesPerTrace(int nSamplesPerTrace)
    {
        this.nSamples = nSamplesPerTrace;
    }

    public int getSampleFormat()
    {
        if(overrideSampleFormat != StsSEGYFormat.NONE)
            return overrideSampleFormat;
        else if( headerSampleFormat != StsSEGYFormat.NONE)
            return headerSampleFormat;
        else
            return StsSEGYFormat.defaultSampleFormat;
    }

    public int getHeaderSampleFormat()
    {
        return headerSampleFormat;
    }

    public int getOverrideSampleFormat()
    {
        return overrideSampleFormat;
    }

    public void setHeaderSampleFormat(int headerSampleFormat)
    {
        if(headerSampleFormat == StsSEGYFormat.NONE)
            this.headerSampleFormat = StsSEGYFormat.defaultSampleFormat;
        else
            this.headerSampleFormat = headerSampleFormat;
        setBytesPerSampleFromFormat();
    }

    public void setOverrideSampleFormat(int sampleFormat)
    {
        this.overrideSampleFormat = sampleFormat;
        setBytesPerSampleFromFormat();
        dataSet.resetStatus();
//        isSampleFormatOverride = true;
    }

    private boolean isSampleFormatBad(int sampleFormat)
    {
        return sampleFormat < 1 || sampleFormat > 8;
    }

    public String getSampleFormatString()
    {
        return sampleFormatStrings[getSampleFormat()];
    }
/*
    public void setSampleFormatString(String sampleFormatString)
    {
//        sampleFormat = StsSEGYFormat.getSampleFormatFromString(sampleFormatString);
        setHeaderSampleFormat(StsSEGYFormat.getSampleFormatFromString(sampleFormatString));
    }
*/
    public long getSegyFileSize()
    {
        return segyFileSize;
    }

    public void setSegyFileSize(long segyFileSize)
    {
        this.segyFileSize = segyFileSize;
    }

    public int getNTotalTraces()
    {
        return nTotalTraces;
    }

    public long getNTotalSamples()
    {
        return ((long)nTotalTraces)* nSamples;
    }

    public void setNTotalTraces(int nTotalTraces)
    {
        this.nTotalTraces = nTotalTraces;
    }

    public RandomAccessFile getRandomAccessSegyFile()
    {
        return randomAccessSegyFile;
    }

    public void setRandomAccessSegyFile(RandomAccessFile randomAccessSegyFile)
    {
        this.randomAccessSegyFile = randomAccessSegyFile;
    }

    public float getVerticalScalar()
    {
        return verticalScalar;
    }

    public void setVerticalScalar(float verticalScalar)
    {
        this.verticalScalar = verticalScalar;
    }

    public float getHorizontalScalar()
    {
        return horizontalScalar;
    }

    public void setHorizontalScalar(float horizontalScalar)
    {
        this.horizontalScalar = horizontalScalar;
    }

    public int getNumVolumes()
    {
        return numVolumes;
    }

    public void setNumVolumes(int numVolumes)
    {
        this.numVolumes = numVolumes;
    }

    public boolean isMultiVolumeFile()
    {
        return isMultiVolumeFile;
    }

    public void setMultiVolumeFile(boolean multiVolumeFile)
    {
        this.isMultiVolumeFile = multiVolumeFile;
    }

    static public StsSegyData[] getSegyDatasets(StsSeismicBoundingBox[] volumes)
    {
        if(volumes == null) return new StsSegyData[0];
        int nVolumes = volumes.length;
        StsSegyData[] segyDataSets = new StsSegyData[nVolumes];
        for(int n = 0; n < nVolumes; n++)
            segyDataSets[n] = volumes[n].segyData;
        return segyDataSets;
    }
}
