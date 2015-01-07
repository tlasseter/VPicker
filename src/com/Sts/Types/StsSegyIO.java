package com.Sts.Types;

import com.Sts.IO.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

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

public class StsSegyIO
{
    String name;
    StsSegyData segyData;
	StsSEGYFormat segyFormat;
    int fileHeaderSize;
	int traceHeaderSize;
	int bytesPerTrace;
	int nTotalTraces;
    int bytesPerSample;
    int samplesPerTrace;
    int sampleBytesPerTrace;
    boolean isLittleEndian;
    int sampleFormat;

//    double dataAvg;
    public float dataMin = Float.MAX_VALUE;
	public float dataMax = -Float.MAX_VALUE;
    public float scale;
    public float userNull = Float.MAX_VALUE;

    private RandomAccessFile randomAccessSegyFile;

    StsSEGYFormatRec dataScaleRec;
    StsSEGYFormatRec traceTypeRec;

    public float[] histogramSamples;
    int maxNHistogramSamples;
    int histogramSampleInterval;
    int histogramIntervalCount;
    int nHistogramSamples;
    boolean scanHistogram = false;
	SampleProcess sampleProcess;

    private float floatSample;
	private double total = 0.0;
	double stepIncrement = 1.0;
	protected byte[] byteData;
    protected float[] floatData;
    protected float[] goodFloatData;
    protected byte[] traceHeader;
	//protected boolean fmtOkay = true;

    private int nTracesToScan;
	private long ioffs[];

    private int[] traceTypes = new int[100];
	private int typeCount = 0;

    protected double progressValue = 0.0;
	private boolean canceled = false;
    StsProgressPanel progressPanel;
    
    protected int analysisResult = ANALYSIS_OK;
    public static final int ANALYSIS_OK = 0;
	public static final int ANALYSIS_INVALID_FORMAT = 1;
	public static final int ANALYSIS_INVALID_ENDIAN = 2;
	public static final int ANALYSIS_IO_EXCEPTION = 3;
	public static final int ANALYSIS_UNKNOWN_EXCEPTION = 4;
    public static final int ANALYSIS_CANCELED = -1;
    private String progressDescription = "";

	static public final byte nullByte = StsParameters.nullByte;
    static public final byte nullMinByte = StsParameters.nullMinByte;
    static public final byte nullMaxByte = StsParameters.nullMaxByte;

    static public final int minNSamplesToScan = 10000;
    static public final int minNTracesToScan = 10;

    public StsSegyIO(StsSegyData segyData)
    {
        this.segyData = segyData;
        initialize();
    }

    // copy values over here for convenience  
    private void initialize()
	{
        name = segyData.dataSet.getName();
        segyFormat = segyData.segyFormat;
        fileHeaderSize = segyData.fileHeaderSize;
	    traceHeaderSize = segyData.traceHeaderSize;
	    bytesPerTrace = segyData.bytesPerTrace;
	    nTotalTraces = segyData.nTotalTraces;
        samplesPerTrace = segyData.nSamples;
        bytesPerSample = segyData.bytesPerSample;
        isLittleEndian = segyData.isLittleEndian;
        sampleBytesPerTrace = samplesPerTrace * bytesPerSample;
        sampleFormat = segyData.getSampleFormat();
        randomAccessSegyFile = segyData.randomAccessSegyFile;

        dataMin = Float.MAX_VALUE;
	    dataMax = -Float.MAX_VALUE;
//        scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
		if(sampleFormat == StsSEGYFormat.FLOAT8 || sampleFormat == StsSEGYFormat.FLOAT16)
			dataScaleRec = segyFormat.getTraceRecFromUserName("B2SCALCO");
        traceTypeRec = segyFormat.getTraceRecFromUserName("TRC_TYPE");
        initializeProcessSample();
        scanHistogram = false;
        userNull = segyData.dataSet.userNull;
    }

    public void initializeDataRange(float dataMin, float dataMax)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
    }

    public void initializeScanHistogramSamples(long nTracesToScan)
    {
        scanHistogram = true;
        maxNHistogramSamples = (int)Math.min(nTracesToScan* samplesPerTrace, minNSamplesToScan);
        histogramSampleInterval = 1;
        histogramSamples = new float[maxNHistogramSamples];
        nHistogramSamples = 0;
        histogramIntervalCount = 0;
    }

    public void initializeVolumeHistogramSamples()
    {
        scanHistogram = true;
        maxNHistogramSamples = minNSamplesToScan;
        histogramSampleInterval = Math.max(1, (int)((float)nTotalTraces*samplesPerTrace/minNSamplesToScan));
        histogramSamples = new float[maxNHistogramSamples];
        nHistogramSamples = 0;
    }

    protected void accumulateHistogramSamples(float value)
    {
        if(nHistogramSamples >= maxNHistogramSamples) return;
        if(histogramIntervalCount == histogramSampleInterval)
        {
            histogramSamples[nHistogramSamples++] = value;
            histogramIntervalCount = 1;
        }
        else
            histogramIntervalCount++;
    }

    protected void accumulateHistogramSamples(float[] values, int nValues)
    {
        int nCopy = Math.min(nValues, maxNHistogramSamples - nHistogramSamples);
        if(nCopy <= 0) return;
        System.arraycopy(values, 0, histogramSamples, nHistogramSamples, nCopy);
        nHistogramSamples += nCopy;
    }

    public boolean processTrace(byte[] segyTraceBytes, StsMappedFloatBuffer rowFloatBuffer, int nTrace)
	{
		int n = 0, pos = 0, s = 0;

		try
		{
			int offset = traceHeaderSize;
            for(s = 0; s < samplesPerTrace; s++)
			{
				sampleProcess.processSample(segyTraceBytes, offset, rowFloatBuffer);
				offset += bytesPerSample;
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyIO.processTrace() failed at  trace " + nTrace + " sample " + s + " pos " + pos, e, StsException.WARNING);
			return false;
		}
	}

	public void processSampleBytes(byte[] segyTraceBytes, int offset, StsMappedFloatBuffer rowFloatBuffer)
	{
		sampleProcess.processSample(segyTraceBytes, offset, rowFloatBuffer);
	}

	private void initializeProcessSample()
	{
        switch(sampleFormat)
		{
			case StsSEGYFormat.IBMFLT:
				sampleProcess = new sampleProcessIBMfloatSample();
				break;
			case StsSEGYFormat.IEEEFLT:
				sampleProcess = new sampleProcessIEEEfloatSample();
				break;
			case StsSEGYFormat.INT2:
				sampleProcess = new sampleProcessInt2ToBytes();
				break;
			case StsSEGYFormat.INT4:
				sampleProcess = new sampleProcessInt4ToBytes();
				break;
			case StsSEGYFormat.BYTE:
				sampleProcess = new sampleProcessByteToBytes();
				break;
			case StsSEGYFormat.FLOAT8:
				sampleProcess = new sampleProcessFloatByteToBytes();
				break;
			case StsSEGYFormat.FLOAT16:
				sampleProcess = new sampleProcessFloatInt2ToBytes();
				break;
			default:
				new StsMessage(null, StsMessage.WARNING, "Sample processing failed: unsupported sample format " + sampleFormat);
				break;
		}
	}

	public boolean processTrace(byte[] segyTraceBytes, int nSamples, int nTrace, byte[] bytes)
	{
		try
		{
			float[] floats = new float[nSamples];
			if(!processTrace(segyTraceBytes, nTrace, floats)) return false;
			for(int n = 0; n < nSamples; n++)
				bytes[n] = StsMath.unsignedIntToUnsignedByte254(Math.round((floats[n] - dataMin) * scale));
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

	public boolean processTrace(byte[] segyTraceBytes, int nTrace, float[] floats)
	{
		try
		{
			int offset = traceHeaderSize;
            for(int n = 0; n < samplesPerTrace; n++, offset += bytesPerSample)
				floats[n] = sampleProcess.processSample(segyTraceBytes, offset);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

	public boolean processTrace(byte[] segyTraceBytes, int nSamples, int nTrace, float[] floats)
	{
		try
		{
			int offset = traceHeaderSize;
            for(int n = 0; n < nSamples; n++, offset += bytesPerSample)
				floats[n] = sampleProcess.processSample(segyTraceBytes, offset);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
			return false;
		}
	}

    /** Converts trace floats to scaled trace bytes */
    public boolean processTrace(float[] floats, int nSamples, int nTrace, byte[] bytes)
    {
        try
        {
            for(int n = 0; n < nSamples; n++)
                bytes[n] = StsMath.unsignedIntToUnsignedByte254(Math.round((floats[n] - dataMin) * scale));
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
            return false;
        }
    }

    abstract class SampleProcess
    {
		abstract float processSample(byte[] bytes, int offset);

		void processSample(byte[] bytes, int offset, StsMappedFloatBuffer rowFloatBuffer)
		{
			float sampleFloat = processSample(bytes, offset);
			rowFloatBuffer.put(sampleFloat);
//			if(inlineBuffer != null)checkSetFloatToByte(inlineBuffer, sampleFloat);
			accumulateHistogramSamples(sampleFloat);
		}

		public boolean processTrace(byte[] segyTraceBytes, int offset, int nSamples, int nTrace, double[] doubles)
		{
			try
			{
                for(int n = 0; n < nSamples; n++, offset += bytesPerSample)
					doubles[n] = sampleProcess.processSample(segyTraceBytes, offset);
				return true;
			}
			catch(Exception e)
			{
				StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
				return false;
			}
		}
		public boolean processTrace(byte[] segyTraceBytes, int offset, int nSamples, int nTrace, float[] floats)
		{
			try
			{
                for(int n = 0; n < nSamples; n++, offset += bytesPerSample)
					floats[n] = sampleProcess.processSample(segyTraceBytes, offset);
				return true;
			}
			catch(Exception e)
			{
				StsException.outputException("StsSegyVolume.processTrace() failed at  trace " + nTrace, e, StsException.WARNING);
				return false;
			}
		}
	}

	class sampleProcessIBMfloatSample extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			return StsMath.convertIBMFloatBytes(bytes, offset, isLittleEndian);
		}
	}

	class sampleProcessIEEEfloatSample extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			return Float.intBitsToFloat(StsMath.convertIntBytes(bytes, offset, isLittleEndian));
		}
	}

	class sampleProcessInt2ToBytes extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			return(float)StsMath.convertBytesToShort(bytes, offset, isLittleEndian);
		}
	}

	class sampleProcessInt4ToBytes extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			System.out.println("Do we need to add support for a 4 byte integer?????? Why?");
			return StsParameters.nullValue;
		}
	}

	class sampleProcessByteToBytes extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			return(float)bytes[offset];
		}
	}

	class sampleProcessFloatInt2ToBytes extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			int scaledInt = (int)StsMath.convertBytesToShort(bytes, offset, isLittleEndian);
			int dataScale = (int)dataScaleRec.getHdrValue(bytes, isLittleEndian);
			return(float)scaledInt / (float)(Math.pow(2.0f, dataScale));
		}
	}

	class sampleProcessFloatByteToBytes extends SampleProcess
    {
		float processSample(byte[] bytes, int offset)
		{
			int dataScale = (int)dataScaleRec.getHdrValue(bytes, isLittleEndian);
			return(float)bytes[offset] / (float)(Math.pow(2.0f, dataScale));
		}
	}

	private StsAnalyzeTraces analyzeTraces = null;

	public int analyzeTraces(RandomAccessFile randomAccessSegyFile, float scanPercentage, StsProgressPanel progressPanel)
	{
		try
		{
            this.progressPanel = progressPanel;
            nTracesToScan = (int)(nTotalTraces * (scanPercentage / 100.0f));
            if(nTracesToScan < minNTracesToScan)
                nTracesToScan = Math.min(minNTracesToScan, nTotalTraces);
            
            if( nTracesToScan <= 0) return ANALYSIS_INVALID_FORMAT;

            analyzeTraces = constructTraceAnalyzer();

			analysisResult = analyzeTraces.analyze(progressPanel);
			if (analysisResult == ANALYSIS_OK)
			{
				dataMin = analyzeTraces.getDataMin();
				dataMax = analyzeTraces.getDataMax();
			}
            else
            {
                dataMin = -StsParameters.largeFloat;
                dataMax = StsParameters.largeFloat;
            }
            return analysisResult;
		}
		catch(Exception ex)
		{
			progressDescription = ex.getMessage();
			return ANALYSIS_UNKNOWN_EXCEPTION;
		}
		finally
		{
			analyzeTraces = null;
		}
	}
/*
	public int analyzeTraces(RandomAccessFile randomAccessSegyFile, int nTracesToScan, float scanPercentage)
	{
		try
		{
            this.nTracesToScan = nTracesToScan;

			if(nTracesToScan <= 0)
				nTracesToScan = Math.max(((int)((float)nTotalTraces * (scanPercentage / 100.0f))), 10); // Default to ~1% of traces

			analyzeTraces = constructTraceAnalyzer();

			analyzeTraces.analyze();
			analysisResult = analyzeTraces.getAnalysisResult();
			if (analysisResult == ANALYSIS_OK)
			{
				dataMin = analyzeTraces.getDataMin();
				dataMax = analyzeTraces.getDataMax();
			}
			return analysisResult;
		}
		catch(Exception ex)
		{
			progressDescription = ex.getMessage();
			return ANALYSIS_UNKNOWN_EXCEPTION;
		}
		finally
		{
			analyzeTraces = null;
		}
	}
*/
    static public int getBytesPerSample(int sampleFormat)
    {
        switch(sampleFormat)
        {
            case StsSEGYFormat.BYTE: // format type 8: 8 bit integer
                return 1;
            case StsSEGYFormat.INT2: // format type 3: 2 byte integer
                return 2;
            case StsSEGYFormat.INT4: // format type 2: 4 byte integer
                 // FALL_THROUGH
            case StsSEGYFormat.IBMFLT: // format type 1: 4 byte IBM floating point
                // FALL_THROUGH
            case StsSEGYFormat.IEEEFLT: // format type 5: 4 byte IEEE Float
                return 4;
            case StsSEGYFormat.FLOAT8: // format type 6: LGC - Float8 (int8 with scalar applied)
                return 1;
            case StsSEGYFormat.FLOAT16: // format type 7: LGC - Float8 (int16 with scalar applied)
                return 2;
            case StsSEGYFormat.FIXED: // format type 4: 4 byte fixed point w/gain code (obsolete)
                return 4;
            default:
 //               System.err.println("Format " + sampleFormat + " unrecognized");
                return 0;
        }
    }

	StsAnalyzeTraces constructTraceAnalyzer()
	{
		StsAnalyzeTraces analyze = null;
		switch(sampleFormat)
		{
			case StsSEGYFormat.BYTE: // 8 bit integer
				analyze = new AnalyzeByte();
				break;
			case StsSEGYFormat.INT2: // 2 byte integer
				analyze = new AnalyzeInt2();
				break;
			case StsSEGYFormat.INT4: // 4 byte integer
				analyze = new AnalyzeInt4();
				break;
			case StsSEGYFormat.IBMFLT: // 4 byte IBM floating point
				analyze = new AnalyzeIBMFloat();
				break;
			case StsSEGYFormat.IEEEFLT: // 4 byte IEEE Float
				analyze = new AnalyzeIEEEFloat();
				break;
			case StsSEGYFormat.FLOAT8: // 8 bit integer with scalar
				analyze = new AnalyzeFloat8(traceHeaderSize);
				break;
			case StsSEGYFormat.FLOAT16: // 16 bit integer with scalar
				analyze = new AnalyzeFloat16(traceHeaderSize);
				break;
		}
        return analyze;
	}

    abstract class StsAnalyzeTraces
    {
        protected StsAnalyzeTraces()
        {
        }

        private void initialise()
        {
            initializeScanHistogramSamples(nTracesToScan);

            traceTypes = null;
            traceTypes = new int[100];
            typeCount = 0;

            sampleBytesPerTrace = samplesPerTrace*bytesPerSample;
            bytesPerTrace = sampleBytesPerTrace + traceHeaderSize;
            byteData = new byte[sampleBytesPerTrace];
            floatData = new float[samplesPerTrace];
            goodFloatData = new float[samplesPerTrace];
            traceHeader = new byte[traceHeaderSize];
//            data = new byte[nSampleBytesPerTrace + traceHeaderSize];
            ioffs = new long[nTracesToScan];
            long ioff = fileHeaderSize;
            // Build the offset array to scan
            if((nTotalTraces > 0 && nTotalTraces <= nTracesToScan))
            {
                nTracesToScan = nTotalTraces;
                for(int i = 0; i < nTracesToScan; i++)
                {
                    ioffs[i] = ioff + (long)i * bytesPerTrace;
                }
            }
            else
            {
                Random random = new Random();
                for(int i = 0; i < nTracesToScan; i++)
                {
                    ioffs[i] = ioff + (long)(random.nextFloat() * nTotalTraces) * bytesPerTrace;
                }
                Arrays.sort(ioffs); //Probably faster if sorted
            }
        }

        public int analyze(StsProgressPanel progressPanel)
        {
            try
            {
                initialise();

                StsIntervalProgressBar progressBar = progressPanel.progressBar;
                progressBar.setIntervalCount(nTracesToScan);
                progressBar.setDescriptionAndLevel("Analyzing traces for volume " + name, StsProgressBar.INFO);
//                double traceAverageSum = 0.0;
                int nNonNullSamples = 0;
                int nTotalSamples = 0;
                for (int t = 0; t < nTracesToScan; t++)
                {
                    randomAccessSegyFile.seek(ioffs[t]);
                    randomAccessSegyFile.read(traceHeader);
                    randomAccessSegyFile.read(byteData);

                    accumulateTypes();

                    if (progressBar.isCanceled())
                    {
                        dataMin = Float.MAX_VALUE;
	                    dataMax = -Float.MAX_VALUE;
                        return ANALYSIS_CANCELED;
                    }
                    getSamples();
                    int nGoodSamples = 0;
                    for(int s = 0; s < samplesPerTrace; s++)
                    {
                        float floatSample = floatData[s];
                        if (floatSample != -Float.MAX_VALUE && floatSample != Float.MAX_VALUE && floatSample != userNull)
                        {
                            if (floatSample < dataMin)
                                dataMin = floatSample;
                            else if (floatSample > dataMax)
                                dataMax = floatSample;
                            goodFloatData[nGoodSamples++] = floatSample;
                            nNonNullSamples++;
                        }

                        if (canceled || analysisResult != ANALYSIS_OK)//!fmtOkay)
                        {
                            //fmtOkay = false;
                            return analysisResult;
                        }
                        nTotalSamples++;
                    }
                    accumulateHistogramSamples(goodFloatData, nGoodSamples);
                    if(analysisResult != ANALYSIS_OK) return analysisResult;
                    progressBar.incrementCount();
                }

                traceTypes = (int[])StsMath.trimArray(traceTypes, typeCount);
                float percentNull = 100*((float)nTotalSamples - nNonNullSamples)/nTotalSamples;
                // System.out.println("nTotalSamples: " + nTotalSamples + " percent null: " + percentNull);
                return analysisResult;
            }
            catch(IOException ioe)
            {
                analysisResult = ANALYSIS_IO_EXCEPTION;
                //fmtOkay = false;
                StsException.systemError("StsAnalyzeTraces.analyze() failed to read randomAccessFile");
                return analysisResult;
            }
            catch(Exception e)
            {
                analysisResult = ANALYSIS_UNKNOWN_EXCEPTION;
                //fmtOkay = false;
                StsException.outputException("StsAnalyzeTraces.analyze() failed.", e, StsException.WARNING);
                return analysisResult;
            }
        }

        private void getCheckFloatSample(float floatSample)
        {
            if (floatSample != -Float.MAX_VALUE && floatSample != Float.MAX_VALUE && floatSample != userNull)
            {
                if (floatSample < dataMin)
                    dataMin = floatSample;
                else if (floatSample > dataMax)
                    dataMax = floatSample;
            }
        }
/*
        public void analyze()
        {
            try
            {
                initialise();

                stepIncrement = 100.0 / (nTracesToScan);
                progressValue = 0.0;
                progressDescription = "";
                for (int t = 0; t < nTracesToScan; t++)
                {
                    randomAccessSegyFile.seek(ioffs[t]);
                    randomAccessSegyFile.read(traceHeader);
                    randomAccessSegyFile.read(byteData);

                    accumulateTypes();

                    if (canceled || analysisResult != ANALYSIS_OK)//!fmtOkay)
                    {
                        //fmtOkay = false;
                        return;
                    }
                    getSamples();
                    int s = 0;
                    if(t == 0)
                    {
                        floatSample = floatData[0];
                        if ((floatSample == -Float.MAX_VALUE) || (floatSample == Float.MAX_VALUE))
                        {
                            floatSample = 0.0f;
                            floatData[0] = 0.0f;
                        }
                        dataMin = floatSample;
                        dataMax = floatSample;
                        s = 1;
                    }
                    for(; s < samplesPerTrace; s++)
                    {
                        floatSample = floatData[s];
                        if ((floatSample == -Float.MAX_VALUE) || (floatSample == Float.MAX_VALUE))
                        {
                            floatSample = 0.0f;
                            floatData[s] = 0.0f;
                        }
                        if (floatSample < dataMin)
                            dataMin = floatSample;
                        else if (floatSample > dataMax)
                            dataMax = floatSample;
                    }
                    accumulateHistogramSamples(floatData);
                    progressValue += stepIncrement;
                    progressDescription = " Scanned Trace " + t + " of " + nTracesToScan + " to be scanned.";

                    // Truncate the type area
                }
                traceTypes = (int[])StsMath.trimArray(traceTypes, typeCount);

                progressDescription = "";
                progressValue = 100.0;
            }
            catch(IOException ioe)
            {
                analysisResult = ANALYSIS_IO_EXCEPTION;
                //fmtOkay = false;
                StsException.systemError("StsAnalyzeTraces.analyze() failed to read randomAccessFile");
            }
            catch(Exception e)
            {
                analysisResult = ANALYSIS_UNKNOWN_EXCEPTION;
                //fmtOkay = false;
                StsException.outputException("StsAnalyzeTraces.analyze() failed.", e, StsException.WARNING);
            }
        }

        public void resetProgressBar()
        {
            progressDescription = "";
            progressValue = 0.0;
        }
*/
        public float getDataMin()
        {
            return dataMin;
        }

        public float getDataMax()
        {
            return dataMax;
        }

        public void setDataMin(float dataMin_)
        {
            dataMin = dataMin_;
        }

        public void setDataMax(float dataMax_)
        {
            dataMax = dataMax_;
        }

        public int[] getTraceTypes()
        {
            return traceTypes;
        }

        protected abstract float[] getSamples();
        protected abstract float getSample(int index);

        public double getTotal()
        {
            return total;
        }

        /*public boolean getFmtOkay()
        {
            return fmtOkay;
        }*/

        public void accumulateTypes()
        {
            int tid = -1;
            boolean found = false;

            // Way too many types to worry about so skip it.
            if(typeCount == traceTypes.length)
                return;

            // Collect all the traceTypes
            tid = (int)traceTypeRec.getHdrValue(traceHeader, isLittleEndian);
            // Could just be random trash
            if((tid > 20) || (tid < 1))
                return;

            for(int n = 0; n < typeCount; n++)
            {
                if(traceTypes[n] == tid)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
                traceTypes[typeCount++] = tid;
        }

        public void run()
        {
        }

        public double getProgressValue()
        {
            return progressValue;
        }

        public String getProgressDescription()
        {
            return progressDescription;
        }

        public int getAnalysisResult()
        {
            return analysisResult;
        }

        public boolean getFmtOkay()
        {
            return analysisResult == ANALYSIS_OK;
        }
    }

    class AnalyzeByte extends StsAnalyzeTraces
    {
        AnalyzeByte()
        {
            super();
            bytesPerSample = 1;
        }

        protected float getSample(int index)
        {
            return(float)byteData[index];
        }

        protected float[] getSamples()
        {
            for(int n = 0; n < samplesPerTrace; n++)
                floatData[n] = (float)byteData[n];
            return floatData;
        }
    }

    class AnalyzeInt2 extends StsAnalyzeTraces
    {
        AnalyzeInt2()
        {
            super();
            bytesPerSample = 2;
        }

        protected float getSample(int index)
        {
            return(float)StsMath.shortValue(index, byteData, isLittleEndian);
        }

        protected float[] getSamples()
        {
            if(!StsMath.convertShortBytesToFloats(byteData, floatData,  0, samplesPerTrace, isLittleEndian)) return null;
            return floatData;
        }
    }

    class AnalyzeInt4 extends StsAnalyzeTraces
    {
        AnalyzeInt4()
        {
            super();
            bytesPerSample = 4;
        }

        protected float getSample(int index)
        {
            return(float)StsMath.intValue(index, byteData, isLittleEndian);
        }

        protected float[] getSamples()
        {
            if(!StsMath.convertIntBytes(byteData, floatData, samplesPerTrace, isLittleEndian))
                return null;
            return floatData;
        }
    }

    class AnalyzeIBMFloat extends StsAnalyzeTraces
    {
        float[] fdata;

        AnalyzeIBMFloat()
        {
            super();
            bytesPerSample = 4;
        }

        protected float getSample(int index)
        {
            int offset = bytesPerSample*index;
            return StsMath.convertIBMFloatBytes(byteData, offset);
        }

        protected float[] getSamples()
        {
           if(!StsMath.verifyIBMFloatType(byteData, isLittleEndian))
            {
                progressDescription = " Error with IBM Float format. Try IEEE Float.";
                analysisResult = ANALYSIS_INVALID_FORMAT;
                //fmtOkay = false;
                return null;
            }
            StsMath.convertIBMFloatBytes(byteData, floatData, samplesPerTrace, isLittleEndian);
            return floatData;
        }
    }

    class AnalyzeIEEEFloat extends StsAnalyzeTraces
    {
        float[] fdata;

        AnalyzeIEEEFloat()
        {
            super();
            bytesPerSample = 4;
        }

        protected float getSample(int index)
        {
            return Float.intBitsToFloat(StsMath.intValue(index, byteData, isLittleEndian));
        }

        protected float[] getSamples()
        {
            StsMath.convertIEEEBytesToFloats(byteData, floatData, samplesPerTrace, isLittleEndian);
            return floatData;
        }
    }

    class AnalyzeFloat8 extends StsAnalyzeTraces
    {
        AnalyzeFloat8(int traceHeaderSize)
        {
            super();
            bytesPerSample = 1;
        }

        private float getDataScale()
        {
            int intScale = (int)dataScaleRec.getHdrValue(traceHeader, isLittleEndian);
            return 1/(float)(Math.pow(2.0f, (float)intScale));
        }

        protected float getSample(int index)
        {
            return(float)byteData[index]*getDataScale();
        }

        protected float[] getSamples()
        {
            float dataScale = getDataScale();
            for(int n = 0; n < samplesPerTrace; n++)
                floatData[n] = (float)byteData[n]*dataScale;
            return floatData;
        }
    }

    class AnalyzeFloat16 extends StsAnalyzeTraces
    {
        AnalyzeFloat16(int traceHeaderSize)
        {
            super();
            bytesPerSample = 2;
        }

        private float getDataScale()
        {
            int intScale = (int)dataScaleRec.getHdrValue(traceHeader, isLittleEndian);
            return 1/(float)(Math.pow(2.0f, (float)intScale));
        }

        protected float getSample(int index)
        {
            return(float)StsMath.shortValue(index, byteData, isLittleEndian) * getDataScale();
        }

        protected float[] getSamples()
        {
           float dataScale = getDataScale();
           StsMath.convertShortBytesToFloats(byteData, floatData, 0, samplesPerTrace, isLittleEndian);
           for(int n = 0; n < samplesPerTrace; n++)
                floatData[n] *= dataScale;
           return floatData;
        }
    }
}