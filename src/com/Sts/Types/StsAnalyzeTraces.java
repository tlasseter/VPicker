package com.Sts.Types;

import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
public abstract class StsAnalyzeTraces
{
	public static final int ANALYSIS_OK = 0;
	public static final int ANALYSIS_INVALID_FORMAT = 1;
	public static final int ANALYSIS_INVALID_ENDIAN = 2;
	public static final int ANALYSIS_IO_EXCEPTION = 3;
	public static final int ANALYSIS_UNKNOWN_EXCEPTION = 4;

    private String name;
    private float sample;
	private double total = 0.0;
	double stepIncrement = 1.0;
	protected byte[] data;
	protected byte[] traceHeader;
	//protected boolean fmtOkay = true;
	protected boolean isLittleEndian;
	private RandomAccessFile randomAccessSegyFile;
	private int nTracesToScan;
	protected int indexStart;
	protected int indexStep;
	private long ioffs[];
	//private boolean scanHistogram;
	protected int traceHeaderSize;
	private int fileHeaderSize;
	private int nTotalTraces;
	private int nSampleBytesPerTrace;
	protected int bytesPerSample;
	protected StsSEGYFormatRec dataScaleRec;
	protected StsSEGYFormatRec traceTypeRec;
	protected int nSamplesPerTrace;

	private float dataMin = Float.MAX_VALUE;
	private float dataMax = -Float.MAX_VALUE;

	private int[] traceTypes = new int[100];
	private int typeCount = 0;

    protected double progressValue = 0.0;
	protected String progressDescription = "";
	private boolean canceled = false;

	private int bytesPerTrace;

	private int dataCnt[] = new int[255];
	private int ttlHistogramSamples = 0;
	protected int analysisResult = 0;

	protected StsAnalyzeTraces()
	{
	}

	private void initialise()
	{
		traceTypes = null;
		traceTypes = new int[100];
		typeCount = 0;

		ttlHistogramSamples = 0;
		for(int i = 0; i < 255; i++)
		{
			dataCnt[i] = 0;
		}

		data = new byte[nSampleBytesPerTrace];
		traceHeader = new byte[traceHeaderSize];
		data = new byte[nSampleBytesPerTrace + traceHeaderSize];
		ioffs = new long[nTracesToScan];

		long ioff = fileHeaderSize;
//		ioff = ioff + traceHeaderSize;

		bytesPerTrace = nSampleBytesPerTrace + traceHeaderSize;

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

	public static StsAnalyzeTraces constructor(RandomAccessFile randomAccessSegyFile,
											   int fileHeaderSize,
											   int nTracesToScan,
											   int nTotalTraces,
											   int nSamplesPerTrace,
											   int bytesPerSample,
											   int traceHeaderSize,
											   boolean isLittleEndian,
											   int sampleFormat,
											   StsSEGYFormatRec dataScaleRec,
											   StsSEGYFormatRec traceTypeRec,
                                               String name)
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
		analyze.nSamplesPerTrace = nSamplesPerTrace;
		analyze.randomAccessSegyFile = randomAccessSegyFile;
		analyze.fileHeaderSize = fileHeaderSize;
		analyze.nTracesToScan = nTracesToScan;
		analyze.nTotalTraces = nTotalTraces;
		analyze.nSampleBytesPerTrace = nSamplesPerTrace * bytesPerSample;
		analyze.bytesPerSample = bytesPerSample;
		analyze.traceHeaderSize = traceHeaderSize;
		analyze.isLittleEndian = isLittleEndian;
		analyze.dataScaleRec = dataScaleRec;
		analyze.traceTypeRec = traceTypeRec;
        analyze.name = name;
        analyze.initializeDataRange();
        return analyze;
	}

    private void initializeDataRange()
    {
        dataMin = Float.MAX_VALUE;
        dataMax = -Float.MAX_VALUE;
    }

    public void analyze(StsProgressPanel progressPanel)
	{
		try
		{
			initialise();
			preOuterLoop();

            int progressIncrement = Math.max(nTracesToScan/100, 1);
            int incrementCounter = 0;
            int progressValue = 0;
            progressPanel.initialize(nTracesToScan);
            progressPanel.setDescriptionAndLevel("Analyzing traces for volume " + name, StsProgressBar.INFO);
            for (int i = 0; i < nTracesToScan; i++)
			{
				randomAccessSegyFile.seek(ioffs[i]);
				randomAccessSegyFile.read(traceHeader);
				randomAccessSegyFile.read(data);

				accumulateTypes();
				preInnerLoop();

				if (progressPanel.isCanceled() || analysisResult != ANALYSIS_OK)//!fmtOkay)
				{
					//fmtOkay = false;
					return;
				}

				for (int j = indexStart; j < (nSamplesPerTrace * indexStep) + indexStart; j += indexStep)
				{
					sample = getSample(j);
					if ((sample == -Float.MAX_VALUE) || (sample == Float.MAX_VALUE))
						sample = 0.0f;
					accumulateHistogram(sample);
					if (sample < dataMin)
					{
						dataMin = sample;
					}
					if (sample > dataMax)
					{
						dataMax = sample;
					}
					total = total + sample;
					if (canceled || analysisResult != ANALYSIS_OK)//!fmtOkay)
					{
						//fmtOkay = false;
						return;
					}
				}
                incrementCounter++;
                if(incrementCounter == progressIncrement)
                {
                    progressValue += progressIncrement;
                    progressPanel.setValue(progressValue);
                    incrementCounter = 0;
                }
			}
			traceTypes = (int[])StsMath.trimArray(traceTypes, typeCount);

			progressPanel.finished();
		}
		catch(IOException ioe)
		{
			analysisResult = this.ANALYSIS_IO_EXCEPTION;
			//fmtOkay = false;
			StsException.systemError("StsAnalyzeTraces.analyze() failed to read randomAccessFile");
		}
		catch(Exception e)
		{
			analysisResult = this.ANALYSIS_UNKNOWN_EXCEPTION;
			//fmtOkay = false;
			StsException.outputException("StsAnalyzeTraces.analyze() failed.", e, StsException.WARNING);
		}
	}

    public void analyze()
	{
		try
		{
			initialise();
			preOuterLoop();

			stepIncrement = 100.0 / (nTracesToScan);
			progressValue = 0.0;
			progressDescription = "";
			for (int i = 0; i < nTracesToScan; i++)
			{
				randomAccessSegyFile.seek(ioffs[i]);
				randomAccessSegyFile.read(traceHeader);
				randomAccessSegyFile.read(data);

				accumulateTypes();
				preInnerLoop();

				if (canceled || analysisResult != ANALYSIS_OK)//!fmtOkay)
				{
					//fmtOkay = false;
					return;
				}

				for (int j = indexStart; j < (nSamplesPerTrace * indexStep) + indexStart; j += indexStep)
				{
					sample = getSample(j);
					if ((sample == -Float.MAX_VALUE) || (sample == Float.MAX_VALUE))
						sample = 0.0f;
					accumulateHistogram(sample);
					if (sample < dataMin)
					{
						dataMin = sample;
					}
					if (sample > dataMax)
					{
						dataMax = sample;
					}
					total = total + sample;
					if (canceled || analysisResult != ANALYSIS_OK)//!fmtOkay)
					{
						//fmtOkay = false;
						return;
					}
				}
				progressValue += stepIncrement;
				progressDescription = " Scanned Trace " + i + " of " + nTracesToScan + " to be scanned.";

				// Truncate the type area
			}
			traceTypes = (int[])StsMath.trimArray(traceTypes, typeCount);

			progressDescription = "";
			progressValue = 100.0;
		}
		catch(IOException ioe)
		{
			analysisResult = this.ANALYSIS_IO_EXCEPTION;
			//fmtOkay = false;
			StsException.systemError("StsAnalyzeTraces.analyze() failed to read randomAccessFile");
		}
		catch(Exception e)
		{
			analysisResult = this.ANALYSIS_UNKNOWN_EXCEPTION;
			//fmtOkay = false;
			StsException.outputException("StsAnalyzeTraces.analyze() failed.", e, StsException.WARNING);
		}
	}

    public void resetProgressBar()
    {
        progressDescription = "";
        progressValue = 0.0;
    }
	private void accumulateHistogram(float value)
	{
		float scaledFloat = 254 * (value - dataMin) / (dataMax - dataMin);
		int index = (int)StsMath.minMax(Math.round(scaledFloat), 0, 254);
		dataCnt[index] += 1;
		ttlHistogramSamples++;
	}

	public float[] calculateHistogram()
	{
		float[] dataHist = new float[255];
		for(int i = 0; i < 255; i++)
			dataHist[i] = (float)((float)dataCnt[i] / (float)ttlHistogramSamples) * 100.0f;
		return dataHist;
	}

	public int[] getHistogramData()
	{
		return dataCnt;
	}

	public int getTtlHistogramSamples()
	{
		return ttlHistogramSamples;
	}

	public float getDataMin()
	{
		return dataMin;
	}

	public float getDataMax()
	{
		return dataMax;
	}

	public void setDataMin(float dataMin)
	{
		this.dataMin = dataMin;
	}

	public void setDataMax(float dataMax)
	{
		this.dataMax = dataMax;
	}

	public int[] getTraceTypes()
	{
		return traceTypes;
	}

	protected void preOuterLoop()
	{
	}

	protected void preInnerLoop()
	{

	}

	protected abstract float getSample(int index);

	public double getTotal()
	{
		return total;
	}

	public double getAverage()
	{
		return(total / (double)(nTracesToScan * nSamplesPerTrace));
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
		indexStart = traceHeaderSize;
		indexStep = 1;
	}

	protected float getSample(int index)
	{
		return(float)data[index];
	}
}

class AnalyzeInt2 extends StsAnalyzeTraces
{
	AnalyzeInt2()
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 2;
	}

	protected float getSample(int index)
	{
		return(float)StsMath.shortValue(index, data, isLittleEndian);
	}
}

class AnalyzeInt4 extends StsAnalyzeTraces
{
	AnalyzeInt4()
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 4;
	}

	protected float getSample(int index)
	{
		return(float)StsMath.intValue(index, data, isLittleEndian);
	}
}

class AnalyzeIBMFloat extends StsAnalyzeTraces
{
	float[] fdata;

	AnalyzeIBMFloat()
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 1;
	}

	protected void preOuterLoop()
	{
		fdata = new float[nSamplesPerTrace];
	}

	protected void preInnerLoop()
	{
		if(!StsMath.verifyIBMFloatType(data, isLittleEndian))
		{
			progressDescription = " Error with IBM Float format. Try IEEE Float.";
			analysisResult = this.ANALYSIS_INVALID_FORMAT;
			//fmtOkay = false;
			return;
		}
		StsMath.convertIBMFloatBytes(data, fdata, nSamplesPerTrace, isLittleEndian);
	}

	protected float getSample(int index)
	{
		return(float)fdata[index];
	}
}

class AnalyzeIEEEFloat extends StsAnalyzeTraces
{
	float[] fdata;

	AnalyzeIEEEFloat()
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 4;
	}

	protected void preOuterLoop()
	{
		fdata = new float[nSamplesPerTrace];
	}

	protected float getSample(int index)
	{
		return Float.intBitsToFloat(StsMath.intValue(index, data, isLittleEndian));
	}
}

class AnalyzeFloat8 extends StsAnalyzeTraces
{
	int dataScale;

	AnalyzeFloat8(int traceHeaderSize)
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 1;
	}

	protected void preOuterLoop()
	{
//		int nSampleBytesPerTrace = nSamplesPerTrace * bytesPerSample;
//		data = new byte[nSampleBytesPerTrace + traceHeaderSize];
	}

	protected void preInnerLoop()
	{
		dataScale = (int)dataScaleRec.getHdrValue(traceHeader, isLittleEndian);
	}

	protected float getSample(int index)
	{
		return(float)data[index] / (float)(Math.pow(2.0f, (float)dataScale));
	}
}

class AnalyzeFloat16 extends StsAnalyzeTraces
{
	int dataScale;

	AnalyzeFloat16(int traceHeaderSize)
	{
		super();
		indexStart = traceHeaderSize;
		indexStep = 2;
	}

	protected void preOuterLoop()
	{
//		int nSampleBytesPerTrace = nSamplesPerTrace * bytesPerSample;
//		data = new byte[nSampleBytesPerTrace + traceHeaderSize];
	}

	protected void preInnerLoop()
	{
		dataScale = (int)dataScaleRec.getHdrValue(traceHeader, isLittleEndian);
	}

	protected float getSample(int index)
	{
		return(float)StsMath.shortValue(index, data, isLittleEndian) /
			(float)(Math.pow(2.0f, (float)dataScale));
	}
}
