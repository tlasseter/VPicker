package com.Sts.SeismicAttributes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSeismicVelocityConstructor extends StsSeismicVolumeConstructor
{
	StsSeismicVelocityModel velocityModel;
	StsSeismicVolume inputVelocityVolume;
	float[] inputTraceFloats;
//    byte[] inputTraceBytes;
    byte[] outputTraceBytes;
	float[] outputTraceFloats;
	int nSurfaces = 0;
	StsVelocityGrid[] velocityGrids;
	StsModelSurface[] timeSurfaces;
//    double[][][] surfaceCorrectionFactors;
	float[][][] surfaceTimes;
	float[][][] surfaceDepths;
	double[][][] surfaceVelocities;
	float[] oneWayConstantIntervalVelocities;
	float oneWayConstantBottomIntervalVelocity;
	boolean hasConstantLayerIntervalVelocities = false;
	float timeMin, timeMax, timeInc;
	float timeDatum, depthDatum;
	int nSlices;
	/** row offset of volume in seismicVelocityModel grid. */
	public int volumeRowOffset;
	/** col offset of volume in seismicVelocityModel grid. */
	public int volumeColOffset;
	int nModelRows, nModelCols;
	double scaleFactor, scaleOffset;
	/** factor multiplying velocities in input velocity volume. Since seismic data and picked horizons are two-way,
	 *  if input velocity volume is one-way, factor is 0.5 to reduce one-way to two-way velocities.
	 */
	double scaleMultiplier = 1.0;

    /** assume the max reasonable two-way velocity is 10 m/msec; convert this to project units velocity for comparision */
    float maxReasonableVelocity;

    private float depthMax = 0.0f;
//	static public float depthInc = 0.0f;
	static public final String AVG_VELOCITY = "avg_velocity";
    static public final String INT_VELOCITY = "intv_velocity";
	StsTimer timer = null;
	boolean runTimer = false;

	private StsSeismicVelocityConstructor(StsModel model, StsSeismicVelocityModel velocityModel, StsProgressPanel panel) throws StsException
	{
		this.velocityModel = velocityModel;
        this.model = model;
        this.panel = panel;
		nSurfaces = velocityModel.getNSurfaces();
		velocityGrids = velocityModel.getVelocityGrids();
		timeSurfaces = velocityModel.getTimeSurfaces();
		inputVelocityVolume = velocityModel.getInputVelocityVolume();
		timeDatum = velocityModel.timeDatum;
		depthDatum = velocityModel.depthDatum;
		scaleMultiplier = velocityModel.scaleMultiplier;
        maxReasonableVelocity = model.getProject().convertVelocityToProjectUnits(10.0f, StsParameters.VEL_M_PER_MSEC);
        depthMax = 0.0f;
		hasConstantLayerIntervalVelocities = velocityModel.hasConstantLayerIntervalVelocities;
        displayProcessingBlock = false;
        if (inputVelocityVolume != null)
		{
			if(!constructFromInputVelocityVolume(model, velocityModel, AVG_VELOCITY, StsParameters.SAMPLE_TYPE_VEL_AVG)) throw new StsException(StsException.WARNING,  "Failed to construct seismic velocity volume.");
		}
		else
		{
			if(!constructVelocityVolume(model, velocityModel, AVG_VELOCITY,StsParameters.SAMPLE_TYPE_VEL_AVG)) throw new StsException(StsException.WARNING,  "Failed to construct seismic velocity volume.");
		}
    }
	private StsSeismicVelocityConstructor(StsModel model, StsSeismicVelocityModel velocityModel, StsProgressPanel panel, String mode, byte type) throws StsException
	{
		this.velocityModel = velocityModel;
		this.model = model;
		this.panel = panel;
		nSurfaces = velocityModel.getNSurfaces();
		velocityGrids = velocityModel.getVelocityGrids();
		timeSurfaces = velocityModel.getTimeSurfaces();
		inputVelocityVolume = velocityModel.getInputVelocityVolume();
		timeDatum = velocityModel.timeDatum;
		depthDatum = velocityModel.depthDatum;
		scaleMultiplier = velocityModel.scaleMultiplier;
		maxReasonableVelocity = model.getProject().convertVelocityToProjectUnits(10.0f, StsParameters.VEL_M_PER_MSEC);
		depthMax = 0.0f;
		hasConstantLayerIntervalVelocities = velocityModel.hasConstantLayerIntervalVelocities;
		displayProcessingBlock = false;
		if (inputVelocityVolume != null)
		{
			if(!constructFromInputVelocityVolume(model, velocityModel, mode, type)) throw new StsException(StsException.WARNING,  "Failed to construct seismic velocity volume.");
		}
		else
		{
			if(!constructVelocityVolume(model, velocityModel, mode, type)) throw new StsException(StsException.WARNING,  "Failed to construct seismic velocity volume.");
		}
    }
	static public StsSeismicVolume constructor(StsModel model, StsSeismicVelocityModel velocityModel, StsProgressPanel panel)
	{
		try
		{
			StsSeismicVolume outputVolume = StsSeismicVolumeConstructor.getExistingVolume(model, velocityModel.stemname + "." + AVG_VELOCITY);
			if (outputVolume != null)
			{
				outputVolume.delete();
			}

			StsSeismicVelocityConstructor seismicVelocityConstructor = new StsSeismicVelocityConstructor(model, velocityModel, panel, AVG_VELOCITY, StsParameters.SAMPLE_TYPE_VEL_AVG);
			outputVolume = seismicVelocityConstructor.getVolume();
			if (!outputVolume.initialize(model))
			{
				outputVolume.delete();
				return null;
			}
			return outputVolume;
		}
		catch (Exception e)
		{
			StsMessage.printMessage("StsSeismicVelocityConstructor.constructor() failed.");
			e.printStackTrace();
			return null;
		}
	}
	static public StsSeismicVolume intervalconstructor(StsModel model, StsSeismicVelocityModel velocityModel, StsProgressPanel panel)
		{
			try
			{
				StsSeismicVolume outputVolume = StsSeismicVolumeConstructor.getExistingVolume(model, velocityModel.stemname + "." + INT_VELOCITY);
				if (outputVolume != null)
				{
					outputVolume.delete();
				}

				StsSeismicVelocityConstructor seismicVelocityConstructor = new StsSeismicVelocityConstructor(model, velocityModel, panel, INT_VELOCITY, StsParameters.SAMPLE_TYPE_VEL_INTERVAL);
				outputVolume = seismicVelocityConstructor.getVolume();
				if (!outputVolume.initialize(model))
				{
					outputVolume.delete();
					return null;
				}
				return outputVolume;
			}
			catch (Exception e)
			{
				StsMessage.printMessage("StsSeismicVelocityConstructor.constructor() failed.");
				e.printStackTrace();
				return null;
			}
	}
	private boolean constructFromInputVelocityVolume(StsModel model, StsSeismicVelocityModel velocityModel, String mode, byte type)
	{
//		StsSeismicVolume inputVelocityVolume = velocityModel.getInputVelocityVolume();
		initialize(inputVelocityVolume, mode);
		timeMin = inputVelocityVolume.zMin;
		timeMax = inputVelocityVolume.zMax;
		timeInc = inputVelocityVolume.zInc;
//        surfaceCorrectionFactors = new double[nSurfaces][][];
		surfaceVelocities = new double[nSurfaces][][];
		surfaceTimes = new float[nSurfaces][][];
		surfaceDepths = new float[nSurfaces][][];
		for (int n = 0; n < nSurfaces; n++)
		{
//            surfaceCorrectionFactors[n] = velocityGrids[n].getErrorCorrections();
			surfaceVelocities[n] = velocityGrids[n].getVelocities();
			surfaceTimes[n] = timeSurfaces[n].getPointsZ();
			surfaceDepths[n] = timeSurfaces[n].getAdjPointsZ();
		}
		if (runTimer)
		{
			timer = new StsTimer();
		}
		initializeScalingFromVolume(inputVelocityVolume);

        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, velocityModel, dataMin, dataMax, true, false, inputVelocityVolume.stemname, mode, "rw");
		outputVolume.setZDomain(StsParameters.TD_TIME);
		outputVolume.volumeType = type;
		outputVolume.setType(type);
		outputVolume.vertUnitsString = velocityModel.vertUnitsString;
		outputVolume.oneOrTwoWayVelocity = StsParameters.TWO_WAY_VELOCITY;
		outputVolume.velocityUnits = model.getProject().getVelocityUnits();
		// jbw rewrite header with added info
		outputVolume.writeHeaderFile();
//        initializeVolume(model, inputVelocityVolume, dataMin, dataMax, StsParameters.VELOCITY, StsProject.TD_TIME);
		nSlices = inputVelocityVolume.nSlices;
        inputDataIsFloat = inputVelocityVolume.isDataFloat;
//        if(inputDataIsFloat)
        inputTraceFloats = new float[nSlices];
//        else
//            inputTraceBytes = new byte[nInputSlices];

        outputTraceFloats = new float[nSlices];

		nModelRows = velocityModel.nRows;
		nModelCols = velocityModel.nCols;
        createOutputVolume(mode);
//		attributeVolume.setZDomain(StsParameters.TD_TIME);
        return true;
    }

	private void initializeScalingFromVolume(StsSeismicVolume velocityVolume)
	{
		dataMin = (float) scaleMultiplier * velocityVolume.dataMin;
		dataMax = (float) scaleMultiplier * velocityVolume.dataMax;
		scaleFactor = 254 / (dataMax - dataMin);
		scaleOffset = -dataMin * scaleFactor;
	}

	private void initializeScalingFromIntervalVelocities(StsSeismicVelocityModel velocityVolume)
	{
        StsVelocityGrid[] velocityGrids = velocityVolume.intervalVelocityGridArray;
        dataMin = velocityGrids[0].velMin;
        dataMax = velocityGrids[0].velMax;
        for(int n = 1; n < velocityGrids.length; n++)
        {
            dataMin = Math.min(dataMin, velocityGrids[n].velMin);
            dataMax = Math.max(dataMax, velocityGrids[n].velMax);
        }
		scaleFactor = 254 / (dataMax - dataMin);
		scaleOffset = -dataMin * scaleFactor;
	}

    public void initializeVolumeOutput()
    {
        outputVolume.setDataRange(dataMin, dataMax);
        velocityModel.setDataRange(dataMin, dataMax);
        super.initializeVolumeOutput();
    }

    private boolean constructVelocityVolume(StsModel model, StsSeismicVelocityModel velocityModel, String mode, byte type)
	{
		initialize(mode);
		timeMin = velocityModel.timeMin;
		timeMax = velocityModel.timeMax;
		timeInc = velocityModel.timeInc;
//		timeMin = velocityModel.zMin;
//		timeMax = velocityModel.zMax;
//		timeInc = velocityModel.zInc;

		surfaceVelocities = new double[nSurfaces][][];
		surfaceTimes = new float[nSurfaces][][];
		surfaceDepths = new float[nSurfaces][][];
		for (int n = 0; n < nSurfaces; n++)
		{
			surfaceVelocities[n] = velocityGrids[n].getVelocities();
			surfaceTimes[n] = timeSurfaces[n].getPointsZ();
			surfaceDepths[n] = timeSurfaces[n].getAdjPointsZ();
		}
		oneWayConstantIntervalVelocities = velocityModel.getOneWayConstantIntervalVelocities();
		oneWayConstantBottomIntervalVelocity = velocityModel.getOneWayConstantBottomIntervalVelocity();
//		initializeScalingFromIntervalVelocities(velocityModel);


        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, velocityModel, dataMin, dataMax, true, false, velocityModel.stemname, mode, "rw");
//        initializeVolume(model, velocityModel, dataMin, dataMax, StsParameters.VELOCITY, StsProject.TD_TIME);
		outputVolume.setZDomain(StsParameters.TD_TIME);
		outputVolume.volumeType = type;
		outputVolume.setType(type);
		outputVolume.vertUnitsString = velocityModel.vertUnitsString;
		outputVolume.oneOrTwoWayVelocity = StsParameters.TWO_WAY_VELOCITY;
		outputVolume.velocityUnits = model.getProject().getVelocityUnits();
		outputVolume.isDataFloat = true;

		// jbw rewrite header with added info
		outputVolume.writeHeaderFile();

		nSlices = velocityModel.nSlices;

        outputTraceFloats = new float[nSlices];
		/*
			if (panel != null)
			{
				panel.initializeProgressBar(nInputRows * nInputCols);
			}
		 */
		volumeRowOffset = 0;
		volumeColOffset = 0;
		nModelRows = velocityModel.nRows;
		nModelCols = velocityModel.nCols;
        createOutputVolume(mode);
        return true;
//       initializeVolume(model);
	}

    public boolean isOutputDataFloat()
    {
		return true;
    }

    /*
		private void initializeVolume(StsModel model)
		{
			String stemname = "intervalVelocity";

			attributeVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, stemname);
			if (attributeVolume != null) attributeVolume.delete();

			attributeVolume = new StsSeismicVolume(false);

	 attributeVolume.stsDirectory = model.getProject().getDataFullDirString();
			attributeVolume.setName(stemname);
			attributeVolume.setFilenames();
			attributeVolume.deleteExistingFiles();

			attributeVolume.setIsDataFloat(isDataFloat);

			attributeVolume.initializeToBoundingBox(model.getGridDefinition());
	 attributeVolume.xMax = attributeVolume.xMin + attributeVolume.xInc*(nInputCols-1);
	 attributeVolume.yMax = attributeVolume.yMin + attributeVolume.yInc*(nInputRows-1);
			timeMax = model.project.getTimeMax();
			timeMin = model.project.getTimeMin();
			timeInc = model.project.getTimeInc();
			int dtInc = StsMath.below((velocityModel.t0 - timeMin)/timeInc);
			timeMin += dtInc*timeInc;
			nInputSlices = Math.round((timeMax - timeMin)/timeInc) + 1;

			// PostStack3d is native in time units, so set z range values to time.
			attributeVolume.zMin = timeMin;
			attributeVolume.zMax = timeMax;
			attributeVolume.zInc = timeInc;

			float[] dataRange = velocityModel.computeVelocityRange();
			attributeVolume.setDataMin(dataRange[0]);
			attributeVolume.setDataMax(dataRange[1]);

			attributeVolume.writeHeaderFile();
			attributeVolume.initializeSpectrum();
			attributeVolume.addToModel();
			attributeVolume.refreshObjectPanel();

			nInputRows = attributeVolume.getNRows();
			nInputCols = attributeVolume.getNCols();
			nInputSlices = attributeVolume.getNSlices();

			if (isDataFloat) attributeVolume.createMappedRowBuffer("rw");

//        readWriteVolume();

			 createVolume();
			 attributeVolume.classInitialize(model);
		}
	 */

    public void initializeVolumeInput()
    {
        initializeGridRange();

        nInputBlockLastRow = -1;
        nInputRowsPerBlock = 1;
        nInputBlocks = nInputRows;
        nOutputSamplesPerInputBlock = memoryAllocation.nOutputSamplesPerInputBlock;

        dataMin = StsParameters.largeFloat;
        dataMax = -StsParameters.largeFloat;
        traceFloats = new float[nInputSlices];

        panel.appendLine("Creating volume: " + outputVolume.getName());
        panel.appendLine("Volume size: nInputRows " + nInputRows + " nInputCols " + nInputCols + " nSamples " + nInputSlices);
        panel.appendLine("Input processing will be in " + nInputBlocks + " blocks with " + nInputRowsPerBlock + " rows per block.");
        panel.initializeIntervals(2, 100);
        panel.setIntervalCount(nInputBlocks);

        outputFloatBuffer = outputVolume.createMappedFloatRowBuffer();
        outputPosition = 0;
        if(inputVolumes == null) return;
        int nInputVolumes = inputVolumes.length;
        this.inputIsFloat = inputVolumesAreFloat();
        if (nInputVolumes > 0)
        {
            inputBuffers = new StsMappedBuffer[nInputVolumes];
            for (int n = 0; n < nInputVolumes; n++)
                inputBuffers[n] = inputVolumes[n].createMappedRowBuffer("r");
        }
        inputPosition = 0;
    }

    protected void initializeGridRange()
    {
        nInputRows = memoryAllocation.nOutputRows;
        nInputCols = memoryAllocation.nOutputCols;
        nInputSlices = memoryAllocation.nOutputSlices;

        nOutputRows = nInputRows;
        nOutputCols = nInputCols;
        nOutputSlices = nInputSlices;
//        inputRowInc = 1;
//        inputColInc = 1;
//        inputSliceInc = 1;
    }

    public boolean initializeBlockInput(int nBlock)
    {
        nInputBlockFirstRow = nInputBlockLastRow + 1;
        nInputBlockLastRow += nInputRowsPerBlock;

        if (nInputBlockLastRow > nInputRows - 1)
        {
            nInputBlockLastRow = nInputRows - 1;
            nInputRowsPerBlock = nInputBlockLastRow - nInputBlockFirstRow + 1;
            nInputBlockTraces = nInputRowsPerBlock * nInputCols;
            nOutputSamplesPerInputBlock = nInputRowsPerBlock * nInputCols * nInputSlices;
        }

        if (debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " nInputSamplesPerInputBlock: " + nOutputSamplesPerInputBlock);
        if(!outputFloatBuffer.map(outputPosition, nOutputSamplesPerInputBlock)) return false;
        outputPosition += nOutputSamplesPerInputBlock;

        return true;
    }

	public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers) { return doProcessInputBlock(nBlock, inputBuffers, null);}

	public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers, String mode)
	{
		if (runTimer) timer.start();
        try
        {
            if (inputVelocityVolume != null)
			{
				if ((mode != null) && mode.equals(INT_VELOCITY))
					processInputVelocityIntervalVolumeBlock(nBlock, mode);
				else
				   processInputVelocityVolumeBlock(nBlock, mode);
			}
            else
			{
				if ((mode != null) && mode.equals(INT_VELOCITY))
					processIntervalVelocityIntervalVolumeBlock(nBlock, mode);
				else
				    processIntervalVelocityVolumeBlock(nBlock, mode);
			}
            if (runTimer) timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "doProcessBlock", e);
            return false;
        }
    }
    final private void processInputVelocityVolumeBlock(int nBlock, String mode)
	{
		int nTrace = 0;

        float outputY = outputVolume.getYCoor(nInputBlockFirstRow);
        float outputYInc = outputVolume.yInc;
        for (int outputRow = nInputBlockFirstRow; outputRow <= nInputBlockLastRow; outputRow++, outputY += outputYInc)
		{
            float outputX = outputVolume.getXCoor(0);
            float outputXInc = outputVolume.xInc;
            for (int outputCol = 0; outputCol < nInputCols; outputCol++, nTrace++, outputX += outputXInc)
			{
				if (isCanceled())
				{
					return;
				}
				if (panel != null && (++nTracesDone % 1000) == 0)
				{
                    panel.appendLine("Completed processing " + nTracesDone + " traces.");
				}
                if(!getInterpolatedInputTrace(outputX, outputY, inputTraceFloats))
                    processIntervalVelocityVolumeTrace(outputRow, outputCol, nTrace, outputFloatBuffer, outputVolume);
                else
                {
					if(inputVelocityVolume.getType() == StsParameters.SAMPLE_TYPE_VEL_INSTANT)
					{  // jbw input interval velocity volume can give the average volume better structure
						double scaleMultiplier2 = scaleMultiplier; // units fix /(39.37d/12.d);
						scaleConvertIntervalVelocity(inputTraceFloats, scaleMultiplier2, timeInc, inputVelocityVolume.getType());
						processVelocityIntervalAverageVolumeTrace(outputRow, outputCol, nTrace, outputFloatBuffer, outputVolume, scaleMultiplier2, timeInc);
					}
					else
					{
						scaleConvertVelocity(inputTraceFloats, scaleMultiplier, timeInc, inputVelocityVolume.getType());
						processVelocityVolumeTrace(outputRow, outputCol, nTrace, outputFloatBuffer, outputVolume);
					}
                }
                if (debug && nTrace % 1000 == 0)
				{
					System.out.println("    processed trace: " + nTrace);
				}
			}
        }
    }
	final private void processInputVelocityIntervalVolumeBlock(int nBlock, String mode)
	{
	int nTrace = 0;

	float outputY = outputVolume.getYCoor(nInputBlockFirstRow);
	float outputYInc = outputVolume.yInc;
	for (int outputRow = nInputBlockFirstRow; outputRow <= nInputBlockLastRow; outputRow++, outputY += outputYInc)
	{
		float outputX = outputVolume.getXCoor(0);
		float outputXInc = outputVolume.xInc;
		for (int outputCol = 0; outputCol < nInputCols; outputCol++, nTrace++, outputX += outputXInc)
		{
			if (isCanceled())
			{
				return;
			}
			if (panel != null && (++nTracesDone % 1000) == 0)
			{
				panel.appendLine("Completed processing " + nTracesDone + " traces.");
			}
			if(!getInterpolatedInputTrace(outputX, outputY, inputTraceFloats))
				processIntervalVelocityIntervalVolumeTrace(outputRow, outputCol, nTrace, outputFloatBuffer, outputVolume);
			else
			{
				//scaleMultiplier /= (39.37d/12.d);
				double scaleMultiplier2 = scaleMultiplier ; // units fix /(39.37d/12.d);
				scaleConvertIntervalVelocity(inputTraceFloats, scaleMultiplier2, timeInc, inputVelocityVolume.getType());
				processVelocityIntervalVolumeTrace(outputRow, outputCol, nTrace, outputFloatBuffer, outputVolume,  scaleMultiplier2, timeInc);
			}
			if (debug && nTrace % 1000 == 0)
			{
				System.out.println("    processed trace: " + nTrace);
			}
		}
	}
}

/*
    final private void processInputVelocityVolumeBlock(int nBlock)
	{
		StsMappedBuffer inputBuffer = inputBuffers[0];
		inputBuffer.position(0);
		int nTrace = 0;
 //       FloatBuffer floatBuffer = outputRowFloatBuffer.asFloatBuffer();
        for (int row = nInputBlockFirstRow; row <= nInputBlockLastRow; row++)
		{
			for (int col = 0; col < nInputCols; col++, nTrace++)
			{
				if (isCanceled())
				{
					return;
				}
				if (panel != null && (++nTracesDone % 1000) == 0)
				{
                    panel.appendLine("Completed processing " + nTracesDone + " traces.");
				}
				inputBuffer.get(inputTraceFloats);
                scaleConvertVelocity(inputTraceFloats, scaleMultiplier, timeInc, inputVelocityVolume.getType());
                processVelocityVolumeTrace(row, col, nTrace, outputFloatBuffer, outputVolume);
				if (debug && nTrace % 1000 == 0)
				{
					System.out.println("    processed trace: " + nTrace);
				}
			}
//            panel.incrementCount();
        }
//       panel.appendLine("Model Construction Complete");
	}
*/
    private boolean getInterpolatedInputTrace(float outputX, float outputY, float[] inputTraceFloats)
    {
		//jbw
        float rowF = inputVelocityVolume.getRowCoor(outputY);
        float colF = inputVelocityVolume.getColCoor(outputX);

        return velocityModel.getVelocityTraceValues(inputVelocityVolume, rowF, colF, inputTraceFloats);
    }

    private void scaleConvertVelocity(float[] inputTraceFloats, double scaleMultiplier, float timeInc, byte velocityType)
    {
        if(velocityType == StsParameters.SAMPLE_TYPE_VEL_INSTANT)
            StsMath.convertInstantToAvgVelocity(inputTraceFloats, scaleMultiplier, timeInc);
        else
            StsMath.scale(inputTraceFloats, (float)scaleMultiplier);
    }

	private void scaleConvertIntervalVelocity(float[] inputTraceFloats, double scaleMultiplier, float timeInc, byte velocityType)
	{
		 if(velocityType == StsParameters.SAMPLE_TYPE_VEL_AVG)
	        StsMath.convertAvgToInstantVelocity(inputTraceFloats, scaleMultiplier, timeInc);
		 else
			 StsMath.scale(inputTraceFloats, (float)scaleMultiplier);
	}

    final private void processIntervalVelocityVolumeBlock(int nBlock, String mode)
	{
		int nTrace = 0;
 //       FloatBuffer floatBuffer = outputRowFloatBuffer.asFloatBuffer();
        for (int row = nInputBlockFirstRow; row <= nInputBlockLastRow; row++)
		{
			for (int col = 0; col < nInputCols; col++, nTrace++)
			{
				if (isCanceled())
				{
					return;
				}
				if (panel != null && (++nTracesDone % 1000) == 0)
				{
                    panel.appendLine("Completed processing " + nTracesDone + " traces.");
				}
				processIntervalVelocityVolumeTrace(row, col, nTrace, outputFloatBuffer, outputVolume);
				if (debug && nTrace % 1000 == 0)
				{
					System.out.println("    processed trace: " + nTrace);
				}
			}
        }
//        panel.appendLine("Model Construction Complete");
	}

	final private void processIntervalVelocityIntervalVolumeBlock(int nBlock, String mode)
	{
		int nTrace = 0;
 //       FloatBuffer floatBuffer = outputRowFloatBuffer.asFloatBuffer();
		for (int row = nInputBlockFirstRow; row <= nInputBlockLastRow; row++)
		{
			for (int col = 0; col < nInputCols; col++, nTrace++)
			{
				if (isCanceled())
				{
					return;
				}
				if (panel != null && (++nTracesDone % 1000) == 0)
				{
					panel.appendLine("Completed processing " + nTracesDone + " traces.");
				}
				processIntervalVelocityIntervalVolumeTrace(row, col, nTrace, outputFloatBuffer, outputVolume);
				if (debug && nTrace % 1000 == 0)
				{
					System.out.println("    processed trace: " + nTrace);
				}
			}
		}
//        panel.appendLine("Model Construction Complete");
	}

	final private boolean processVelocityVolumeTrace(int volumeRow, int volumeCol, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume outputVolume)
	{
		int n = -1, s = -1;
		boolean debug = false;
		try
		{
			int modelRow = volumeRow + volumeRowOffset;
			int modelCol = volumeCol + volumeColOffset;
			modelRow = StsMath.minMax(modelRow, 0, nModelRows - 1);
			modelCol = StsMath.minMax(modelCol, 0, nModelCols - 1);

			if (StsSeismicVelocityModel.debug == true && volumeRow == StsSeismicVelocityModel.debugVolumeRow && volumeCol == StsSeismicVelocityModel.debugVolumeCol)
			{
				debug = false;
				System.out.println("processVelocityVolumeTrace mainDebug. file " + rowFloatBuffer.pathname + " volume row " + volumeRow + " volume col " + volumeCol);
				if (rowFloatBuffer != null)
				{
					long position = rowFloatBuffer.getBufferPosition();
					long correctPosition = volumeRow * nInputCols * nSlices + volumeCol * nSlices;
					System.out.println("    start write position " + position + " correct position " + correctPosition);
				}
			}

//            double[] traceErrorCorrections = getTraceErrorCorrections(offsetRow, offsetCol);
			float[] traceTimeValues = getTraceTimeSurfaces(modelRow, modelCol);
			float[] traceDepthValues = getTraceDepthSurfaces(modelRow, modelCol);
//            double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);

			float botTime = timeDatum;
			double botError = 1.0f;
			float botDepth;
			double botAvgVolumeVelocity;
			double botAvgSurfaceVelocity;
			double errorCorrection = 1.0;

			int botSlice = -1;

			double topError;
			int topSlice;

			for (s = 0; s < nSurfaces; s++)
			{
				topError = botError;
				botTime = traceTimeValues[s];
				botDepth = traceDepthValues[s];
				float botSliceF = outputVolume.getBoundedSliceCoor(botTime);
				botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botSliceF);
				botAvgSurfaceVelocity = (botDepth - depthDatum) / (botTime - timeDatum);
				botError = botAvgSurfaceVelocity / botAvgVolumeVelocity;
                if(debug) System.out.println("Surface: " + timeSurfaces[s].getName() + " surfaceAvgVelocity " + botAvgSurfaceVelocity +
                                            " botAvgVolumeVelocity " + botAvgVolumeVelocity + " error correction: " + botError);
                topSlice = botSlice + 1;
				botSlice = (int) botSliceF;

				for (n = topSlice; n <= botSlice; n++)
				{
                    float f = 0.0f;
                    if(botSlice != topSlice) f = ( (float) (n - topSlice)) / (botSlice - topSlice);
                    errorCorrection = topError + f * (botError - topError);
					outputTraceFloats[n] = (float)errorCorrection*inputTraceFloats[n];
					if (dataMin > outputTraceFloats[n]);
						dataMin = outputTraceFloats[n];
					if (dataMax < outputTraceFloats[n]);
						dataMax = outputTraceFloats[n];

					if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted avg Velocity " + outputTraceFloats[n]);
				}
			}

			// errorCorrection is the same for the velocity points below bottom markerSurface
			// botError stays the same
			topSlice = botSlice + 1;
			botSlice = nSlices - 1;
            boolean extrapolateError = false;
            for (n = topSlice; n <= botSlice; n++)
			{
				float avgVelocity = (float)errorCorrection*inputTraceFloats[n];
				if (!extrapolateError && (avgVelocity < 0.0f || avgVelocity > maxReasonableVelocity))
					extrapolateError = true;
				if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted avg Velocity " + avgVelocity);
				outputTraceFloats[n] = avgVelocity;
				if (dataMin > avgVelocity)
					 dataMin = avgVelocity;
				 if (dataMax < avgVelocity)
					 dataMax = avgVelocity;
			}
            if(extrapolateError)
                System.out.println("ERROR. Extrapolated Avg Velocity exceeded " + maxReasonableVelocity + " at row " + volumeRow + " col " + volumeCol);
            float newDepthMax = depthDatum + outputTraceFloats[nSlices-1] * (timeMax - timeDatum);

            setDepthMax(newDepthMax, volumeRow, volumeCol);
        /*
            if(newDepthMax < 0.0f || newDepthMax > 30000.0f)
            {
                System.out.println("ERROR. Max depth is NaN at row " + volumeRow + " col " + volumeCol);
            }
            else
                depthMax = Math.max(depthMax, newDepthMax);
        */
   //         StsMath.scale(outputTraceFloats, (float)(1.0/scaleMultiplier));
            rowFloatBuffer.put(outputTraceFloats);
            //adjustDataRange(outputTraceFloats);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException(
				"StsSeismicVelocityConstructor.processVelocityVolumeTrace() failed at  trace " + nTrace +
				" sample " + n + " surface " + s,
				e, StsException.WARNING);
                e.printStackTrace();
			return false;
		}
	}



	// jbw interval velocity seismic volume input, to create average velocity volume
final private boolean processVelocityIntervalVolumeTrace(int volumeRow, int volumeCol, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume outputVolume, double multiplier, float timeInc)
{
int n = -1, s = -1;
//System.out.println("experimental processVelocityIntervalAverageVolumeTrace");
boolean debug = false;
//	if (volumeRow < 10 && volumeCol < 10)
//			debug=true;
try
{
	int modelRow = volumeRow + volumeRowOffset;
	int modelCol = volumeCol + volumeColOffset;
	modelRow = StsMath.minMax(modelRow, 0, nModelRows - 1);
	modelCol = StsMath.minMax(modelCol, 0, nModelCols - 1);

	if (StsSeismicVelocityModel.debug == true && volumeRow == StsSeismicVelocityModel.debugVolumeRow && volumeCol == StsSeismicVelocityModel.debugVolumeCol)
	{
		debug = false;
		System.out.println("processVelocityVolumeTrace mainDebug. file " + rowFloatBuffer.pathname + " volume row " + volumeRow + " volume col " + volumeCol);
		if (rowFloatBuffer != null)
		{
			long position = rowFloatBuffer.getBufferPosition();
			long correctPosition = volumeRow * nInputCols * nSlices + volumeCol * nSlices;
			System.out.println("    start write position " + position + " correct position " + correctPosition);
		}
	}

//            double[] traceErrorCorrections = getTraceErrorCorrections(offsetRow, offsetCol);
	float[] traceTimeValues = getTraceTimeSurfaces(modelRow, modelCol);
	float[] traceDepthValues = getTraceDepthSurfaces(modelRow, modelCol);
//            double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);
	float[] inputTraceAvg = new float [inputTraceFloats.length];
	System.arraycopy(inputTraceFloats,0,inputTraceAvg,0,inputTraceFloats.length);
	StsMath.convertInstantToAvgVelocity(inputTraceAvg, 1.0d, timeInc);
	float botTime = timeDatum;
	double botError = 1.0f;
	float botDepth;
	double botAvgVolumeVelocity;
	double botAvgSurfaceVelocity;
	double errorCorrection = 1.0;

	int botSlice = -1;

	double topError;
	int topSlice;

	for (s = 0; s < nSurfaces; s++)
	{
		topError = botError;
		botTime = traceTimeValues[s];
		botDepth = traceDepthValues[s];
		float botSliceF = outputVolume.getBoundedSliceCoor(botTime);
		botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botSliceF, inputTraceAvg);
		botAvgSurfaceVelocity = (botDepth - depthDatum) / (botTime - timeDatum);
		botError = botAvgSurfaceVelocity / botAvgVolumeVelocity;
		if(debug) System.out.println("Surface: " + timeSurfaces[s].getName() + " surfaceAvgVelocity " + botAvgSurfaceVelocity +
									" botAvgVolumeVelocity " + botAvgVolumeVelocity + " error correction: " + botError);
		topSlice = botSlice + 1;
		botSlice = (int) botSliceF;
		if (debug)
				System.out.println("surface "+s+"topError "+topError+" bot "+botError);
		for (n = topSlice; n <= botSlice; n++)
		{
			float f = 0.0f;
			if(botSlice != topSlice) f = ( (float) (n - topSlice)) / (botSlice - topSlice);
			errorCorrection = topError + f * (botError - topError);
			outputTraceFloats[n] = (float)errorCorrection*inputTraceFloats[n];
			if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted INST Velocity " + outputTraceFloats[n]);
		}
	}

	// errorCorrection is the same for the velocity points below bottom markerSurface
	// botError stays the same
	topSlice = botSlice + 1;
	botSlice = nSlices - 1;
	boolean extrapolateError = false;
	for (n = topSlice; n <= botSlice; n++)
	{

		float avgVelocity = (float)errorCorrection*inputTraceFloats[n];
		if (!extrapolateError && (avgVelocity < 0.0f || avgVelocity > maxReasonableVelocity))
			extrapolateError = true;
		if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted avg Velocity " + avgVelocity);
		outputTraceFloats[n] = avgVelocity;
	}
	if(extrapolateError)
		System.out.println("ERROR. Extrapolated Avg Velocity exceeded " + maxReasonableVelocity + " at row " + volumeRow + " col " + volumeCol);
	//float newDepthMax = depthDatum + outputTraceFloats[nCroppedSlices-1] * (timeMax - timeDatum);
	//setDepthMax(newDepthMax, volumeRow, volumeCol);

	//StsMath.convertInstantToAvgVelocity(outputTraceFloats, 1.0d, timeInc);
	for (n = nSlices-1; n >0; n--)
	{

		if (dataMin > outputTraceFloats[n]) dataMin = (float)outputTraceFloats[n];
		if (dataMax < outputTraceFloats[n]) dataMax = (float)outputTraceFloats[n];
	}
	//float newDepthMax = depthDatum + outputTraceFloats[nCroppedSlices-1] * (timeMax - timeDatum);
	rowFloatBuffer.put(outputTraceFloats);
	//setDepthMax(newDepthMax, volumeRow, volumeCol);

	//adjustDataRange(outputTraceFloats);
	return true;
}
catch (Exception e)
{
	StsException.outputException(
		"StsSeismicVelocityConstructor.processVelocityVolumeTrace() failed at  trace " + nTrace +
		" sample " + n + " surface " + s,
		e, StsException.WARNING);
		e.printStackTrace();
	return false;
}
}




	final private boolean processVelocityIntervalVolumeTrace2(int volumeRow, int volumeCol, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume outputVolume, double multiplier, float timeInc)
	{
	int n = -1, s = -1;
	//System.out.println("experimental processVelocityIntervalVolumeTrace  -- interval interval"+volumeRow+" "+volumeCol);
	boolean debug = false;
//	if (volumeRow < 10 && volumeCol < 10)
//			debug=true;
	try
	{
		int modelRow = volumeRow + volumeRowOffset;
		int modelCol = volumeCol + volumeColOffset;
		modelRow = StsMath.minMax(modelRow, 0, nModelRows - 1);
		modelCol = StsMath.minMax(modelCol, 0, nModelCols - 1);

		debug=false;

		if (StsSeismicVelocityModel.debug == true && volumeRow == StsSeismicVelocityModel.debugVolumeRow && volumeCol == StsSeismicVelocityModel.debugVolumeCol)
		{

			System.out.println("processVelocityVolumeTrace mainDebug. file " + rowFloatBuffer.pathname + " volume row " + volumeRow + " volume col " + volumeCol);
			if (rowFloatBuffer != null)
			{
				long position = rowFloatBuffer.getBufferPosition();
				long correctPosition = volumeRow * nInputCols * nSlices + volumeCol * nSlices;
				System.out.println("    start write position " + position + " correct position " + correctPosition);
			}
		}

//            double[] traceErrorCorrections = getTraceErrorCorrections(offsetRow, offsetCol);
		float[] traceTimeValues = getTraceTimeSurfaces(modelRow, modelCol);
		float[] traceDepthValues = getTraceDepthSurfaces(modelRow, modelCol);
//            double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);
		float[] inputTraceAvg = new float [inputTraceFloats.length];
		System.arraycopy(inputTraceFloats,0,inputTraceAvg,0,inputTraceFloats.length);
		StsMath.convertInstantToAvgVelocity(inputTraceAvg, 1.0d, timeInc);



		float botTime = timeDatum;
		double botError = 1.0f;
		float botDepth;
		float topTime = timeDatum;
		float topDepth = depthDatum;
		double botAvgVolumeVelocity, topAvgVolumeVelocity;
		double botAvgSurfaceVelocity, topAvgSurfaceVelocity;
		double errorCorrection = 1.0;

		int botSlice = -1;

		double topError = 1.0;
		int topSlice = 0;

		for (s = 0; s < nSurfaces; s++)
		{

			botTime = traceTimeValues[s];
			botDepth = traceDepthValues[s];
			float botSliceF = outputVolume.getBoundedSliceCoor(botTime);

			botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botSliceF, inputTraceAvg);
			botAvgSurfaceVelocity = (botDepth - depthDatum) / (botTime - timeDatum);


			botError = botAvgSurfaceVelocity / botAvgVolumeVelocity;
			if(debug) System.out.println("Surface: " + timeSurfaces[s].getName() + " surfaceAvgVelocity " + botAvgSurfaceVelocity +
										" botAvgVolumeVelocity " + botAvgVolumeVelocity + " error correction: " + botError);

			botSlice = (int) botSliceF;
			topSlice = botSlice + 1;
			for (n = topSlice; n <= botSlice; n++)
			{
				float f = 0.0f;
				if(botSlice != topSlice) f = ( (float) (n - topSlice)) / (botSlice - topSlice);
				errorCorrection = topError + f * (botError - topError);
				outputTraceFloats[n] = (float)errorCorrection*inputTraceFloats[n];
				if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted INST Velocity " + outputTraceFloats[n]);

			}
			topTime = botTime;
			topDepth = botDepth;
			topError = botError;
		}

		// errorCorrection is the same for the velocity points below bottom markerSurface
		// botError stays the same

		topSlice = botSlice + 1;
		botSlice = nSlices - 1;
		boolean extrapolateError = false;
		for (n = topSlice; n <= botSlice; n++)
		{
			outputTraceFloats[n] = (float)botError*inputTraceFloats[n];
		}

	    for (n = nSlices-1; n >0; n--)
	    {
			if (outputTraceFloats[n] > 0)
			{
				if(dataMin > outputTraceFloats[n])dataMin = (float)outputTraceFloats[n];
				if(dataMax < outputTraceFloats[n])dataMax = (float)outputTraceFloats[n];
			}
	    }


		rowFloatBuffer.put(outputTraceFloats);
		//adjustDataRange(outputTraceFloats);
		return true;
	}
	catch (Exception e)
	{
		StsException.outputException(
			"StsSeismicVelocityConstructor.processVelocityVolumeTrace() failed at  trace " + nTrace +
			" sample " + n + " surface " + s,
			e, StsException.WARNING);
			e.printStackTrace();
		return false;
	}
}

// jbw interval velocity seismic volume input, to create average velocity volume
final private boolean processVelocityIntervalAverageVolumeTrace(int volumeRow, int volumeCol, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume outputVolume, double multiplier, float timeInc)
{
int n = -1, s = -1;
//System.out.println("experimental processVelocityIntervalAverageVolumeTrace");
boolean debug = false;
	if (volumeRow < 10 && volumeCol < 10)
			debug=true;
try
{
	int modelRow = volumeRow + volumeRowOffset;
	int modelCol = volumeCol + volumeColOffset;
	modelRow = StsMath.minMax(modelRow, 0, nModelRows - 1);
	modelCol = StsMath.minMax(modelCol, 0, nModelCols - 1);

	if (StsSeismicVelocityModel.debug == true && volumeRow == StsSeismicVelocityModel.debugVolumeRow && volumeCol == StsSeismicVelocityModel.debugVolumeCol)
	{
		debug = false;
		System.out.println("processVelocityVolumeTrace mainDebug. file " + rowFloatBuffer.pathname + " volume row " + volumeRow + " volume col " + volumeCol);
		if (rowFloatBuffer != null)
		{
			long position = rowFloatBuffer.getBufferPosition();
			long correctPosition = volumeRow * nInputCols * nSlices + volumeCol * nSlices;
			System.out.println("    start write position " + position + " correct position " + correctPosition);
		}
	}

//            double[] traceErrorCorrections = getTraceErrorCorrections(offsetRow, offsetCol);
	float[] traceTimeValues = getTraceTimeSurfaces(modelRow, modelCol);
	float[] traceDepthValues = getTraceDepthSurfaces(modelRow, modelCol);
//            double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);
	float[] inputTraceAvg = new float [inputTraceFloats.length];
	System.arraycopy(inputTraceFloats,0,inputTraceAvg,0,inputTraceFloats.length);
	StsMath.convertInstantToAvgVelocity(inputTraceAvg, 1.0d, timeInc);
	float botTime = timeDatum;
	double botError = 1.0f;
	float botDepth;
	double botAvgVolumeVelocity;
	double botAvgSurfaceVelocity;
	double errorCorrection = 1.0;

	int botSlice = -1;

	double topError;
	int topSlice;

	for (s = 0; s < nSurfaces; s++)
	{
		topError = botError;
		botTime = traceTimeValues[s];
		botDepth = traceDepthValues[s];
		float botSliceF = outputVolume.getBoundedSliceCoor(botTime);
		botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botSliceF, inputTraceAvg);
		botAvgSurfaceVelocity = (botDepth - depthDatum) / (botTime - timeDatum);
		botError = botAvgSurfaceVelocity / botAvgVolumeVelocity;
		if(debug) System.out.println("Surface: " + timeSurfaces[s].getName() + " surfaceAvgVelocity " + botAvgSurfaceVelocity +
									" botAvgVolumeVelocity " + botAvgVolumeVelocity + " error correction: " + botError);
		topSlice = botSlice + 1;
		botSlice = (int) botSliceF;
		if (debug)
				System.out.println("surface "+s+"topError "+topError+" bot "+botError);
		for (n = topSlice; n <= botSlice; n++)
		{
			float f = 0.0f;
			if(botSlice != topSlice) f = ( (float) (n - topSlice)) / (botSlice - topSlice);
			errorCorrection = topError + f * (botError - topError);
			outputTraceFloats[n] = (float)errorCorrection*inputTraceFloats[n];
			if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted INST Velocity " + outputTraceFloats[n]);
		}
	}

	// errorCorrection is the same for the velocity points below bottom markerSurface
	// botError stays the same
	topSlice = botSlice + 1;
	botSlice = nSlices - 1;
	boolean extrapolateError = false;
	for (n = topSlice; n <= botSlice; n++)
	{

		float avgVelocity = (float)errorCorrection*inputTraceFloats[n];
		if (!extrapolateError && (avgVelocity < 0.0f || avgVelocity > maxReasonableVelocity))
			extrapolateError = true;
		if (debug && (n == topSlice || n == botSlice)) System.out.println("    slice " + n + " adjusted avg Velocity " + avgVelocity);
		outputTraceFloats[n] = avgVelocity;
	}
	if(extrapolateError)
		System.out.println("ERROR. Extrapolated Avg Velocity exceeded " + maxReasonableVelocity + " at row " + volumeRow + " col " + volumeCol);
	//float newDepthMax = depthDatum + outputTraceFloats[nCroppedSlices-1] * (timeMax - timeDatum);
	//setDepthMax(newDepthMax, volumeRow, volumeCol);

	StsMath.convertInstantToAvgVelocity(outputTraceFloats, 1.0d, timeInc);
	for (n = nSlices-1; n >0; n--)
	{

		if (dataMin > outputTraceFloats[n]) dataMin = (float)outputTraceFloats[n];
		if (dataMax < outputTraceFloats[n]) dataMax = (float)outputTraceFloats[n];
	}
	float newDepthMax = depthDatum + outputTraceFloats[nSlices-1] * (timeMax - timeDatum);
	rowFloatBuffer.put(outputTraceFloats);
	setDepthMax(newDepthMax, volumeRow, volumeCol);
	model.getProject().checkChangeZRange(depthDatum, newDepthMax, StsProject.TD_DEPTH);
	//adjustDataRange(outputTraceFloats);
	return true;
}
catch (Exception e)
{
	StsException.outputException(
		"StsSeismicVelocityConstructor.processVelocityVolumeTrace() failed at  trace " + nTrace +
		" sample " + n + " surface " + s,
		e, StsException.WARNING);
		e.printStackTrace();
	return false;
}
}


    private void setDepthMax(float newDepthMax, int row, int col)
    {
        if(newDepthMax > depthMax)
        {
            depthMax = newDepthMax;
//            StsException.systemDebug(this, "setDepthMax", "new depth max: " + depthMax + " row: " + row + " col: " + col);
        }
    }

    /*
		final private boolean processVelocityVolumeTrace(int row, int col,
			int nTrace, StsMappedByteBuffer inlineBuffer,
			StsMappedFloatBuffer inlineFloatBuffer, StsSeismicVolume outputVolume)
		{
			int n = 0;

			try
			{
				int offsetRow = row + modelRowOffset;
				int offsetCol = col + modelColOffset;
				offsetRow = StsMath.minMax(offsetRow, 0, nModelRows - 1);
				offsetCol = StsMath.minMax(offsetCol, 0, nModelCols - 1);
//            double[] traceErrorCorrections = getTraceErrorCorrections(offsetRow, offsetCol);
				float[] traceTimeValues = getTraceTimeSurfaces(offsetRow, offsetCol);
				float[] traceDepthValues = getTraceDepthSurfaces(offsetRow, offsetCol);
//            double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);

				float topTime = timeDatum;
				float botTime = traceTimeValues[0];
				float topDepth = depthDatum;
				float botDepth = traceDepthValues[0];
				double topError = 1.0f;
				float botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botTime);
				double botAvgSurfaceVelocity = (botDepth - depthDatum)/(botTime - timeDatum);
				double botError = botAvgSurfaceVelocity/botAvgVolumeVelocity;
				float t = timeMin;
				int nBotSurface = 0;
				double errorCorrection = 1.0f;
				byte b = -1;
				for (n = 0; n < nInputSlices; n++, t += timeInc)
				{
					if (t > botTime && nBotSurface < nSurfaces - 1)
					{
						topTime = botTime;
						topDepth = botDepth;
						topError = botError;
						nBotSurface++;
						botTime = traceTimeValues[nBotSurface];
						botDepth = traceDepthValues[nBotSurface];
						botAvgVolumeVelocity = getTraceAvgVelocityAtTime(botTime);
						botAvgSurfaceVelocity = (botDepth - depthDatum)/(botTime - timeDatum);
						botError = botAvgSurfaceVelocity/botAvgVolumeVelocity;
					}

					float avgVelocity = inputTraceFloats[n];
					if (t <= botTime)
					{
						float f = (t - topTime) / (botTime - topTime);
						errorCorrection = topError + f * (botError - topError);
					}
					avgVelocity *= errorCorrection;
					int scaledValue = (int) (avgVelocity * scaleFactor + scaleOffset);
					scaledValue = StsMath.minMax(scaledValue, 0, 254);
					outputVolume.accumulateHistogram(scaledValue);
					b = StsMath.UnsignedIntToSignedByte(scaledValue);
					outputTraceBytes[n] = b;
					outputTraceFloats[n] = (float)avgVelocity;
				}
				inlineBuffer.put(outputTraceBytes);
				if (inlineFloatBuffer != null) inlineFloatBuffer.put(outputTraceFloats);

				return true;
			}
			catch (Exception e)
			{
				StsException.outputException(
					"StsSeismicVelocityConstructor.processVelocityVolumeTrace() failed at  trace " + nTrace +
					" sample " + n,
					e, StsException.WARNING);
				return false;
			}
		}
	 */

	private float getTraceAvgVelocityAtTime(float sliceF)
	{
		int slice = (int) sliceF;
		slice = StsMath.minMax(slice, 0, nSlices - 2);
		float dSlice = sliceF - slice;
		return inputTraceFloats[slice] * (1.0f - dSlice) + inputTraceFloats[slice + 1] * dSlice;
	}
	private float getTraceAvgVelocityAtTime(float sliceF, float[]trace)
	{
		int slice = (int) sliceF;
		slice = StsMath.minMax(slice, 0, nSlices - 2);
		float dSlice = sliceF - slice;
		return trace[slice] * (1.0f - dSlice) + trace[slice + 1] * dSlice;
	}

	// build an average velocity trace from interval velocity information
	/*
		final private boolean processIntervalVelocityVolumeTrace(int row, int col,
			int nTrace, StsMappedByteBuffer inlineBuffer,
			StsMappedFloatBuffer inlineFloatBuffer, StsSeismicVolume volume)
		{
			int n = 0;

			try
			{
				int offsetRow = row + modelRowOffset;
				int offsetCol = col + modelColOffset;
				offsetRow = StsMath.minMax(offsetRow, 0, nModelRows - 1);
				offsetCol = StsMath.minMax(offsetCol, 0, nModelCols - 1);

				float[] traceTimeValues = getTraceTimeSurfaces(offsetRow, offsetCol);
				double[] traceVelocities = getTraceVelocities(offsetRow, offsetCol);

				int nBotSurface = 0;
				float botTime = traceTimeValues[0];
				double intervalVelocity = traceVelocities[nBotSurface];
				int scaledAvgVelocity = (int) (intervalVelocity * scaleFactor + scaleOffset);
				scaledAvgVelocity = StsMath.minMax(scaledAvgVelocity, 0, 254);
				byte scaledByte = StsMath.UnsignedIntToSignedByte(scaledAvgVelocity);

				float t = timeMin;
				byte b = -1;
				float dt = 0; // change in t from t = timeMin datum
				float dz = 0; // change in z at t = timeMin datum
				double v1 = traceVelocities[0];
				scaledAvgVelocity = (int) (v1 * scaleFactor + scaleOffset);
				scaledAvgVelocity = StsMath.minMax(scaledAvgVelocity, 0, 254);
				scaledByte = StsMath.UnsignedIntToSignedByte(scaledAvgVelocity);
				volume.accumulateHistogram(scaledAvgVelocity);
				outputTraceBytes[0] = scaledByte;
				outputTraceFloats[0] = (float)v1;
				t += timeInc;
				for (n = 1; n < nInputSlices; n++, t += timeInc)
				{
					if (t > botTime)
					{
						if (nBotSurface < nSurfaces - 1)
						{
							nBotSurface++;
							botTime = traceTimeValues[nBotSurface];
							intervalVelocity = traceVelocities[nBotSurface];
						}
					}
					float v0 = (float)v1;
					v1 = intervalVelocity;
					dz += (v0 + v1)*timeInc/2;
					dt += timeInc;
					float avgVelocity = dz/dt;
					scaledAvgVelocity = (int) (avgVelocity * scaleFactor + scaleOffset);
					scaledAvgVelocity = StsMath.minMax(scaledAvgVelocity, 0, 254);
					volume.accumulateHistogram(scaledAvgVelocity);
					scaledByte = StsMath.UnsignedIntToSignedByte(scaledAvgVelocity);
					outputTraceBytes[n] = scaledByte;
					outputTraceFloats[n] = avgVelocity;
				}
				inlineBuffer.put(outputTraceBytes);
				if (inlineFloatBuffer != null) inlineFloatBuffer.put(outputTraceFloats);
				return true;
			}
			catch (Exception e)
			{
				StsException.outputException(
					"StsSeismicVelocityConstructor.processTrace() failed at  trace " + nTrace +
					" sample " + n,
					e, StsException.WARNING);
				return false;
			}
		}
	 */
	final private boolean processIntervalVelocityVolumeTrace(int row, int col, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume volume)
	{
		int n = 0;

		try
		{
			int offsetRow = row + volumeRowOffset;
			int offsetCol = col + volumeColOffset;
			offsetRow = StsMath.minMax(offsetRow, 0, nModelRows - 1);
			offsetCol = StsMath.minMax(offsetCol, 0, nModelCols - 1);
			/** surface times at this row-col */
			float[] traceSurfaceTimeValues = getTraceTimeSurfaces(offsetRow, offsetCol);
			/** two-way interval velocities at this row-col between surface pairs: datum & first horizon down to last two horizons */
			double[] traceIntervalVelocities = getTraceIntervalVelocities(offsetRow, offsetCol);
			/** surface depths at this row-col */
			float[] traceSurfaceDepthValues = getTraceDepthSurfaces(offsetRow, offsetCol);
			double topDepth = depthDatum;
			float topTime = timeDatum;
			float botTime = traceSurfaceTimeValues[0];
			double botDepth = traceSurfaceDepthValues[0];
			double intervalVelocity = traceIntervalVelocities[0];
			int nBotSurface = 0;

			double avgVelocity;

			float t = timeMin;
			for (n = 0; n < nSlices; n++, t += timeInc)
			{
				if(t == timeDatum)
					avgVelocity = intervalVelocity;
				else
				{
					if(t > botTime)
					{
						if(nBotSurface < nSurfaces-1)
						{
							while(true)
							{
								topTime = botTime;
								topDepth = botDepth;
								nBotSurface++;

								botTime = traceSurfaceTimeValues[nBotSurface];
								if (t <= botTime)
								{
									botDepth = traceSurfaceDepthValues[nBotSurface];
									intervalVelocity = traceIntervalVelocities[nBotSurface];
									break;
								}
							}
						}
						else if(nBotSurface == nSurfaces-1)
						{
							if(oneWayConstantBottomIntervalVelocity != 0.0f)
								intervalVelocity = oneWayConstantBottomIntervalVelocity/2;
							topTime = botTime;
							topDepth = botDepth;
							botTime = timeMax;
							nBotSurface++;
						}
					}
					avgVelocity = (topDepth + intervalVelocity * (t - topTime) - depthDatum) / (t - timeDatum);
				}
				if (dataMin > avgVelocity) dataMin = (float)avgVelocity;
				if (dataMax < avgVelocity) dataMax = (float)avgVelocity;
				outputTraceFloats[n] = (float) avgVelocity;
			}
			depthMax = Math.max(depthMax, depthDatum + outputTraceFloats[nSlices-1] * (timeMax - timeDatum));

			rowFloatBuffer.put(outputTraceFloats);
            if(StsSeismicVelocityModel.debug && row == StsSeismicVelocityModel.debugVolumeRow && col == StsSeismicVelocityModel.debugVolumeCol)
                System.out.println("processIntervalVelocityTrace(" + row + ", " + col + ") outputTraceFloats[0] = " + outputTraceFloats[0]);
            adjustDataRange(outputTraceFloats);
            return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "processTrace", "failed at  trace " + nTrace + " sample " + n, e);
			return false;
		}
	}
	final private boolean processIntervalVelocityIntervalVolumeTrace(int row, int col, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume volume)
	{
		int n = 0;

		try
		{
			int offsetRow = row + volumeRowOffset;
			int offsetCol = col + volumeColOffset;
			offsetRow = StsMath.minMax(offsetRow, 0, nModelRows - 1);
			offsetCol = StsMath.minMax(offsetCol, 0, nModelCols - 1);
			/** surface times at this row-col */
			float[] traceSurfaceTimeValues = getTraceTimeSurfaces(offsetRow, offsetCol);
			/** two-way interval velocities at this row-col between surface pairs: datum & first horizon down to last two horizons */
			double[] traceIntervalVelocities = getTraceIntervalVelocities(offsetRow, offsetCol);
			/** surface depths at this row-col */
			float[] traceSurfaceDepthValues = getTraceDepthSurfaces(offsetRow, offsetCol);
			float botTime = traceSurfaceTimeValues[0];
			double intervalVelocity = traceIntervalVelocities[0];
			int nBotSurface = 0;

			float t = timeMin;
			for (n = 0; n < nSlices; n++, t += timeInc)
			{
				if(t > botTime)
				{
					if(nBotSurface < nSurfaces-1)
					{
						while(true)
						{
							nBotSurface++;

							botTime = traceSurfaceTimeValues[nBotSurface];
							if (t <= botTime)
							{
								intervalVelocity = traceIntervalVelocities[nBotSurface];
								break;
							}
						}
					}
					else if(nBotSurface == nSurfaces-1)
					{
						if(oneWayConstantBottomIntervalVelocity != 0.0f)
							intervalVelocity = oneWayConstantBottomIntervalVelocity/2;
						botTime = timeMax;
						nBotSurface++;
					}
				}
				if (dataMin > intervalVelocity) dataMin = (float)intervalVelocity;
				if (dataMax < intervalVelocity) dataMax = (float)intervalVelocity;

				outputTraceFloats[n] = (float) intervalVelocity;
			}

			rowFloatBuffer.put(outputTraceFloats);
			if(StsSeismicVelocityModel.debug && row == StsSeismicVelocityModel.debugVolumeRow && col == StsSeismicVelocityModel.debugVolumeCol)
				System.out.println("processIntervalVelocityTrace(" + row + ", " + col + ") outputTraceFloats[0] = " + outputTraceFloats[0]);
			adjustDataRange(outputTraceFloats);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "processIntervalVelocityIntervalVolumeTrace", "failed at  trace " + nTrace + " sample " + n, e);
			return false;
		}
	}

	final private boolean processIntervalVelocityIntervalVolumeTraceX(int row, int col, int nTrace, StsMappedFloatBuffer rowFloatBuffer, StsSeismicVolume volume)
	{
		int n = 0;

		try
		{
			int offsetRow = row + volumeRowOffset;
			int offsetCol = col + volumeColOffset;
			offsetRow = StsMath.minMax(offsetRow, 0, nModelRows - 1);
			offsetCol = StsMath.minMax(offsetCol, 0, nModelCols - 1);

			float[] traceTimeValues = getTraceTimeSurfaces(offsetRow, offsetCol);
			double[] traceIntervalVelocities = getTraceIntervalVelocities(offsetRow, offsetCol);
			float[] traceDepthValues = getTraceDepthSurfaces(offsetRow, offsetCol);
			double topDepth = depthDatum;
			int nBotSurface = 0;
			float topTime = timeDatum;
			float botTime = traceTimeValues[0];
			double botDepth = traceDepthValues[0];
			double intervalVelocity = traceIntervalVelocities[0];

			float t = timeMin;

			for (n = 0; n < nSlices; n++, t += timeInc)
			{

				if (t > botTime)
				{
					if (nBotSurface < nSurfaces -1)
					{
						topTime = botTime;
						topDepth = botDepth;
						nBotSurface++;
						botTime = traceTimeValues[nBotSurface];
                        if (t > botTime)     // jbw zero thickness
                            continue;
						botDepth = traceDepthValues[nBotSurface];
						intervalVelocity = traceIntervalVelocities[nBotSurface];
						if (intervalVelocity > 0)
							if (dataMin > intervalVelocity) dataMin = (float)intervalVelocity;
						if (dataMax < intervalVelocity) dataMax = (float)intervalVelocity;

					}

				}
				if (dataMin > intervalVelocity) dataMin = (float)intervalVelocity;
				if (dataMax < intervalVelocity) dataMax = (float)intervalVelocity;

				outputTraceFloats[n] = (float) intervalVelocity;
			}
			// nooo  depthMax = Math.max(depthMax, (depthDatum + outputTraceFloats[nCroppedSlices-1] * (timeMax - timeDatum)));

			rowFloatBuffer.put(outputTraceFloats);
			if(StsSeismicVelocityModel.debug && row == StsSeismicVelocityModel.debugVolumeRow && col == StsSeismicVelocityModel.debugVolumeCol)
				System.out.println("processIntervalVelocityTrace(" + row + ", " + col + ") outputTraceFloats[0] = " + outputTraceFloats[0]);
			adjustDataRange(outputTraceFloats);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "processTrace", "failed at  trace " + nTrace + " sample " + n, e);
			return false;
		}
	}


    public void finalizeVolumeOutput()
    {
        super.finalizeVolumeOutput();
        model.getProject().checkSetDepthMax(depthMax);
    }

    private void adjustDataRange(float[] outputTraceFloats)
    {
        float newDataMin = outputTraceFloats[0];
        if(newDataMin != StsParameters.nullValue && newDataMin < dataMin)
            dataMin = newDataMin;
        float newDataMax = outputTraceFloats[nSlices-1];
        if(newDataMax != StsParameters.nullValue && newDataMax > dataMax)
            dataMax = newDataMax;
    }

    /*
		private double[] getTraceErrorCorrections(int row, int col)
		{
			double[] errorCorrections = new double[nSurfaces];
			for (int n = 0; n < nSurfaces; n++)
			{
				errorCorrections[n] = this.surfaceCorrectionFactors[n][row][col];
			}
			return errorCorrections;
		}
	 */
	private float[] getTraceTimeSurfaces(int row, int col)
	{
		float[] timeSurfaces = new float[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
			timeSurfaces[n] = this.surfaceTimes[n][row][col];
		}
		return timeSurfaces;
	}

	private float[] getTraceDepthSurfaces(int row, int col)
	{
		float[] depthSurfaces = new float[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
			depthSurfaces[n] = this.surfaceDepths[n][row][col];
		}
		return depthSurfaces;
	}

	private double[] getTraceIntervalVelocities(int row, int col)
	{
		double[] velocities = new double[nSurfaces];
		for (int n = 0; n < nSurfaces; n++)
		{
			velocities[n] = this.surfaceVelocities[n][row][col];
		}
		return velocities;
	}
}
