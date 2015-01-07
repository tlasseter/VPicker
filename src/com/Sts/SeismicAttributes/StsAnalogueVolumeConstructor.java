package com.Sts.SeismicAttributes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Actions.Wizards.AnalogueCube.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsAnalogueVolumeConstructor extends StsSeismicVolumeConstructor
{
	boolean isDataFloat = false;
	boolean isComplex = false;
	boolean useByteCubes = true;
	byte spiralMethod = StsAnalogueCubeWizard.RUNNING_AVERAGE;
	String name = "Analogue Cube";
	StsSeismicVolume sourceVolume;
	StsSeismicVolume sourceVolumeH;
	StsBoxSubVolume sourceSubVolume;
	StsSeismicVolume targetVolume;
	StsSeismicVolume targetVolumeH;
	StsSubVolume targetSubVolume;

	int sourceRowMin, sourceRowMax;
	int sourceColMin, sourceColMax;
	int sourceSliceMin, sourceSliceMax;
	int sourceRowCenter, sourceColCenter, sourceSliceCenter;
	int nSourceRows, nSourceCols, nSourceSlices;

//    int sourceColStart, sourceColEnd;
//    int sourceSliceStart, sourceSliceEnd;

	int nTargetRows, nTargetCols, nTargetSlices;

//    int targetSliceStart, targetSliceEnd;

	int targetRowMin, targetRowMax, targetColMin, targetColMax, targetSliceMin, targetSliceMax;

	int maxNSpirals = 0;

	StsGridBoundingBox sourceBoundingBox;
	StsGridBoundingBox targetBoundingBox;

	double[] targetValues, targetValuesH;
	float correlThreshold = -1.0f;

	double[][][] sourceCube, sourceCubeH;

	byte[] transparentRowPlane = null;
	float[] nullFloatRowPlane = null;
	byte[] transparentTrace = null;
	float[] nullTrace = null;

	RunParameters parameters = null; // used only by fullRun;  dryRun must construct separate parameters instances for each cursor being run

	int progressBarMax = 0;
	int progressBarValue = 0;

	static int minNumProcessed = 9;
	static int nToSkip = 20;
	static int minNToSkip = 2;

	StsTimer timer = null;
	boolean runTimer = false;
	boolean blockDone = false;

	static boolean fastTest = false;
	static boolean debug = false;

	final static int XDIR = StsCursor3d.XDIR;
	final static int YDIR = StsCursor3d.YDIR;
	final static int ZDIR = StsCursor3d.ZDIR;

	final static public String volumeFilePrefix = "seis.vol.";

	private StsAnalogueVolumeConstructor(StsModel model, StsSeismicVolume sourceVolume, StsSeismicVolume targetVolume,
							  StsBoxSubVolume sourceSubVolume, StsSubVolume targetSubVolume, boolean isDataFloat, byte type,
							  byte spiralMethod, StsProgressPanel panel, boolean createFiles, String stemname)
	{

		if (runTimer) timer = new StsTimer();
		this.panel = panel;
		inputVolumes = new StsSeismicVolume[] { targetVolume };
		volumeName = ANALOGUE;
		this.model = model;
		this.sourceVolume = sourceVolume;
		this.targetVolume = targetVolume;
		this.sourceSubVolume = sourceSubVolume;
		this.targetSubVolume = targetSubVolume;
		this.isDataFloat = isDataFloat;
		this.spiralMethod = spiralMethod;

		isComplex = type == StsAnalogueCubeWizard.COMPLEX;
		outputVolume = StsSeismicVolume.initializeAttributeVolume(model, targetVolume, -1.0f, 1.0f, true, createFiles, stemname, volumeName, "rw");

		StsSpectrumClass spectrumClass = model.getSpectrumClass();
		outputVolume.getColorscale().setSpectrum(spectrumClass.getSpectrum("Default Edge Spectrum"));

		sourceBoundingBox = sourceSubVolume.getGridBoundingBox();
		if (targetSubVolume != null)
		{
			targetBoundingBox = targetSubVolume.getGridBoundingBox();
		}
		else
		{
			targetBoundingBox = targetVolume.getGridBoundingBox(false);

		}
		sourceRowMin = sourceBoundingBox.rowMin;
		sourceRowMax = sourceBoundingBox.rowMax;
		sourceColMin = sourceBoundingBox.colMin;
		sourceColMax = sourceBoundingBox.colMax;
		sourceSliceMin = sourceBoundingBox.sliceMin;
		sourceSliceMax = sourceBoundingBox.sliceMax;

		nSourceRows = sourceBoundingBox.getNRows();
		nSourceCols = sourceBoundingBox.getNCols();
		nSourceSlices = sourceBoundingBox.getNSlices();

		sourceRowCenter = sourceSubVolume.getRowCenter() - sourceRowMin;
		sourceColCenter = sourceSubVolume.getColCenter() - sourceColMin;
		sourceSliceCenter = sourceSubVolume.getSliceCenter() - sourceSliceMin;

		if (debug)
		{
			int nSkipTraces = Math.max(nSourceRows * nSourceCols / nToSkip, minNToSkip);
			System.out.println("After min threshold, will compute only every " + nSkipTraces + " trace.");
		}

		if (fastTest)
			maxNSpirals = 0;
		else
			maxNSpirals = StsMath.max4(sourceRowCenter, sourceColCenter, nSourceRows - sourceRowCenter - 1,
									   nSourceCols - sourceColCenter - 1);

		sourceCube = new double[nSourceRows][nSourceCols][nSourceSlices];

	    targetVolume.clearCache();
		if(isDataFloat) targetVolume.setupReadRowFloatBlocks();

		if(sourceVolume != targetVolume)
		{
			sourceVolume.clearCache();
			if(isDataFloat) sourceVolume.setupReadRowFloatBlocks();
		}
		if (type == StsAnalogueCubeWizard.COMPLEX)
		{
			sourceCubeH = new double[nSourceRows][nSourceCols][nSourceSlices];
			sourceVolumeH = getHilbertVolume(sourceVolume);
			sourceVolumeH.clearCache();
			if(isDataFloat) sourceVolumeH.setupReadRowFloatBlocks();

			if(sourceVolume != targetVolume)
			{
				targetVolumeH = getHilbertVolume(targetVolume);
				if(isDataFloat) targetVolumeH.setupReadRowFloatBlocks();
			}
			else
				targetVolumeH = sourceVolumeH;

			model.setCurrentObjectDisplayAndToolbar(outputVolume);
		}
		nTargetRows = targetBoundingBox.getNRows();
		nTargetCols = targetBoundingBox.getNCols();
		nTargetSlices = targetBoundingBox.getNSlices();
		targetRowMin = targetBoundingBox.rowMin;
		targetRowMax = targetBoundingBox.rowMax;
		targetColMin = targetBoundingBox.colMin;
		targetColMax = targetBoundingBox.colMax;
		targetSliceMin = targetBoundingBox.sliceMin;
		targetSliceMax = targetBoundingBox.sliceMax;
		targetValues = new double[nTargetSlices];
		if (isComplex) targetValuesH = new double[nTargetSlices];
//		fillSourceCube(sourceVolume, sourceCube);
//		if (isComplex) fillSourceCube(sourceVolumeH, sourceCubeH);
	}

	public static StsAnalogueVolumeConstructor constructor(StsModel model, StsSeismicVolume sourceVolume,
												StsSeismicVolume targetVolume, StsBoxSubVolume sourceSubVolume,
												StsSubVolume targetSubVolume, boolean isDataFloat, byte type,
												byte spiralMethod, StsProgressPanel panel, boolean createFiles, String stemname)
	{
		try
		{
			return new StsAnalogueVolumeConstructor(model, sourceVolume, targetVolume, sourceSubVolume, targetSubVolume,
										 isDataFloat, type, spiralMethod, panel, createFiles, stemname);
		}
		catch (Exception e)
		{
			StsException.outputException("StsAnalogueVolumeConstructor.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public boolean checkCreateFiles()
	{
		outputVolume.close();
		outputVolume.createVolumeFilenames(volumeFilePrefix + name);
		outputVolume.writeHeaderFile();
		try
		{
			return outputVolume.initialize(model, "rw");
		}
		catch(Exception e)
		{
			outputVolume.delete();
			StsMessage.printMessage("Failed to classInitialize volume " + outputVolume.getName());
			return false;
		}
	}

	public void fillSourceCube()
	{
		fillSourceCube(sourceVolume, sourceCube);
		if(isComplex) fillSourceCube(sourceVolumeH, sourceCubeH);
	}

	private StsSeismicVolume getHilbertVolume(StsSeismicVolume inputVolume)
	{
		// check if volume already loaded

		String hilbertStemname = new String(inputVolume.stemname + "." + StsSeismicVolumeConstructor.HILBERT);
		StsSeismicVolume hilbertVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, hilbertStemname);
		if(hilbertVolume == null)
		{
			if(inputVolume.derivedFileOK(hilbertStemname))
            {
                hilbertVolume = StsSeismicVolume.checkLoadFromStemname(model, inputVolume.stsDirectory, hilbertStemname, true);
		    }
        }
        if(hilbertVolume != null && !inputVolume.sameAs(hilbertVolume))
		{
			hilbertVolume.delete();
			hilbertVolume = null;
		}

		// if volume files don't exist or are out of date, create and load

		if (hilbertVolume == null)
		{
			panel.appendLine("Creating Hilbert Volume for volume " + inputVolume.getName());
			StsHilbertTransformConstructor hilbert = StsHilbertTransformConstructor.constructor(model, inputVolume, true, panel);

			// since constructing the hilbert makes it the currentObject, we need to reset it to
			// the initial input object
			model.setCurrentObject(inputVolume);
			if (hilbert == null)
			{
				return null;
			}
			hilbertVolume = hilbert.getVolume();
			if(!hilbertVolume.initialize(model))
			{
				hilbertVolume.delete();
				return null;
			}
			panel.appendLine("Completed Hilbert Volume construction.");
		}
		return hilbertVolume;
	}

	public void initializeRun()
	{
		this.panel = panel;
		panel.initialize(nTargetRows);
	}

	public void dryRun(byte[] planeData, Integer dirObject, Integer nPlaneObject)
	{
		int row, col, slice;
		int n;
		double result;

		RunParameters parameters = new RunParameters();

		if (runTimer)
		{
			timer.start();

		}
		int dirNo = dirObject.intValue();
		int nPlane = nPlaneObject.intValue();
		int progressBarAdd = 0;
		if(outputVolume == null) return;
        byte zDomain = sourceVolume.getZDomain();
		switch(dirNo)
		{
			case XDIR:
				progressBarAdd = nInputRows * nInputSlices;
				progressBarMax += progressBarAdd;
				panel.setMaximum(progressBarMax);
				n = 0;
				col = nPlane;

				if (col < targetColMin || col > targetColMax)
				{
					for (row = 0; row < nInputRows; row++)
					{
						if (row % 10 == 0)
						{
							progressBarValue += nInputSlices * Math.min(10, (nInputRows - row));
							panel.setValue(progressBarValue);
						}
						for (slice = 0; slice < nInputSlices; slice++)
						{
							planeData[n++] = (byte)-1;
						}
					}
				}
				else
				{
					if(targetSubVolume != null)
					{
						float xCoor = targetVolume.getXCoor(col);
						byte[] subVolumeColPlane = computeSubVolumePlane(XDIR, xCoor, nInputRows * nInputSlices, zDomain);

						for (row = 0; row < nInputRows; row++)
						{
							if (row % 10 == 0)
							{
								progressBarValue += nInputSlices * Math.min(10, (nInputRows - row));
								panel.setValue(progressBarValue);
							}
							for (slice = 0; slice < nInputSlices; slice++)
							{
								if (subVolumeColPlane[n] == 0)
									planeData[n++] = (byte)-1;
								else
								{
									result = crossCorrelate(row, col, slice, parameters, dirNo);
									planeData[n++] = doubleToByte(result);
								}
							}
						}
					}
					else
					{
						for (row = 0; row < nInputRows; row++)
						{
							if (row % 10 == 0)
							{
								progressBarValue += nInputSlices * Math.min(10, (nInputRows - row));
								panel.setValue(progressBarValue);
							}
							for (slice = 0; slice < nInputSlices; slice++)
							{
								result = crossCorrelate(row, col, slice, parameters, dirNo);
								planeData[n++] = doubleToByte(result);
							}
						}
					}
				}
				break;

			case YDIR:
				progressBarAdd = nInputCols * nInputSlices;
				progressBarMax += progressBarAdd;
				panel.setMaximum(progressBarMax);
				n = 0;
				row = nPlane;
				if (row < targetRowMin || row > targetRowMax)
				{
					for (col = 0; col < nInputCols; col++)
					{
						if (col % 10 == 0)
						{
							progressBarValue += nInputSlices * Math.min(10, (nInputCols - col));
							panel.setValue(progressBarValue);
						}
						for (slice = 0; slice < nInputSlices; slice++)
						{
							planeData[n++] = (byte)-1;
						}
					}
				}
				else
				{
					if(targetSubVolume != null)
					{
						float yCoor = targetVolume.getYCoor(row);
						byte[] subVolumeRowPlane = computeSubVolumePlane(YDIR, yCoor, nInputCols * nInputSlices, zDomain);

						for (col = 0; col < nInputCols; col++)
						{
							if (col % 10 == 0)
							{
								progressBarValue += nInputSlices * Math.min(10, (nInputCols - col));
								panel.setValue(progressBarValue);
							}
							for (slice = 0; slice < nInputSlices; slice++)
							{
								if (subVolumeRowPlane[n] == 0)
									planeData[n++] = (byte)-1;
								else
								{
									result = crossCorrelate(nPlane, col, slice, parameters, dirNo);
									planeData[n++] = doubleToByte(result);
								}
							}
						}
					}
					else
					{
						for (col = 0; col < nInputCols; col++)
						{
							if (col % 10 == 0)
							{
								progressBarValue += nInputSlices * Math.min(10, (nInputCols - col));
								panel.setValue(progressBarValue);
							}
							for (slice = 0; slice < nInputSlices; slice++)
							{
								result = crossCorrelate(nPlane, col, slice, parameters, dirNo);
								planeData[n++] = doubleToByte(result);
							}
						}
					}
				}
				break;
			case ZDIR:
				progressBarAdd = nInputCols * nInputRows;
				progressBarMax += progressBarAdd;
				panel.setMaximum(progressBarMax);
				n = 0;
				slice = nPlane;
				if (nPlane < targetSliceMin || nPlane > targetSliceMax)
				{
					for (row = 0; row < nInputRows; row++)
					{
						if (row % 10 == 0)
						{
							progressBarValue += nInputCols * Math.min(10, (nInputRows - row));
							panel.setValue(progressBarValue);
						}
						for (col = 0; col < nInputCols; col++)
						{
							planeData[n++] = (byte)-1;
						}
					}
				}
				else
				{
					if(targetSubVolume != null)
					{
						float zCoor = targetVolume.getZCoor(slice);
						byte[] subVolumeSlicePlane = computeSubVolumePlane(ZDIR, zCoor, nInputRows * nInputCols, zDomain);

						for (row = 0; row < nInputRows; row++)
						{
							if (row % 10 == 0)
							{
								progressBarValue += nInputCols * Math.min(10, (nInputRows - row));
								panel.setValue(progressBarValue);
							}
							for (col = 0; col < nInputCols; col++)
							{
								if (subVolumeSlicePlane[n] == 0)
									planeData[n++] = (byte)-1;
								else
								{
									result = crossCorrelate(row, col, slice, parameters, dirNo);
									planeData[n++] = doubleToByte(result);
								}
							}
						}
					}
					else
					{
						for (row = 0; row < nInputRows; row++)
						{
							if (row % 10 == 0)
							{
								progressBarValue += nInputCols * Math.min(10, (nInputRows - row));
								panel.setValue(progressBarValue);
							}
							for (col = 0; col < nInputCols; col++)
							{
								result = crossCorrelate(row, col, nPlane, parameters, YDIR);
								planeData[n++] = doubleToByte(result);
							}
						}
					}
				}
				break;
		}
		progressBarMax -= progressBarAdd;
		panel.setMaximum(progressBarMax);
		progressBarValue -= progressBarAdd;
		panel.setValue(progressBarValue);

		if (runTimer)
		{
			timer.stopPrint("dry run real for dirNo: " + dirNo + " plane " + nPlane);
		}
		targetVolume.printMemorySummary();
		parameters.printSummary("Dry run for dirNo " + dirNo + " plane " + nPlane + ". Processed ");
//        model.setCurrentObjectDisplayAndToolbar(attributeVolume);
//        model.win3d.getCursor3d().clearTextureDisplays();
//		model.win3dDisplay();
	}

	static private final byte doubleToByte(double result)
	{
		return StsMath.unsignedIntToUnsignedByte( (int) ( (result + 1.0) * 127));
	}

	private void updateDryRunPanel()
	{
		try
		{
//            DateFormat dateFormat = DateFormat.getDateTimeInstance();
//            panel.appendLine(dateFormat.format(new Date(System.currentTimeMillis())) + " :Processing trace #" + value + " of " + (nValues + 1) + " traces.");
//            panel.setProgressPercent( (int) ( ( (float) value / (float) (nValues + 1)) * 100));
			panel.setValue(panelValue);
		}
		catch (Exception e)
		{}
	}

	public boolean fullRun()
	{
		parameters = null;
		createVolume(null);
		if(parameters != null) parameters.printSummary("Full run. Processed ");
		return !canceled;
	}

	static int panelValue = 0;

	private void updateFullRunPanel()
	{
		panel.setValue(panelValue);
	/*
		try
		{
			SwingUtilities.invokeAndWait(new Runnable()
			{
				public void run()
				{
					{
						DateFormat dateFormat = DateFormat.getDateTimeInstance();
						panel.appendLine(dateFormat.format(new Date(System.currentTimeMillis())) +
										 " :Processing inline " + panelValue + " of " + nTargetRows + " inlines.");
						panel.setProgressPercent( (int) ( ( (float) panelValue / (float) (nTargetRows)) * 100));
						if (blockDone)
						{
							panel.appendLine(dateFormat.format(new Date(System.currentTimeMillis())) +
											 " :Writing block to disk, may take a few minutes...");
						}
					}
				}
			});
		}
		catch (Exception e)
		{}
	*/
    }

	public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers)
	{
		int nTrace = -1;
		int row = -1;
		int col = -1;
		int slice = -1;
		blockDone = false;
		int count = 0;

		try
		{
			if (parameters == null)
			{
				parameters = new RunParameters();
			}
			if (runTimer)
			{
				timer.start();
			}

			if (nInputBlockFirstRow < targetRowMin || nInputBlockLastRow > targetRowMax)
			{
				if (transparentRowPlane == null) createTransparentRowPlane();
				if(nullFloatRowPlane == null) createNullFloatRowPlane();
			}
			if (0 < targetColMin || nInputCols - 1 > targetColMax)
			{
				if (transparentTrace == null) createTransparentTrace();
				if (nullTrace == null) createNullFloatTrace();
			}
            for (row = nInputBlockFirstRow; row <= nInputBlockLastRow; row++)
			{
				if(isCanceled()) return false;

				if (row < targetRowMin || row > targetRowMax)
				{
					outputFloatBuffer.put(nullFloatRowPlane);
					count += nInputSlices * nInputCols;
					continue;
				}

				int targetRowStart = Math.max(row - sourceRowCenter, targetRowMin);
				int targetRowEnd = Math.min(row + sourceRowCenter, targetRowMax);
//				targetVolume.checkAllocateRowPlanes(targetRowStart, targetRowEnd, sourceRowMin, sourceRowMax);
//				if (isComplex)
//					targetVolumeH.checkAllocateRowPlanes(targetRowStart, targetRowEnd, sourceRowMin, sourceRowMax);
				panelValue = row + 1 - targetRowMin;
				updateFullRunPanel();
                byte zDomain = sourceVolume.getZDomain();
				if (targetSubVolume != null)
				{
					for (col = 0; col < targetColMin; col++)
					{
						count += nInputSlices;
					}

                    for (col = 0; col < targetColMin; col++)
                        outputFloatBuffer.put(nullTrace);

					float yCoor = targetVolume.getYCoor(row);
					byte[] subVolumeRowPlane = computeSubVolumePlane(YDIR, yCoor, nInputCols * nInputSlices, zDomain);

					int n = targetColMin * nInputSlices;
					for (col = targetColMin; col <= targetColMax; col++)
					{
						for (slice = 0; slice < nInputSlices; slice++, n++)
						{
							if (subVolumeRowPlane[n] == 0)
							{
								outputFloatBuffer.put(StsParameters.nullValue);
								count++;
							}
							else
							{
								double value = crossCorrelate(row, col, slice, parameters, YDIR);
								int unsignedInt = (int) ( (value + 1.0) * 127);
								outputVolume.accumulateHistogram(unsignedInt);
								outputFloatBuffer.put( (float) value);
								count++;
							}
						}
					}
					for (col = targetColMax + 1; col < nInputCols; col++)
					{
						count += nInputSlices;
					}
                    for (col = targetColMax + 1; col < nInputCols; col++)
                        outputFloatBuffer.put(nullTrace);
				}
				else
				{
					for (col = 0; col < nInputCols; col++)
					{
						for (slice = 0; slice < nInputSlices; slice++)
						{
							double value = crossCorrelate(row, col, slice, parameters, YDIR);
							int unsignedInt = (int) ( (value + 1.0) * 127);
							outputVolume.accumulateHistogram(unsignedInt);
							outputFloatBuffer.put( (float) value);
							count++;
						}
					}
				}
			}
			blockDone = true;
 //           panel.incrementCount();
			updateFullRunPanel();

			if (runTimer)
			{
				timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");
			}
			model.viewObjectChanged(this, outputVolume);
			model.win3dDisplay();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsAnalogueVolumeConstructor.doProcessBlock() failed for " + "row: " + row + "col: " + col +
										 " slice: " + slice + "sampleCount: " + count + "\n" + " nInputBlockFirstRow: " +
                    nInputBlockFirstRow + " nInputBlockLastRow: " + nInputBlockLastRow + " targetColMin: " +
										 targetColMin + " targetColMax: " + targetColMax + " targetSliceMin: " +
										 targetSliceMin + " targetSliceMax: " + targetSliceMax, e, StsException.WARNING);
			return false;
		}
	}

	private byte[] computeSubVolumePlane(int dir, float coor, int nBytes, byte zDomainData)
	{
		byte[] subVolumePlane = new byte[nBytes];
		targetSubVolume.addUnion(subVolumePlane, dir, coor, outputVolume, zDomainData);
		return subVolumePlane;
	}

	private void createTransparentRowPlane()
	{
		transparentRowPlane = new byte[nInputCols * nInputSlices];
		for (int n = 0; n < nInputCols * nInputSlices; n++)
		{
			transparentRowPlane[n] = (byte) - 1;
		}
	}

	private void createNullFloatRowPlane()
	{
		nullFloatRowPlane = new float[nInputCols * nInputSlices];
		for (int n = 0; n < nInputCols * nInputSlices; n++)
		{
			nullFloatRowPlane[n] = StsParameters.nullValue;
		}
	}

	private void createTransparentTrace()
	{
		transparentTrace = new byte[nInputSlices];
		for (int n = 0; n < nInputSlices; n++)
		{
			transparentTrace[n] = (byte) - 1;
		}
	}

	private void createNullFloatTrace()
	{
		nullTrace = new float[nInputSlices];
		for (int n = 0; n < nInputSlices; n++)
		{
			nullTrace[n] = StsParameters.nullValue;
		}
	}

	private void fillSourceCube(StsSeismicVolume volume, double[][][] sourceCube)
	{
		for (int row = sourceRowMin, r = 0; row <= sourceRowMax; row++, r++)
		{
			for (int col = sourceColMin, c = 0; col <= sourceColMax; col++, c++)
			{
				volume.getTraceValues(row, col, sourceSliceMin, sourceSliceMax, YDIR, useByteCubes, sourceCube[r][c]);
			}
		}
	}

	private double crossCorrelate(int row, int col, int slice, RunParameters parameters, int dirNo)
	{
		int sourceSlice = -1, targetSlice = -1;
		int sourceColStart, sourceColEnd;
		int sourceSliceStart, sourceSliceEnd;
		int targetSliceStart, targetSliceEnd;
		SpiralIterator spiralIterator = null;
		try
		{
			targetSliceStart = slice - sourceSliceCenter;
			targetSliceEnd = targetSliceStart + nSourceSlices - 1;

			if (targetSliceStart < targetSliceMin)
			{
				sourceSliceStart = targetSliceMin - targetSliceStart;
				targetSliceStart = targetSliceMin;
			}
			else
			{
				sourceSliceStart = 0;
			}
			if (targetSliceEnd > targetSliceMax)
			{
				sourceSliceEnd = nSourceSlices - 1 - (targetSliceEnd - targetSliceMax);
				targetSliceEnd = targetSliceMax;
			}
			else
			{
				sourceSliceEnd = nSourceSlices - 1;

			}
			int nSlices = sourceSliceEnd - sourceSliceStart + 1;

//        int nInputRows = sourceRowEnd - sourceRowStart + 1;
//        int nInputCols = sourceColEnd - sourceColStart + 1;
//        int nTraces = nInputRows*nInputCols;
			double cor = 0.0;
			int nTraces = 0;
//            mainDebug = (slice == 10);
			if (debug)
			{
				System.out.println("Block center at target row: " + row + " col: " + col + " slice: " + slice);
			}
			spiralIterator = new SpiralIterator(row, col, sourceRowCenter, sourceColCenter);
			int nSkipTraces = Math.max(spiralIterator.nTargetTraces / nToSkip, minNToSkip);
			int nTracesSkipped = 0;
			double traceCor = 0.0;
			boolean hitThreshold = false;
			while (spiralIterator.hasNext())
			{
				parameters.nTotalTraces++;
				if (!spiralIterator.isOK())
				{
					continue;
				}
				nTraces++;
				if (hitThreshold)
				{
					nTracesSkipped++;
					if (nTracesSkipped == nSkipTraces)
					{
						nTracesSkipped = 0;
					}
					if (nTracesSkipped > 0)
					{
						if (debug)
						{
							spiralIterator.nSpiralTracesSkipped++;
						}
						parameters.nTotalTracesSkipped++;
						cor += traceCor;
						continue;
					}
				}

				if (debug)
				{
					spiralIterator.nSpiralTracesComputed++;

				}
				targetVolume.getTraceValues(spiralIterator.targetRow, spiralIterator.targetCol, targetSliceStart,
									   targetSliceEnd, dirNo, useByteCubes, targetValues);
				double[] sourceValues = sourceCube[spiralIterator.sourceRow][spiralIterator.sourceCol];

				double[] sourceValuesH = null;
				if (isComplex)
				{
					targetVolumeH.getTraceValues(spiralIterator.targetRow, spiralIterator.targetCol, targetSliceStart,
											targetSliceEnd, dirNo, useByteCubes, targetValuesH);
					sourceValuesH = sourceCubeH[spiralIterator.sourceRow][spiralIterator.sourceCol];
				}

				traceCor = 0.0;
				double sourceCovar = 0.0;
				double targetCovar = 0.0;
				for (sourceSlice = sourceSliceStart, targetSlice = 0; sourceSlice <= sourceSliceEnd; sourceSlice++,
					 targetSlice++)
				{
					double sourceValue = sourceValues[sourceSlice];
					double targetValue = targetValues[targetSlice];
					float coef = StsMath.simpsonCoefs(targetSlice, nSlices);
					if (!isComplex)
					{
						sourceCovar += coef * sourceValue * sourceValue;
						targetCovar += coef * targetValue * targetValue;
						traceCor += coef * targetValue * sourceValue;
					}
					else
					{
						double sourceValueH = sourceValuesH[sourceSlice];
						double targetValueH = targetValuesH[targetSlice];
						sourceCovar += coef * (sourceValue * sourceValue + sourceValueH * sourceValueH);
						targetCovar += coef * (targetValue * targetValue + targetValueH * targetValueH);
						traceCor += coef * (targetValue * sourceValue + targetValueH * sourceValueH);
					}
				}
				sourceCovar /= nSlices;
				targetCovar /= nSlices;
				traceCor /= nSlices * Math.sqrt(sourceCovar * targetCovar);
				cor += traceCor;
				if (debug)
				{
					spiralIterator.spiralDebugCheck(traceCor);
				}
				if (!hitThreshold && spiralIterator.hitThreshold(traceCor, cor, nTraces))
				{
					hitThreshold = true;
				}
			}
			if (debug)
			{
				spiralIterator.blockDebugCheck();
			}
			if (nTraces == 0)
			{
				return 0.0f;
			}
			else
			{
				return cor / nTraces;
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsAnalogueVolumeConstructor.crossCorrelate() failed.", e, StsException.WARNING);
			return 0.0f;
		}
	}

	class RunParameters
	{
		int nTotalTraces = 0;
		int nTotalTracesSkipped = 0;
		int nTracesProcessed = 0;

		RunParameters()
		{
		}

		void printSummary(String s)
		{
			nTracesProcessed = nTotalTraces - nTotalTracesSkipped;
			if (nTotalTraces > 0)
			{
				System.out.println(s + " " + nTracesProcessed + " out of " + nTotalTraces + " (" +
								   100 * nTracesProcessed / nTotalTraces + "%).");
			}
		}
	}

	class SpiralIterator
	{
		int targetRow, targetCol;
		int sourceRow, sourceCol;
		boolean spiralEnded = true;
		int spiralSize = -1;
		int spiralSide = 3;
		int sideIndex = 0;
		int spiral = 0;
		double spiralCorrelSum = 0;
		int nSpiralTraces = 0;
		int[] rowInc = new int[]
			{0, 1, 0, -1};
		int[] colInc = new int[]
			{1, 0, -1, 0};
		int targetRowStart;
		int targetRowEnd;
		int targetColStart;
		int targetColEnd;
		int nTargetTraces;
		int nCountedTargetTraces = 0;
		boolean spiralComplete = false;
		boolean aboveThreshold = false;
		int nSpiralTracesSkipped;
		int nSpiralTracesComputed;

		SpiralIterator(int row, int col, int sourceRowCenter, int sourceColCenter)
		{
			targetRow = row;
			targetCol = col;
			sourceRow = sourceRowCenter;
			sourceCol = sourceColCenter;
			targetRowStart = Math.max(row - sourceRowCenter, targetRowMin);
			targetRowEnd = Math.min(row + sourceRowCenter, targetRowMax);
			targetColStart = Math.max(col - sourceColCenter, targetColMin);
			targetColEnd = Math.min(col + sourceColCenter, targetColMax);
			nTargetTraces = (targetRowEnd - targetRowStart + 1) * (targetColEnd - targetColStart + 1);
		}

		public boolean hasNext()
		{
			if (spiralSize == -1)
			{
				spiralSize = 1;
				spiralComplete = true;
				return true;
			}
			if (spiralComplete)
			{
				if (debug)
				{
					System.out.println(" nSpiralTracesComputed: " + nSpiralTracesComputed + " nSpiralTracesSkipped: " +
									   nSpiralTracesSkipped);

				}
				if (++spiral > maxNSpirals)
				{
					return false;
				}
				spiralSide = 0;
				sideIndex = 1;

				if (spiralSize == 1)
				{
					sourceRow--;
					targetRow--;
				}
				else
				{
					sourceRow -= 2;
					targetRow -= 2;
				}
				sourceCol--;
				targetCol--;

				spiralSize += 2;
				spiralComplete = false;

				spiralCorrelSum = 0;
				nSpiralTraces = 0;
				aboveThreshold = false;
				nSpiralTracesSkipped = 0;
				nSpiralTracesComputed = 0;

				if (debug)
				{
					printSpiralStart();
				}
				return true;
			}
			else
			{
				sourceRow += rowInc[spiralSide];
				sourceCol += colInc[spiralSide];
				targetRow += rowInc[spiralSide];
				targetCol += colInc[spiralSide];
				sideIndex++;

				if (sideIndex == spiralSize)
				{
					sideIndex = 1;
					spiralSide++;
				}
				else if (sideIndex == spiralSize - 1 && spiralSide == 3)
				{
					spiralComplete = true;

				}
			}
			return true;
		}

		void printSpiralStart()
		{
			System.out.println("    Spiral start at targetRow: " + targetRow + " targetCol: " + targetCol);
		}

		public boolean isOK()
		{
			boolean traceOK = targetRow >= targetRowStart && targetRow <= targetRowEnd && targetCol >= targetColStart &&
				targetCol <= targetColEnd && sourceRow >= 0 && sourceRow < nSourceRows && sourceCol >= 0 &&
				sourceCol < nSourceCols;
			if (traceOK)
			{
				nCountedTargetTraces++;
			}
			return traceOK;
		}

		public boolean hitThreshold(double traceCorrel, double sumCorrel, int nTraces)
		{
			switch (spiralMethod)
			{
				case StsAnalogueCubeWizard.RUNNING_AVERAGE:
					if (debug && sumCorrel / nTraces < correlThreshold)
					{
						System.out.println("    hit RUNNING_AVG threshold for " + nTraces +
										   " total traces at avg value " + sumCorrel / nTraces);
					}
					return sumCorrel / nTraces < correlThreshold;
				case StsAnalogueCubeWizard.SPIRAL_AVERAGE:
					spiralCorrelSum += traceCorrel;
					nSpiralTraces++;
					if (!spiralComplete)
					{
						return false;
					}
					boolean hitThreshold = spiralCorrelSum / nSpiralTraces < correlThreshold;
					if (hitThreshold && debug)
					{
						System.out.println("    hit SPIRAL_AVG threshold for " + nSpiralTraces +
										   " traces at avg value " + spiralCorrelSum / nSpiralTraces);
					}
					return hitThreshold;
				case StsAnalogueCubeWizard.SPIRAL_MAXIMUM:
					if (!aboveThreshold)
					{
						aboveThreshold = traceCorrel >= correlThreshold;
					}
					if (debug && spiralComplete && !aboveThreshold)
					{
						System.out.println("    hit SPIRAL_MAX threshold for " + nSpiralTraces +
										   " traces at trace value " + traceCorrel);
					}
					return spiralComplete && !aboveThreshold;
				default:
					return false;
			}
		}

		void spiralDebugCheck(double traceCor)
		{
			System.out.println("        targetRow: " + targetRow + " targetCol: " + targetCol + " traceCor: " +
							   traceCor);
		}

		void blockDebugCheck()
		{
			if (nCountedTargetTraces != nTargetTraces)
			{
				System.out.println("nCountedTargetTraces: " + nCountedTargetTraces + " and nTargetTraces: " +
								   nTargetTraces + " don't agree.");
			}
		}
	}

	public void setCorrelationThreshold(float correlThreshold)
	{
		this.correlThreshold = correlThreshold;
	}

	public StsBoxSubVolume getSourceSubVolume()
	{
		return sourceSubVolume;
	}

	public StsSeismicVolume getSourceVolume()
	{
		return sourceVolume;
	}

	public StsSeismicVolume getSourceVolumeH()
	{
		return sourceVolumeH;
	}

	public StsSeismicVolume getTargetVolume()
	{
		return targetVolume;
	}

	public StsSeismicVolume getTargetVolumeH()
	{
		return targetVolumeH;
	}

	public StsSubVolume getTargetSubVolume()
	{
		return targetSubVolume;
	}

	public boolean isDataFloat()
	{
		return isDataFloat;
	}

	public boolean isIsComplex()
	{
		return isComplex;
	}
	public void useByteCubes(boolean b)
	{
		useByteCubes = b;
	}
}
