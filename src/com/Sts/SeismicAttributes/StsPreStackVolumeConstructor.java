package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

/** This is just a workhorse class as it's created and used just for computing a preStack attribute volume */
public class StsPreStackVolumeConstructor extends StsVolumeConstructor
{
    protected StsPreStackLineSet3d lineSet3d;
    int windowWidth = 15;
    float maxAmplitude = 1.0f;
    boolean applyAGCPoststack = false;
    static final boolean printIntervalTimes = true;

    public boolean isOutputDataFloat() { return true; }
    public boolean initializeBlockInput() { return true; }
 //   public void finalizeVolumeOutput() { }

    public StsPreStackVolumeConstructor(StsModel model, StsPreStackLineSet3d seismic, String typename, StsPreStackVolume volume, StsProgressDialog dialog)
    {
        initialize(model, seismic, typename, volume, dialog.getProgressPanel());
        this.dialog = dialog;
        createOutputVolume();
    }

    public StsPreStackVolumeConstructor(StsModel model, StsPreStackLineSet3d seismic, String typename, StsPreStackVolume volume, StsProgressPanel panel)
    {
        initialize(model, seismic, typename, volume, panel);
        createOutputVolume();
    }

    private void  initialize(StsModel model, StsPreStackLineSet3d seismic, String typename, StsPreStackVolume volume, StsProgressPanel panel)
    {
        this.model = model;
        this.lineSet3d = seismic;
        outputVolume = volume;
        this.panel = panel;
        volumeName = typename;
        lineSet3d.checkTransparentTrace(nInputSlices);
        outputPosition = 0;
        StsSeismicTimer.printIntervalElapsedTimes(printIntervalTimes);
    }

    public void createOutputVolume()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                if(runCreateOutputVolume())  ((StsPreStackVolume)outputVolume).validateAllPlanes();
            }
        };
        StsToolkit.runRunnable(runnable);
    }

    public void initializeVolumeInput()
    {
        super.initializeVolumeInput();
        outputFloatBuffer = outputVolume.createMappedFloatRowBuffer();
        outputPosition = 0;
        initializeAGC();
    }

    private void initializeAGC()
    {
        StsAGCPreStackProperties agcProperties = lineSet3d.agcProperties;
        applyAGCPoststack = agcProperties.getApplyAGCPoststack();
        if(applyAGCPoststack)
        {
            windowWidth = agcProperties.getWindowWidth(StsFilterProperties.POSTSTACK);
            maxAmplitude = 1.0f;
            //maxAmplitude = Math.max(Math.abs(dataMin), Math.abs(dataMax));
        }
    }


    static public void computeVolumeWithDialog(StsModel model, StsPreStackLineSet3d seismic, byte type, StsPreStackVolume volume)
    {
        String typename = StsPreStackVolume.typenames[type];
        StsProgressDialog dialog = new StsProgressDialog(model.win3d, "Compute Seismic PostStack3d " + typename, false);
        new StsPreStackVolumeConstructor(model, seismic, typename, volume, dialog);
    }

    protected String getFullStemname(StsModel model)
    {
        return model.getProject().getName() + "." + lineSet3d.stemname;
    }

    public boolean allocateMemory()
    {
        memoryAllocation = StsMemAllocVolumeProcess.constructor(outputVolume);
        if(memoryAllocation == null) return false;
        nOutputSamplesPerRow = memoryAllocation.nOutputSamplesPerRow;
        return true;
    }
	public  boolean processBlockInput(int nBlock, String mode ) { return processBlockInput(nBlock); }
    public boolean processBlockInput(int nBlock)
	{
		if (runTimer) timer.getBlockInputTimer.start();

        try
        {
            outputFloatBuffer.map(outputPosition, nOutputSamplesPerInputBlock);
            outputPosition += nOutputSamplesPerInputBlock;

            int i = 0;
            StsSuperGather gather = new StsSuperGather(model, lineSet3d);
            StsPreStackVolume preStackVolume = (StsPreStackVolume)outputVolume;
            boolean isVolumeScaling = outputVolume.isVolumeScaling();
            for (int row = nInputBlockFirstRow; row <= nInputBlockLastRow; row++)
            {
                for (int col = 0; col < nInputCols; col++)
                {
                    if (!gather.initializeSuperGather(row, col))
                    {
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        outputFloatBuffer.put(lineSet3d.floatTransparentTrace);
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                        continue;
                    }
                    float[] floatData = preStackVolume.computeTrace(gather, nInputSlices);
                    // float[] floatData = preStackVolume.computeTraceAGC(gather, nInputSlices, windowWidth, maxAmplitude);
                    if(floatData != null)
                    {
                        StsMath.scaleNormalize(floatData);
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        outputFloatBuffer.put(floatData);
                        // if(isVolumeScaling) adjustOutputDataRange(floatData);
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    else
                    {
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.start();
    //                    outputRowByteBuffer.put(segyVolume.byteTransparentTrace);
    //                    if (segyVolume.isDataFloat)
                        outputFloatBuffer.put(lineSet3d.floatTransparentTrace);
                        if (runTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                }
 //               panel.incrementCount();
            }
            outputFloatBuffer.clear();
            if (runTimer)
            {
                timer.getBlockInputTimer.stopAccumulateIncrementCountPrintInterval("process " + nInputBlockTraces + " traces for block " + nBlock + ":");
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "doProcessBlock", e);
            return false;
        }
    }

    protected void initializeVolumeOutput()
    {
        super.initializeVolumeOutput();
        outputVolume.setDataRange(dataMin, dataMax);
    }

     public void finalizeVolumeOutput()
     {
         if (runTimer) StsSeismicTimer.getMapBlockTimer.start();
         finalizeVolumeOutputBuffers();
         if (runTimer) StsSeismicTimer.getMapBlockTimer.stopAccumulateIncrementCount();
         outputVolume.calculateHistogram();
         outputVolume.writeHeaderFile();
         model.getProject().runCompleteLoading();
         StsSeismicTimer.printIntervalElapsedTimes(!printIntervalTimes);
     }
 }
