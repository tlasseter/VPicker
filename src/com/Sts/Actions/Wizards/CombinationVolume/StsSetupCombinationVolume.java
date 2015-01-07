
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsSetupCombinationVolume extends StsWizardStep implements Runnable
{
    public StsStatusPanel panel;
    private StsHeaderPanel header;
    private StsToolbar subVolumeToolbar = null;
    private StsCombinationVolumeWizard wizard = null;

    RandomAccessFile rFile, sFile, cFile;
    float min, max;
    transient private int dataCnt[] = new int[255];
    transient private int ttlHistogramSamples = 0;
    public float[] dataHist = new float[255]; // Histogram of the data distribution
    StsSegyVolume svol = new StsSegyVolume();
    boolean autoScale = true;

//	final static public String volumeFilePrefix = "seis.vol.";

    public StsSetupCombinationVolume(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        panel = (StsStatusPanel) getContainer();
        this.wizard = (StsCombinationVolumeWizard)wizard;

        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("PostStack3d Combination Definition");
        header.setSubtitle("Produce Combination PostStack3d(s)");
        header.setInfoText(wizardDialog,"(1) Press the Cancel Button to stop the processing.\n" +
                                   "(2) Press the Finish Button once processing is complete.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/VolumeCombination.html");
    }

    public boolean start()
    {
        panel.setTitle("Creating Combination PostStack3d");
        run();
        return true;
    }

    public void run()
    {
        try
        {
            panel.setMaximum(100); // 10 is the scaling: max = n*scaling
            disablePrevious();
            disableFinish();
            panel.setProgress(0.0f);
//            Thread.currentThread().sleep(10);

            checkAddToolbar();
            int type = wizard.getType();
            if(type == wizard.BOX_SET)
            {
                panel.setText("Creating Box Set Summation PostStack3d");
                createBoxSumVolume();
            }
            else if(type == wizard.CUMMULATIVE)
            {
                panel.setText("Creating Cumulative Difference PostStack3d");
                createCummulativeVolume();
            }
            else if(type == wizard.SV_CUMMULATIVE)
            {
                panel.setText("Creating Cumulative Sub-PostStack3d PostStack3d");
                createCummulativeSubVolumeVolume();
            }
            else if(type == wizard.DIFFERENCE)
            {
                panel.setText("Creating Difference PostStack3d");
                createDifferenceVolume();
            }
            else if(type == wizard.MATH)
            {
                panel.setText("Creating Math PostStack3d");
                createMathVolume();
            }
 //           Thread.currentThread().sleep(200);
            panel.setText("Setup Complete");
            panel.setProgress(100.0f);
            success = true;
        }
        catch(Exception e)
        {
            success = false;
        }
        finally
        {
            disableCancel();
            wizard.enableFinish();
            wizard.finish();
        }
    }

    public boolean end()
    {
        return true;
    }

    private void createBoxSumVolume()
    {
        int j, k, h, i, nPoints;
        byte[][] planes = new byte[2][];
        byte[] planeData;
        boolean found = false;
        float[] floatValues;
        int svVal = 1;

        StsBoxSetSubVolume boxSet = wizard.getBoxSet();
        StsObjectRefList boxes = boxSet.getBoxes();
        String name = wizard.getVolumeName();

        // Does the user want to keep the data inside or outside the subvolumes
        if(wizard.isInclusive()) svVal = 1;
        else svVal = 0;

        StsSeismicVolume[] vols = new StsSeismicVolume[boxes.getSize()];
        try
        {
            min = StsParameters.largeFloat;
            max = -StsParameters.largeFloat;
            int numVols = 0;
            for(i=0; i<boxes.getSize(); i++)
            {
                StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(i);
                StsSeismicVolume vol = box.getVolume();
                found = false;
                for(j = 0; j<vols.length; j++)
                {
                    if(vols[j] == null)
                        break;
                    if(vols[j] == vol)
                    {
                        found = true;
                        break;
                    }
                }
                if(found == false)
                {
                    vols[j] = vol;
                    if(vol.getDataMin() < min)
                        min = vol.getDataMin();
                    if(vol.getDataMax() > max)
                        max = vol.getDataMax();
                    numVols++;
                }
            }
            StsMath.trimArray(vols, numVols);
            getOutputFiles(name);
            clearHistogram();

            panel.setText("Processing Slice Data");
            Main.logUsage();

            nPoints = vols[0].nCols * vols[0].nRows;
            floatValues = new float[nPoints];
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                byte zDomain = vols[0].getZDomain();
                for(h=0; h<numVols; h++)
                {
                    planeData = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));

                    byte[] subVolumePlane = new byte[nPoints];
                    boxSet.addUnion(subVolumePlane, StsCursor3d.ZDIR, vols[h].getZCoor(i), vols[h], true, zDomain);

                    // Merge the plane and the subVolume mask
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(subVolumePlane[n] == (byte) svVal)
                            floatValues[n] = vols[h].getScaledValue(planeData[n]);
                    }
                }
                panel.setProgress((float)((float)i/(float)vols[0].nSlices) * 33.0f);
                sFile.write(computeUnscaledByteValues(floatValues, min, max));
            }

            panel.setText("Processing Crossline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nRows;
            floatValues = new float[nPoints];
            for(i = 0; i <= vols[0].nCols; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                byte zDomain = vols[0].getZDomain();
                for(h = 0; h<numVols; h++)
                {
                    planeData = vols[h].readBytePlaneData(StsCursor3d.XDIR, vols[h].getXCoor(i));

                    byte[] subVolumePlane = new byte[nPoints];
                    boxSet.addUnion(subVolumePlane,StsCursor3d.XDIR, vols[h].getXCoor(i), vols[h], true, zDomain);

                    for (int n = 0; n < nPoints; n++)
                    {
                        if (subVolumePlane[n] == svVal)
                            floatValues[n] = vols[h].getScaledValue(planeData[n]);
                    }
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nCols) * 33.0f) + 33.0f);
                cFile.write(computeUnscaledByteValues(floatValues, min, max));
            }

            panel.setText("Processing Inline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nCols;
            floatValues = new float[nPoints];
            for(i = 0; i <= vols[0].nRows; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                byte zDomain = vols[0].getZDomain();
                for(h = 0; h<numVols; h++)
                {
                    planeData = vols[h].readBytePlaneData(StsCursor3d.YDIR, vols[h].getYCoor(i));

                    byte[] subVolumePlane = new byte[nPoints];
                    boxSet.addUnion(subVolumePlane, StsCursor3d.YDIR, vols[h].getYCoor(i), vols[h], true, zDomain);

                    for (int n = 0; n < nPoints; n++)
                    {
                        if (subVolumePlane[n] == svVal)
                            floatValues[n] = vols[h].getScaledValue(planeData[n]);
                    }
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nRows) * 33.0f) + 67.0f);
                rFile.write(computeUnscaledByteValues(floatValues, min, max));
            }
            calculateHistogram();

            StsSegyVolume newSegY = buildSegYObject(name, vols[0]);
            StsParameterFile.writeObjectFields(model.getProject().getRootDirString() + "seis3d.txt." + name, newSegY, StsSegyVolume.class, StsBoundingBox.class);

            closeOutputFiles();

            // Load the seismic volume
            StsSeismicVolume seismicVolume = StsSeismicVolume.checkLoadFromStemname(model, model.getProject().getRootDirString(), name , true);
            if(seismicVolume == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING,"Failed to add volume to project: " + name);
                return;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("createDifferenceVolume failed.", e, StsException.WARNING);
        }
    }

    public float[] calcDataRange()
    {
        int j, k, h, i;
        byte[][] planes = new byte[2][];
        byte[] plane;
        int nPoints;
        float[] floatValues;
        float[] range = new float[2];

        StsSeismicVolume[] vols = wizard.cummDiffVolume.getVolumes();
        range[0] = StsParameters.largeFloat;
        range[1] = -StsParameters.largeFloat;

        nPoints = vols[0].nCols * vols[0].nRows;
        floatValues = new float[nPoints];

        // Process 10% of traces
        for (i = vols[0].nSlices/4; i <= vols[0].nSlices/4 + vols[0].nSlices*.10; i++)
        {
            for(j=0; j<nPoints; j++)
                floatValues[j] = StsParameters.nullValue;
            for(h = 0; h<vols.length-(vols.length%2); h=h+2)
            {
                planes[0] = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                planes[1] = vols[h+1].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                for (int n = 0; n < nPoints; n++)
                {
                    if(!isAValueNull(planes[0][n], planes[1][n]))
                    {
                        float value0 = vols[h].getScaledValue(planes[0][n]);
                        float value1 = vols[h+1].getScaledValue(planes[1][n]);
                        if(floatValues[n] == StsParameters.nullValue)
                            floatValues[n] = value0 - value1;
                        else
                            floatValues[n] = floatValues[n] + (value0 - value1);
                    }
                }
            }
            for(j=0; j<nPoints; j++)
            {
                if((floatValues[j] > max) || (floatValues[j] < min))
                            floatValues[j] = StsParameters.nullValue;
            }
            range = recalcScale(floatValues, range);
        }
        return range;
    }

    private void createCummulativeVolume()
    {
        int j, k, h, i;
        byte[][] planes = new byte[2][];
        byte[] plane;
        int nPoints;
        float[] floatValues;

        StsSeismicVolume[] vols = wizard.cummDiffVolume.getVolumes();
        min = wizard.cummDiffVolume.getDataMin();
        max = wizard.cummDiffVolume.getDataMax();
        String name = wizard.getVolumeName();

        try
        {
            getOutputFiles(name);
            clearHistogram();

            panel.setText("Processing Slice Data");
            Main.logUsage();

            nPoints = vols[0].nCols * vols[0].nRows;
            floatValues = new float[nPoints];
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                for(h = 0; h<vols.length-(vols.length%2); h=h+2)
                {
                    planes[0] = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                    planes[1] = vols[h+1].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(!isAValueNull(planes[0][n], planes[1][n]))
                        {
                            float value0 = vols[h].getScaledValue(planes[0][n]);
                            float value1 = vols[h+1].getScaledValue(planes[1][n]);
                            if(floatValues[n] == StsParameters.nullValue)
                                floatValues[n] = value0 - value1;
                            else
                                floatValues[n] = floatValues[n] + (value0 - value1);
                        }
                    }
                }
                for(j=0; j<nPoints; j++)
                {
                    if((floatValues[j] > max) || (floatValues[j] < min))
                                floatValues[j] = StsParameters.nullValue;
                }
                panel.setProgress((float)((float)i/(float)vols[0].nSlices) * 33.0f);
                sFile.write(computeUnscaledByteValues(floatValues, min, max));
            }

            panel.setText("Processing Crossline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nRows;
            floatValues = new float[nPoints];
            for(i = 0; i <= vols[0].nCols; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                for(h = 0; h<vols.length-(vols.length%2); h=h+2)
                {
                    planes[0] = vols[h].readBytePlaneData(StsCursor3d.XDIR, vols[h].getXCoor(i));
                    planes[1] = vols[h+1].readBytePlaneData(StsCursor3d.XDIR, vols[h].getXCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(!isAValueNull(planes[0][n], planes[1][n]))
                        {
                            float value0 = vols[h].getScaledValue(planes[0][n]);
                            float value1 = vols[h+1].getScaledValue(planes[1][n]);
                            if(floatValues[n] == StsParameters.nullValue)
                                floatValues[n] = value0 - value1;
                            else
                                floatValues[n] = floatValues[n] + (value0 - value1);
                        }
                    }
                }
                for(j=0; j<nPoints; j++)
                {
                    if((floatValues[j] > max) || (floatValues[j] < min))
                                floatValues[j] = StsParameters.nullValue;
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nCols) * 33.0f) + 33.0f);
                cFile.write(computeUnscaledByteValues(floatValues, min, max));
            }

            panel.setText("Processing Inline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nCols;
            floatValues = new float[nPoints];
            for(i = 0; i <= vols[0].nRows; i++)
            {
                for(j=0; j<nPoints; j++)
                    floatValues[j] = StsParameters.nullValue;
                for(h = 0; h<vols.length-(vols.length%2); h=h+2)
                {
                    planes[0] = vols[h].readBytePlaneData(StsCursor3d.YDIR, vols[h].getYCoor(i));
                    planes[1] = vols[h+1].readBytePlaneData(StsCursor3d.YDIR, vols[h+1].getYCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(!isAValueNull(planes[0][n], planes[1][n]))
                        {
                            float value0 = vols[h].getScaledValue(planes[0][n]);
                            float value1 = vols[h+1].getScaledValue(planes[1][n]);
                            if(floatValues[n] == StsParameters.nullValue)
                                floatValues[n] = value0 - value1;
                            else
                                floatValues[n] = floatValues[n] + (value0 - value1);
                        }
                    }
                }
                for(j=0; j<nPoints; j++)
                {
                    if((floatValues[j] > max) || (floatValues[j] < min))
                                floatValues[j] = StsParameters.nullValue;
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nRows) * 33.0f) + 67.0f);
                rFile.write(computeUnscaledByteValues(floatValues, min, max));
            }
            calculateHistogram();

            StsSegyVolume newSegY = buildSegYObject(name, vols[0]);
            StsParameterFile.writeObjectFields(model.getProject().getRootDirString() + "seis3d.txt." + name, newSegY, StsSegyVolume.class, StsBoundingBox.class);

            closeOutputFiles();

            // Load the seismic volume
            StsSeismicVolume seismicVolume = StsSeismicVolume.checkLoadFromStemname(model, model.getProject().getRootDirString(), name, true);
            if(seismicVolume == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add volume to project: " + name);
                return;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("createDifferenceVolume failed.", e, StsException.WARNING);
        }

    }

    private double computeValue(double currentValue, double value)
    {
        switch(wizard.mathVolume.getOperator())
        {
            case StsMathVirtualVolume.ADDITION:
                return currentValue + value;
            case StsMathVirtualVolume.SUBTRACTION:
                return currentValue - value;
            case StsMathVirtualVolume.MULTIPLY:
                return currentValue * value;
            case StsMathVirtualVolume.DIVIDE:
                return currentValue / value;
            case StsMathVirtualVolume.AVERAGE:
                return currentValue + value;
            case StsMathVirtualVolume.MAXIMUM:
                if(currentValue < value)
                    return value;
                else
                    return currentValue;
            case StsMathVirtualVolume.MINIMUM:
                if(currentValue > value)
                    return value;
                else
                    return currentValue;
            default:
                return currentValue;
        }
    }

    private void createMathVolume()
    {
        int j, k, h, i;
        byte[] plane;
        int nPoints;
        float[] floatValues;
        double[] doubleValues;

        StsSeismicVolume[] vols = wizard.mathVolume.getVolumes();

        min = wizard.mathVolume.getDataMin();
        max = wizard.mathVolume.getDataMax();
        String name = wizard.getVolumeName();

        try
        {
            getOutputFiles(name);
            clearHistogram();

            panel.setText("Processing Slice Data");
            Main.logUsage();

            nPoints = vols[0].nCols * vols[0].nRows;
            floatValues = new float[nPoints];
            doubleValues = new double[nPoints];

            // Determine Min and Max for ComboVolume
            float ttlMin = StsParameters.largeFloat;
            float ttlMax = StsParameters.smallFloat;

            // Determine the total minimum and maximum of the volume
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for(j=0; j<nPoints; j++)
                {
                    doubleValues[j] = StsParameters.nullValue;
                    floatValues[j] = StsParameters.nullValue;
                }
                for(h = 0; h<vols.length; h++)
                {
                    min = wizard.mathVolume.getDataMin();
                    if(!vols[h].getIsVisible())
                        continue;

                    plane = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(plane[n] != -1)
                        {
                            double value0 = vols[h].getScaledValue(plane[n]);
                            if((value0 > max) || (value0 < min))
                                continue;

                            if(doubleValues[n] == StsParameters.nullValue)
                                doubleValues[n] = value0;
                            else
                                doubleValues[n] = computeValue(doubleValues[n], value0);
                        }
                    }
                }
                // Calc Min and max
                for(j=0; j<nPoints; j++)
                {
                    if((wizard.mathVolume.getOperator() == StsMathVirtualVolume.AVERAGE) && (doubleValues[j] != StsParameters.nullValue))
                        doubleValues[j] = doubleValues[j]/(double)vols.length;
                    if(doubleValues[j] < Float.MAX_VALUE)
                        floatValues[j] = (float)doubleValues[j];
                    if(floatValues[j] != StsParameters.nullValue)
                    {
                        if (floatValues[j] > ttlMax)
                            ttlMax = floatValues[j];
                        if (floatValues[j] < ttlMin)
                            ttlMin = floatValues[j];
                    }
                }
            }

            if(ttlMax > StsParameters.doubleNullValue)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Maximum value is too large, constrain parameters");
                return;
            }
            // Initialize the pointset
            int volIdx = 0;

            // Output the volumes and pointset
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for(j=0; j<nPoints; j++)
                {
                    doubleValues[j] = StsParameters.nullValue;
                    floatValues[j] = StsParameters.nullValue;
                }
                for(h = 0; h<vols.length; h++)
                {
                    min = wizard.mathVolume.getDataMin();
                    plane = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(plane[n] != -1)
                        {
                            double value0 = vols[h].getScaledValue(plane[n]);
                            if((value0 > max) || (value0 < min))
                                continue;

                            if(doubleValues[n] == StsParameters.nullValue)
                                doubleValues[n] = value0;
                            else
                                doubleValues[n] = computeValue(doubleValues[n], value0);
                        }
                    }
                }
                for(int n=0; n<nPoints; n++)
                {
                    if((wizard.mathVolume.getOperator() == StsMathVirtualVolume.AVERAGE) && (doubleValues[n] != StsParameters.nullValue))
                        doubleValues[n] = doubleValues[n]/(double)vols.length;
                    if(doubleValues[n] < Float.MAX_VALUE)
                        floatValues[n] = (float)doubleValues[n];
                }
                panel.setProgress((float)((float)i/(float)vols[0].nSlices) * 33.0f);
                sFile.write(computeUnscaledByteValues(floatValues, ttlMin, ttlMax));
            }
            for(h=0; h<vols.length; h++)
            {
                vols[h].clearCache();
            }
            floatValues = null;
            doubleValues = null;

            panel.setText("Processing Crossline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nRows;
            floatValues = new float[nPoints];
            doubleValues = new double[nPoints];

            for(i = 0; i <= vols[0].nCols; i++)
            {
                for(j=0; j<nPoints; j++)
                {
                    doubleValues[j] = StsParameters.nullValue;
                    floatValues[j] = StsParameters.nullValue;
                }
                for(h = 0; h<vols.length; h++)
                {
                    min = wizard.mathVolume.getDataMin();

                    // Disabled Volumes are automatically excluded but still create a frame
                    if(!vols[h].getIsVisible())
                        continue;

                    plane = vols[h].readBytePlaneData(StsCursor3d.XDIR, vols[h].getXCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(plane[n] != -1)
                        {
                            double value0 = vols[h].getScaledValue(plane[n]);
                            if((value0 > max) || (value0 < min))
                                continue;
                            if(doubleValues[n] == StsParameters.nullValue)
                                doubleValues[n] = value0;
                            else
                                doubleValues[n] = computeValue(doubleValues[n], value0);
                        }
                    }
                }
                for(int n=0; n<nPoints; n++)
                {
                    if((wizard.mathVolume.getOperator() == StsMathVirtualVolume.AVERAGE) && (doubleValues[n] != StsParameters.nullValue))
                        doubleValues[n] = doubleValues[n]/(double)vols.length;
                    if(doubleValues[n] < Float.MAX_VALUE)
                        floatValues[n] = (float)doubleValues[n];
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nCols) * 33.0f) + 33.0f);
                cFile.write(computeUnscaledByteValues(floatValues, ttlMin, ttlMax));
            }
            for(h=0; h<vols.length; h++)
                vols[h].clearCache();
            floatValues = null;
            doubleValues = null;

            panel.setText("Processing Inline Data");
            Main.logUsage();

            nPoints = vols[0].nSlices * vols[0].nCols;
            floatValues = new float[nPoints];
            doubleValues = new double[nPoints];

            for(i = 0; i <= vols[0].nRows; i++)
            {
                for(j=0; j<nPoints; j++)
                {
                    doubleValues[j] = StsParameters.nullValue;
                    floatValues[j] = StsParameters.nullValue;
                }
                for(h = 0; h<vols.length; h++)
                {
                    min = wizard.mathVolume.getDataMin();
                    if(!vols[0].getIsVisible())
                        continue;

                    plane = vols[h].readBytePlaneData(StsCursor3d.YDIR, vols[h].getYCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(plane[n] != -1)
                        {
                            double value0 = vols[h].getScaledValue(plane[n]);
                            if((value0 > max) || (value0 < min))
                                continue;
                            if(floatValues[n] == StsParameters.nullValue)
                                doubleValues[n] = value0;
                            else
                                doubleValues[n] = computeValue(doubleValues[n], value0);
                        }
                    }
                }
                for(int n=0; n<nPoints; n++)
                {
                    if((wizard.mathVolume.getOperator() == StsMathVirtualVolume.AVERAGE) && (doubleValues[n] != StsParameters.nullValue))
                        doubleValues[n] = doubleValues[n]/(double)vols.length;
                    if(doubleValues[n] < Float.MAX_VALUE)
                        floatValues[n] = (float)doubleValues[n];
                }
                panel.setProgress(((float)((float)i/(float)vols[0].nRows) * 33.0f) + 67.0f);
                rFile.write(computeUnscaledByteValues(floatValues, ttlMin, ttlMax));
            }
            for(h=0; h<vols.length; h++)
            {
                vols[h].clearCache();
            }
            floatValues = null;
            doubleValues = null;

            calculateHistogram();

            try
            {
                StsSeismicVolume combineVolume = new StsSeismicVolume(false);
                StsToolkit.copySubToSuperclass(vols[0], combineVolume,
                                               StsSeismicBoundingBox.class, StsBoundingBox.class, false);
                combineVolume.setName(name);
                combineVolume.setStsDirectory(model.getProject().getRootDirString());
                setFilenames(combineVolume, name);
                combineVolume.setDataMin(ttlMin);
                combineVolume.setDataMax(ttlMax);
                combineVolume.addToModel();
                combineVolume.addToProject(true);
                combineVolume.dataHist = dataHist;
                combineVolume.setDataHistogram();
                combineVolume.initialize(model);
                combineVolume.writeHeaderFile();
                closeOutputFiles();
            }
            catch(Exception e)
            {
                StsException.outputException("StsSetupCombinationVolumecreateSummationVolume() failed.",
                                             e, StsException.WARNING);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("createMathVolume failed.", e, StsException.WARNING);
        }
    }

	private void setFilenames(StsSeismicVolume volume, String name)
	{
		volume.createVolumeFilenames(name);
	}

    private void createCummulativeSubVolumeVolume()
    {

    }
    private void createDifferenceVolume()
    {
        StsSeismicVolume vol1 = wizard.diffVolume.getVolumeOne();
        StsSeismicVolume vol2 = wizard.diffVolume.getVolumeTwo();
        min = wizard.diffVolume.getDataMin();
        max = wizard.diffVolume.getDataMax();
        String name = wizard.getVolumeName();

        float i;
        int j, k;
        byte[][] planes = new byte[2][];
        byte[] plane;

        try
        {
            getOutputFiles(name);
            clearHistogram();

            panel.setText("Processing Slice Data");
            Main.logUsage();

            int nPoints = vol1.nCols * vol1.nRows;
            float[] floatValues = new float[nPoints];

            for(i = 0; i <= vol1.nSlices; i++)
            {
                planes[0] = vol1.readBytePlaneData(StsCursor3d.ZDIR, vol1.getZCoor(i));
                planes[1] = vol2.readBytePlaneData(StsCursor3d.ZDIR, vol1.getZCoor(i));
                for (int n = 0; n < nPoints; n++)
                {
                    if(!isAValueNull(planes[0][n], planes[1][n]))
                    {
                        float value0 = vol1.getScaledValue(planes[0][n]);
                        float value1 = vol2.getScaledValue(planes[1][n]);
                        floatValues[n] = value0 - value1;
                        if((floatValues[n] > max) || (floatValues[n] < min))
                            floatValues[n] = StsParameters.nullValue;
                    }
                    else
                        floatValues[n] = StsParameters.nullValue;
               }
               sFile.write(computeUnscaledByteValues(floatValues, min, max));
               panel.setProgress((float)((float)i/(float)vol1.nSlices) * 33.0f);
            }

            panel.setText("Processing Crossline Data");
            Main.logUsage();

            nPoints = vol1.nSlices * vol1.nRows;
            floatValues = new float[nPoints];
            for(i = 0; i <= vol1.nCols; i++)
            {

                planes[0] = vol1.readBytePlaneData(StsCursor3d.XDIR, vol1.getXCoor(i));
                planes[1] = vol2.readBytePlaneData(StsCursor3d.XDIR, vol1.getXCoor(i));
                for (int n = 0; n < nPoints; n++)
                {
                    if(!isAValueNull(planes[0][n], planes[1][n]))
                    {
                        float value0 = vol1.getScaledValue(planes[0][n]);
                        float value1 = vol2.getScaledValue(planes[1][n]);
                        floatValues[n] = value0 - value1;
                        if((floatValues[n] > max) || (floatValues[n] < min))
                            floatValues[n] = StsParameters.nullValue;
                    }
                    else
                        floatValues[n] = StsParameters.nullValue;
               }
               cFile.write(computeUnscaledByteValues(floatValues, min, max));
               panel.setProgress(((float)((float)i/(float)vol1.nCols) * 33.0f) + 33.0f);
            }

            panel.setText("Processing Inline Data");
            Main.logUsage();

            nPoints = vol1.nSlices * vol1.nCols;
            floatValues = new float[nPoints];
            for(i = 0; i <= vol1.nRows; i++)
            {
                planes[0] = vol1.readBytePlaneData(StsCursor3d.YDIR, vol1.getYCoor(i));
                planes[1] = vol2.readBytePlaneData(StsCursor3d.YDIR, vol1.getYCoor(i));
                for (int n = 0; n < nPoints; n++)
                {
                    if(!isAValueNull(planes[0][n], planes[1][n]))
                    {
                        float value0 = vol1.getScaledValue(planes[0][n]);
                        float value1 = vol2.getScaledValue(planes[1][n]);
                        floatValues[n] = value0 - value1;
                        if((floatValues[n] > max) || (floatValues[n] < min))
                            floatValues[n] = StsParameters.nullValue;
                    }
                    else
                        floatValues[n] = StsParameters.nullValue;
               }
               rFile.write(computeUnscaledByteValues(floatValues, min, max));
               panel.setProgress(((float)((float)i/(float)vol1.nRows) * 33.0f) + 67.0f);
            }
            calculateHistogram();

            StsSegyVolume newSegY = buildSegYObject(name, vol1);
            StsParameterFile.writeObjectFields(model.getProject().getRootDirString() + StsSeismicBoundingBox.group3d + "." + name, newSegY, StsSegyVolume.class, StsBoundingBox.class);

            closeOutputFiles();

            // Load the seismic volume
            StsSeismicVolume seismicVolume = StsSeismicVolume.checkLoadFromStemname(model, model.getProject().getRootDirString(), name, true);
            if(seismicVolume == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING,"Failed to add volume to project: " + name);
                return;
            }

        }
        catch(Exception e)
        {
            StsException.outputException("createDifferenceVolume failed.", e, StsException.WARNING);
        }
    }

    private boolean isAValueNull(byte value0, byte value1)
    {
        return value0 == -1 || value1 == -1;
    }

    private void getOutputFiles(String name)
    {
        String filename;
        try
        {
            String directory = model.getProject().getRootDirString();
            filename = StsSeismicBoundingBox.createRowByteFilename(StsSeismicBoundingBox.group3d, name);
            rFile = new RandomAccessFile(directory + filename, "rw");
            filename = StsSeismicBoundingBox.createColByteFilename(StsSeismicBoundingBox.group3d, name);
            cFile = new RandomAccessFile(directory + filename, "rw");
            filename = StsSeismicBoundingBox.createSliceByteFilename(StsSeismicBoundingBox.group3d, name);
            sFile = new RandomAccessFile(directory + filename, "rw");
        }
        catch(Exception e)
        {
            StsException.outputException("StsSetupMicroseismicVolume: Creation of output files failed", e, StsException.WARNING);
        }
    }

    private StsSegyVolume buildSegYObject(String stemname, StsSeismicVolume vol)
    {
        svol.segyDirectory = vol.getSegyDirectory();
        svol.stsDirectory = vol.stsDirectory;
        svol.segyFilename = vol.segyFilename;
        svol.stemname = stemname;
        svol.setFilenames(stemname);
        svol.setZDomain(StsParameters.TD_TIME);
        if(vol.getIsDepth())
            svol.setZDomain(StsParameters.TD_DEPTH);
        svol.segyLastModified = vol.segyLastModified;
        svol.isXLineCCW = vol.isXLineCCW;
        svol.isRegular = vol.isRegular;

        svol.dataMin = min;
        svol.dataMax = max;

        svol.nRows = vol.nRows; svol.nCols = vol.nCols; svol.nSlices = vol.nSlices;
        svol.xInc = vol.xInc;  svol.yInc = vol.yInc; svol.zInc = vol.zInc;
        svol.rowNumMin = vol.rowNumMin; svol.rowNumMax = vol.rowNumMax; svol.rowNumInc = vol.rowNumInc;
        svol.colNumMin = vol.colNumMin; svol.colNumMax = vol.colNumMax; svol.colNumInc = vol.colNumInc;
        svol.angle = vol.angle;
        svol.xMin = vol.xMin; svol.yMin = vol.yMin; svol.zMin = vol.zMin;
        svol.xMax = vol.xMax; svol.yMax = vol.yMax; svol.zMax = vol.zMax;
        svol.xOrigin = vol.xOrigin; svol.yOrigin = vol.yOrigin;
        svol.originSet = vol.originSet;
        return svol;
    }

    private void closeOutputFiles()
    {
        try
        {
            sFile.close();
            rFile.close();
            cFile.close();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSetupMicroseismicVolume: Closing of output files failed", e, StsException.WARNING);
        }
    }

    private void accumulateHistogram(byte bindex)
    {
        int index = StsMath.signedByteToUnsignedInt(bindex);
        if(index > 254) index = 254;
        if(index < 0) index = 0;
        dataCnt[index] = dataCnt[index] + 1;
        ttlHistogramSamples++;
    }
    private void calculateHistogram()
    {
        for(int i=0; i<255; i++)
            dataHist[i] = (float)((float)dataCnt[i]/(float)ttlHistogramSamples)*100.0f;
        svol.dataHist = dataHist;
    }
    private void clearHistogram()
    {
        for(int i=0; i< 255; i++)
        {
            dataCnt[i] = 0;
            dataHist[i] = 0.0f;
        }
        ttlHistogramSamples = 0;
    }

    private float[] recalcScale(float[] values, float[] range)
    {
        for(int i=0; i<values.length; i++)
        {
            if((values[i] < range[0]) && (values[i] != StsParameters.nullValue))
                range[0] = values[i];
            if((values[i] > range[1]) && (values[i] != StsParameters.nullValue))
                range[1] = values[i];
        }
        return range;
    }

    private byte[] computeUnscaledByteValues(float[] values, float min, float max)
    {
        if(values == null) return null;
        int nValues = values.length;

        float scale = 254 / (max - min);

        byte[] planeData = new byte[nValues];
        for (int n = 0; n < nValues; n++)
        {
            if(values[n] == StsParameters.nullValue)
                planeData[n] = -1;
            else
            {
                float scaledFloat = (values[n] - min) * scale;
                int scaledInt = Math.round(scaledFloat);
                planeData[n] = StsMath.unsignedIntToUnsignedByte(scaledInt);
                accumulateHistogram(planeData[n]);
            }
        }
        return planeData;
    }

}
