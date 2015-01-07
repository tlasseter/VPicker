
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Pointset;

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

public class StsSetupPointSet extends StsWizardStep implements Runnable
{
    public StsStatusPanel panel;
    private StsHeaderPanel header;
    private StsToolbar subVolumeToolbar = null;
    private StsCreatePointSetWizard wizard = null;

    RandomAccessFile rFile, sFile, cFile;
    float min, max;
    transient private int dataCnt[] = new int[255];
    transient private int ttlHistogramSamples = 0;
    public float[] dataHist = new float[255]; // Histogram of the data distribution
    StsSegyVolume svol = new StsSegyVolume();
    boolean autoScale = true;

    public StsSetupPointSet(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        panel = (StsStatusPanel) getContainer();
        this.wizard = (StsCreatePointSetWizard)wizard;

        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("PointSet Definition");
        header.setSubtitle("Produce PointSet");
        header.setInfoText(wizardDialog,"(1) Press the Cancel Button to stop the processing.\n" +
                                   "(2) Press the Finish Button once processing is complete.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/PointSet.html");
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
            if(createPointSet())
            {
                panel.setText("Setup Complete");
                panel.setProgress(100.0f);
                success = true;
            }
            else
            {
                success = false;
            }
        }
        catch(Exception e)
        {
            success = false;
        }
        finally
        {
            if(success)
            {
                disableCancel();
                wizard.enableFinish();
                wizard.finish();
            }
        }
    }

    public boolean end()
    {
        return true;
    }

    private boolean createPointSet()
    {
        int j, k, h, i;
        byte[] plane;
        int nPoints;
        StsPointList pointSet = null;
        boolean append = false;

        StsSeismicVolume[] vols = wizard.definePointSet.getVolumes();
        int[] count = new int[vols.length];
        for(h = 0; h<vols.length; h++)
            count[h] = 0;

        min = wizard.definePointSet.getDataMin();
        max = wizard.definePointSet.getDataMax();
        String psName = wizard.definePointSet.panel.getPointSetName();
        if((model.getObjectWithName(StsPointList.class, psName) != null) && (wizard.definePointSet.getPointSetToAppend() == null))
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Invalid Poinset Name: Already in use");
            wizard.gotoStep(wizard.definePointSet);
            return false;
        }

        if(wizard.definePointSet.getPointSetToAppend() != null)
        {
            append = true;
            pointSet = wizard.definePointSet.getPointSetToAppend();
        }

        try
        {
            clearHistogram();

            panel.setText("Processing Data");
            Main.logUsage();

            nPoints = vols[0].nCols * vols[0].nRows;

            // Determine Min and Max for Pointset
            float psTtlMin = StsParameters.largeFloat;
            float psTtlMax = wizard.definePointSet.getDataMax();

            // Determine the total minimum and maximum of the volume
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for (h = 0; h < vols.length; h++) {
                    min = wizard.definePointSet.getDataMin();
                    if (!vols[h].getIsVisible())
                        continue;

                    // Calculate min and max as percentage if required.
                    if (wizard.definePointSet.getSetType() ==
                        StsPointList.PERCENTAGE) {
                        double calcMin = calcPercentageMin(vols[h], max);
                        if (calcMin > min)
                            min = (float)calcMin;
                    }

                    // Calculate Pointset Min and Max
                    if (min < psTtlMin)
                        psTtlMin = (float) min;

//                    plane = (byte[])vols[h].readPlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
//                    for (int n = 0; n < nPoints; n++)
//                    {
//                        if(plane[n] != -1)
//                        {
//                            double value0 = vols[h].getScaledValue(plane[n]);
//                            if((value0 > max) || (value0 < min))
//                                continue;
//
//                            if (value0 > psTtlMax)
//                                psTtlMax = (float)value0;
//                            if (value0 < psTtlMin)
//                                psTtlMin = (float)value0;
//                        }
//                    }
                }
            }

//            System.out.println("Pointset range (" + psTtlMin + "-" + psTtlMax + ")");
            if(psTtlMax > StsParameters.doubleNullValue)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Maximum value is too large, constrain parameters");
                return false;
            }

            // Initialize the pointset
            int volIdx = 0;
            if(append)
            {
                if(psTtlMin < pointSet.getDataMin())
                    pointSet.setDataMin(psTtlMin);
                if(psTtlMax > pointSet.getDataMax())
                    pointSet.setDataMax(psTtlMax);
                volIdx = pointSet.getNumberVolumes();
                pointSet.addVolumeNames(wizard.definePointSet.panel.getVolumeNames());
            }
            else
            {
                pointSet = new StsPointList(100, psTtlMin, psTtlMax, false);
                pointSet.setName(wizard.definePointSet.panel.getPointSetName());
                pointSet.setVolumeNames(wizard.definePointSet.panel.getVolumeNames());
                if (wizard.definePointSet.getSetType() == StsPointList.PERCENTAGE)
                    pointSet.setPointSetType(pointSet.PERCENTAGE);
                else
                    pointSet.setPointSetType(pointSet.VALUE);
            }

            // Output the volumes and pointset
            for (i = 0; i <= vols[0].nSlices; i++)
            {
                for(h = 0; h<vols.length; h++)
                {
                    min = psTtlMin;
//                    System.out.println("Before Percentage PostStack3d[" + h + "]=" + vols[h] + " with range (" + min + "-" + max + ")");

                    // Disabled Volumes are automatically excluded but still create a frame
                    if(!vols[h].getIsVisible())
                        continue;

                    // Calculate min and max as percentage if required.
                    if(wizard.definePointSet.getSetType() == StsPointList.PERCENTAGE)
                    {
                        double calcMin = calcPercentageMin(vols[h], max);
                        if(calcMin > min)
                            min = (float)calcMin;
                    }
//                    System.out.println("PostStack3d[" + h + "]=" + vols[h] + " with range (" + min + "-" + max + ")");
                    plane = vols[h].readBytePlaneData(StsCursor3d.ZDIR, vols[h].getZCoor(i));
                    for (int n = 0; n < nPoints; n++)
                    {
                        if(plane[n] != -1)
                        {
                            double value0 = vols[h].getScaledValue(plane[n]);
//                            if(value0 > min - 1)
//                                System.out.println("Value=" + value0);
                            if((value0 > max) || (value0 < min))
                                continue;

                            float z = vols[h].getZCoor(i);
                            float x = vols[h].getXCoor((n - ((n/vols[h].nCols) * vols[0].nCols)));
                            float y = vols[h].getYCoor(n/vols[h].nCols);
                            StsPoint point = new StsPoint(5);
                            point.setValues(new float[] {x,y,z,(float)volIdx + h, (float)value0});
                            pointSet.addPoint(point);
//                            System.out.println("Added point with value = " + value0);
                            count[h]++;
                        }
                    }
                }
                panel.setProgress((float)((float)i/(float)vols[0].nSlices));
            }
            for(h=0; h<vols.length; h++)
            {
                vols[h].clearCache();
                if(count[h] == 0)
                {
                    float z = vols[h].getZCoor(0);
                    float x = vols[h].getXCoor(0);
                    float y = vols[h].getYCoor(0);
                    StsPoint point = new StsPoint(5);
                    point.setValues(new float[] {x,y,z,(float)volIdx + h, psTtlMin});
                    pointSet.addPoint(point);
                }
            }
            if(pointSet != null)
            {
                pointSet.trimPointsArray();
                pointSet.buildHistogram();
            }
            calculateHistogram();

            try
            {
                pointSet.addToModel();
            }
            catch(Exception e)
            {
                StsException.outputException("StsSetupCombinationVolumecreateSummationVolume() failed.",
                                             e, StsException.WARNING);
                return false;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("createMathVolume failed.", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    private double calcPercentageMin(StsSeismicVolume vol, float max)
    {
        float count = 0.0f;
        int index = vol.dataHist.length - 1;
        int numSamples = vol.nSlices * vol.nRows * vol.nCols;
        float percentage = wizard.definePointSet.getSetPercentage();
        double interval = (double)(((double)vol.getDataMax() - (double)vol.getDataMin())/(double)vol.dataHist.length);
        if(max != vol.getDataMax())
            index = (int)(((double)((double)max - (double)vol.getDataMin())/interval) - (double)1.0);
        if(index > vol.dataHist.length - 1)
            index = vol.dataHist.length - 1;
        for(int i=index; i>0; i--)
        {
            count += vol.dataHist[i];
            if(count > percentage)
            {
                return (double)max - (interval * (index-i+1));
            }
        }
        return vol.getDataMin();
    }

    private boolean isAValueNull(byte value0, byte value1)
    {
        return value0 == -1 || value1 == -1;
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
