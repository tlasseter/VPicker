
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.io.*;

public class StsTimeSeriesWizard extends StsWizard implements Runnable
{
    public StsTimeSeriesAscii tsAscii = null;
    public StsTimeSeriesDefine tsDefine = null;
    public StsTimeSeriesPlot tsPlot = null;
    public StsTimeSeriesSelect tsSelect = null;

    private StsMovie selectedMovie = null;
    private StsFile asciiFile = null;

    public static final byte ACTUAL_TIME = 0;
    public static final byte ARTIFICIAL_TIME = 1;

    private StsWizardStep[] mySteps =
    {
        tsSelect = new StsTimeSeriesSelect(this),
        tsAscii = new StsTimeSeriesAscii(this),
        tsPlot = new StsTimeSeriesPlot(this),
    };

    public StsTimeSeriesWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        if(((StsMovie[])getModel().getCastObjectList(StsMovie.class)).length == 0)
        {
            new StsMessage(frame, StsMessage.ERROR,"No movies exist, unable to run time series wizard.");
            end();
            return false;
        }

        System.runFinalization();
        System.gc();
        dialog.setTitle("Time Series Animation");
        dialog.getContentPane().setSize(300, 300);

        if(!super.start()) return false;
        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
        if(currentStep == tsPlot)
            tsPlot.destroyChart();
        gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public File getAsciiFile()
    {
        return tsAscii.panel.getAsciiFile();
    }

    public int getTimeIndex()
    {
        return tsAscii.panel.getTimeIndex();
    }

    public int[] getAttributeIndices()
    {
        return tsAscii.panel.getAttributeIndices();
    }

    public int getNumberTokens()
    {
        return tsAscii.panel.getNumberTokens();
    }

    public int getNumberValidRows()
    {
        return tsAscii.panel.getNumberValidRows();
    }

    public long getAsciiStartTime()
    {
        return tsAscii.panel.getStartTime();
    }
    public byte getTimeType()
    {
        return tsAscii.panel.getAsciiTimeType();
    }
    public long getMovieStart()
    {
        return tsSelect.panel.getStartTime();
    }
    public int getMovieIncrement()
    {
        return tsSelect.panel.getIncrement();
    }
    public StsMovie getMovie()
    {
        return tsSelect.panel.getMovie();
    }
    public boolean isMovieStartClockTime()
    {
        return tsSelect.panel.isClockTime();
    }

    static String cleanTimeString(String ts)
    {
        if(ts.length() == 21)
            return ts;
        if(ts.indexOf(".") == -1)
            ts = ts + ".0";
        if(ts.indexOf(".") == 7)
            ts = "0" + ts;
        if((ts.indexOf(".") < 7) || (ts.indexOf(".") > 8))
            return null;
        ts = "2000-01-01 " + ts;
        return ts;
    }

}
