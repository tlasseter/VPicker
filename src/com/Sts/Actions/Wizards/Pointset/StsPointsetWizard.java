
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Pointset;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsPointsetWizard extends StsWizard
{
    public StsPointsetSelect selectPointset = null;
    public StsPointsetLoad loadPointset = null;

    private StsWizardStep[] mySteps =
    {
        selectPointset = new StsPointsetSelect(this),
        loadPointset = new StsPointsetLoad(this)
    };

    public StsPointsetWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Pointset Import");
        dialog.getContentPane().setSize(600, 600);
        return super.start();
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    // Persist any palettes that are not already in database
    public boolean createPointsets()
    {
        float progress, incProgress;

        StsAbstractFile[] files = selectPointset.getSelectedFiles();
        String[] names = selectPointset.getSelectedPointsetNames();
        String dirname = selectPointset.getSelectedDirectory();

        try
        {
            loadPointset.panel.setMaximum(names.length * 10); // 10 is the scaling: max = n*scaling
            progress = 0.0f;
            incProgress = 10.0f;
            disablePrevious();
            loadPointset.panel.setProgress(progress);
            Thread.currentThread().sleep(10);

            StsPointListClass pointsetClass = (StsPointListClass)model.getCreateStsClass(StsPointList.class);
            StsPointList pointset = null;
            StsPoint[] points = null;
            for (int i = 0; i < names.length; i++)
            {
                if (pointsetClass.getPointSet(names[i]) == null)
                {
                    // Add the palette
                    points = readPointset(files[i]);
                    pointset = new StsPointList(names[i], points, null, true);
                    if(pointset != null)
                        loadPointset.panel.setText("Loaded pointset: " + names[i]);
                    else
                        loadPointset.panel.setText("Failed to load pointset: " + names[i]);
                }
                else
                {
                    loadPointset.panel.setText("Pointset: " + names[i] + " already in database.");
                }
                progress += incProgress;
                loadPointset.panel.setProgress(progress);
                disableCancel();
                if(pointset != null)
                    pointset.addToModel();
            }
            enableFinish();
        }
        catch (Exception e)
        {
            StsException.outputException("StsPointsetWizard.createPointsets() failed.", e, StsException.WARNING);
        }
        return true;
    }

    public void finish()
    {
        super.finish();
    }

    public StsPoint[] readPointset(StsAbstractFile file)
    {
        StsAsciiFile asciiFile = null;
        float x, y, z, v, a;
        StsPoint point = null;
        StsPoint[] points = null;
        int incSize = 100;

        asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return null;

        String filename = file.getFilename();
        String[] tokens;
        boolean firstLine = true;
        String line = null;
        int nPoints = 0;
        try
        {
            while((line = asciiFile.readLine()) != null)
            {
                StringTokenizer stok = new StringTokenizer(line, ",");
                if(firstLine)
                {
                    firstLine = false;
                    continue;
                }
                int nTokens = stok.countTokens();
                if(nTokens < 3)
                {
                    String inputLine = asciiFile.getLine();
                    StsMessageFiles.errorMessage("Insufficient entries for line: " + inputLine);
                }
                else if(nTokens == 3)
                {
                    x = Float.parseFloat((String)stok.nextElement());
                    y = Float.parseFloat((String)stok.nextElement());
                    z = Float.parseFloat((String)stok.nextElement());
                    point = new StsPoint(x,y,z);
                }
                else if(nTokens == 4)
                {
                    x = Float.parseFloat((String)stok.nextElement());
                    y = Float.parseFloat((String)stok.nextElement());
                    z = Float.parseFloat((String)stok.nextElement());
                    v = Float.parseFloat((String)stok.nextElement());
                    point = new StsPoint(x,y,z,v);
                }
                else
                {
                    x = Float.parseFloat((String)stok.nextElement());
                    y = Float.parseFloat((String)stok.nextElement());
                    z = Float.parseFloat((String)stok.nextElement());
                    v = Float.parseFloat((String)stok.nextElement());
                    a = Float.parseFloat((String)stok.nextElement());
                    point = new StsPoint(x,y,z,v,a);
                }
                if(point == null)
                    continue;

                points = (StsPoint[])StsMath.arrayAddElement(points, point, nPoints, incSize);
                nPoints++;

                if(points == null)
                {
                    StsException.systemError("Cannot create pointset for: " + filename +
                                             " due to StsMath.arrayAddElement error.");
                    return null;
                }
            }
            points = (StsPoint[])StsMath.trimArray(points, nPoints);
        }
        catch (java.io.IOException e)
        {
            ;
        }
        asciiFile.close();
        return points;
    }

}
