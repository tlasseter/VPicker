
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Super;

import com.Sts.Actions.Wizards.PostStack3d.*;
import com.Sts.Actions.Wizards.PostStack3dLoad.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsSuperWizard extends StsWizard implements Runnable
{
    final static int COLOR = 0;
    final static int SEGY = 1;
    final static int XPLOT = 2;
    final static int VVOL = 3;
    final static int GRID = 4;
    final static int WELL = 5;

    public StsSuperGeneral colorQA = null;
    private String colorQuestion = "Do you have any color palettes to load?";
    private String colorDescription = "Several default color palettes are available to the user and additional " +
        " palettes can be loaded and used. The allowed formats are: R G B, R G B A, or " +
        " I R G B A where I is the color index between 0-255, R is Red, " +
        " G is Green and B is Blue and A is Alpha or Transparency with all values being " +
        " between 0 and 255. The values are delimited by spaces or commas. Any loaded " +
        " palette will be converted to a 255 color palette, interpolating any missing indices.";

    public StsSuperGeneral segyQA = null;
    private String segyQuestion = "Do you have any SegY files to load?";
    private String segyDescription = "S2S requires all SegY Files to be pre-processed to optimize them " +
        "for visualization and analysis. The optimization process will re-organize " +
        "the data into three 8-bit cubes; one inline oriented, one crossline " +
        "oriented and one time slice oriented. A reference to the original " +
        "SegY File will alos be maintain to allow the software to access the " +
        "data in it\'s native resolution as required. The pre-processing step " +
        "will be followed by the loading of the pre-processed volumes. Upon " +
        "completion of this step, the SegY data should be isVisible in the " +
        "application.";

    public StsSuperGeneral xplotQA = null;
    private String xplotQuestion = "Do you want to create a crossplot?";
    private String xplotDescription = "S2S has a substantial crossplotting capability that allows the " +
        "user to plot any slice or sub-volume of data from one volume against the same slice or sub-volume from another " +
        "volume. Once the crossplot is defined it will provide the user with a crossplot toolbar " +
        "that is used to analyze the crossplot results.";

    public StsSuperGeneral vVolQA = null;
    private String vVolQuestion = "Do you want to create a virtual volume?";
    private String vVolDescription = "A virtual volume is a volume that is create from other actual volumes" +
        " or virtual volumes. Virtual volumes are created as the data is accessed, meaning that these volumes" +
        " are never stored on disk. There are four types of virtual volumes: volume math, crossplot, volume blend" +
        " and rgb volumes. PostStack3d math volumes are virtual volumes that are created by applying math or logical" +
        " operators to two volumes or one volume and a scalar. Crossplot virtual volumes simply null all points of" +
        " a selected volume based on crossplot results. Blended volumes involves the merging of two volumes based" +
        " on simple conditions. Rgb volumes are the belnding of three volumes where one is used to set the red value," +
        " one the green value, and one the blue value";

    public StsSuperGeneral gridQA = null;
    private String gridQuestion = "Do you have any gridded surfaces to load?";
    private String gridDescription = "S2S currently supports gridded surface files exported " +
                           "from most other software products. However, the format is " +
                           "fairly restrictive. <<Finish Description>>.";

    public StsSuperGeneral wellQA = null;
    private String wellQuestion = "Do you have any wells and logs to load?";
    private String wellDescription = "S2S currently supports LAS and an S2S developed well and log " +
                           "formats. Since the LAS file format does not support well deviation " +
                           "surveys S2S has introduced a LAS format for well bore path input. " +
                           "Therefore, to load logs from a LAS file an associated LAS like borehole " +
                           "file is also required. <<Finish Description>>.";

    private StsWizardStep[] mySteps =
    {
        colorQA = new StsSuperGeneral(this, COLOR, "Color Palette Load", colorQuestion, colorDescription),
        segyQA = new StsSuperGeneral(this, SEGY, "Segy File Load", segyQuestion, segyDescription),
        xplotQA = new StsSuperGeneral(this, XPLOT, "Crossplot Load or Creation", xplotQuestion, xplotDescription),
        vVolQA = new StsSuperGeneral(this, VVOL, "Virtual PostStack3d Creation", vVolQuestion, vVolDescription),
        gridQA = new StsSuperGeneral(this, GRID, "Gridded Surface Load", gridQuestion, gridDescription),
        wellQA = new StsSuperGeneral(this, WELL, "Well & Log Load", wellQuestion, wellDescription)
    };

    public StsSuperWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Super Wizard");
        dialog.getContentPane().setSize(300, 400);
        if(!super.start()) return false;
        return true;
    }

    private void runWorkflowWizard(Class actionClass, String workflowStepName)
    {
        StsWizard wizard = (StsWizard)actionManager.newWorkflowAction(actionClass, workflowStepName);
        wizard.enableInstructions();
        actionManager.checkStartAction(wizard);
    /*
        if(wizard instanceof Runnable)
        {
            Thread thread = new Thread((Runnable)wizard);
            thread.start();
        }
        else
            wizard.start();
    */
        try { actionManager.threadWait(); }
        catch(Exception e) { }
    }

    /** Launch the SegY process and load Wizards */
    public void hasSegY()
    {
        runWorkflowWizard(StsPostStack3dWizard.class, "Process Seismic");
        runWorkflowWizard(StsVolumeWizard.class, "Load Seismic");
    }

    /** Launch the crossplot Wizard */
    public void wantsXplot()
    {
        runWorkflowWizard(com.Sts.Actions.Wizards.CrossPlot.StsCrossplotWizard.class, "Create/Load Crossplot");
    }

    /** Launch the Virtual PostStack3d Wizard */
    public void wantsVirtualVolume()
    {
        runWorkflowWizard(com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard.class, "Create Virtual PostStack3d");
    }

    /** Launch the Palette Load Wizard */
    public void hasColor()
    {
        runWorkflowWizard(com.Sts.Actions.Wizards.Color.StsPaletteWizard.class, "Load Palettes");
    }

    /** Launch the Surfaces Load Wizard */
    public void hasGrids()
    {
        runWorkflowWizard(com.Sts.Actions.Wizards.Surfaces.StsSurfaceWizard.class, "Load Surfaces");
    }

    /** Launch the Wells Load Wizard */
    public void hasWells()
    {
        runWorkflowWizard(com.Sts.Actions.Wizards.Well.StsWellWizard.class, "Load Wells & Logs");
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
        Thread runnableThread;

        try
        {
            if((currentStep == colorQA) && (colorQA.getAnswer()))
            {
                Runnable colorRunnable = new Runnable()
                {
                    public void run()
                    {
                        hasColor();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(colorRunnable, "Color Thread");
                runnableThread.start();
            }
            else if((currentStep == segyQA) && (segyQA.getAnswer()))
            {
                Runnable segyRunnable = new Runnable()
                {
                    public void run()
                    {
                        hasSegY();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(segyRunnable, "Segy Thread");
                runnableThread.start();
            }
            else if((currentStep == xplotQA) && (xplotQA.getAnswer()))
            {
                Runnable xplotRunnable = new Runnable()
                {
                    public void run()
                    {
                        wantsXplot();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(xplotRunnable, "Xplot Thread");
                runnableThread.start();
            }
            else if((currentStep == vVolQA) && (vVolQA.getAnswer()))
            {
                Runnable vVolRunnable = new Runnable()
                {
                    public void run()
                    {
                        wantsVirtualVolume();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(vVolRunnable, "Vvol Thread");
                runnableThread.start();
            }
            else if((currentStep == gridQA) && (gridQA.getAnswer()))
            {
                Runnable surfaceRunnable = new Runnable()
                {
                    public void run()
                    {
                        hasGrids();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(surfaceRunnable, "Surface Thread");
                runnableThread.start();
            }
            else if((currentStep == wellQA) && (wellQA.getAnswer()))
            {
                Runnable wellRunnable = new Runnable()
                {
                    public void run()
                    {
                        hasWells();
                        gotoNextStep();
                    }
                };
                runnableThread = new Thread(wellRunnable, "Well Thread");
                runnableThread.start();
            }
            else
                gotoNextStep();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSuperWizard.next() failed.",
                e, StsException.WARNING);
        }
    }

    public void finish()
    {
        try
        {
            super.finish();
            actionManager.threadComplete();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSuperWizard.finish() failed.", e, StsException.WARNING);
        }
    }
}
