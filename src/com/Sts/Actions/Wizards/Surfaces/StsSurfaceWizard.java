//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsSurfaceWizard extends StsWizard
{
    int nSurfaces;
    String jarFilename = "grids.jar";

    byte vUnits = StsParameters.DIST_NONE;
    byte binaryVertUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte binaryHorzUnits = StsParameters.DIST_NONE;
    byte tUnits = StsParameters.TIME_NONE;
    byte binaryTimeUnits = StsParameters.TIME_NONE;
    float datumShift = 0.0f;

    int[] order = new int[5];

    boolean concatLines = true;

    public StsSurfaceSelect selectProcessSurfaces = null;
    public StsDefineRows defineRows = null;
    public StsDefineColumns defineCols = null;
    public StsSurfaceLoad processSurfaces = null;

    private StsWizardStep[] mySteps =
    {
        selectProcessSurfaces = new StsSurfaceSelect(this),
        defineRows = new StsDefineRows(this),
        defineCols = new StsDefineColumns(this),
        processSurfaces = new StsSurfaceLoad(this)
    };

    /** file types processed by this loader */
    static final String WEBJAR = StsAbstractFile.WEBJAR;
    static final String JAR = StsAbstractFile.WEBJAR;
    static final String BINARY = StsAbstractFile.binaryFormat;
    static final String SEISMIC = StsSurface.seismicGrp;
    static final String ZMAP = StsSurface.zmapGrp;


    public StsSurfaceWizard(StsActionManager actionManager)
    {
        super(actionManager, 600, 600);

        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Process & Load Surfaces");
        disableFinish();
        return super.start();
    }

    public boolean end()
    {
        if (success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        if (currentStep == defineRows)
        {
            int nFiles = selectProcessSurfaces.getSelectedFiles().length;
            gotoPreviousStep();
            if(nFiles > 0) enableNext(true);
        }
        else
            gotoPreviousStep();        
    }

    public void next()
    {
        if (currentStep == selectProcessSurfaces)
        {
            if (selectProcessSurfaces.getViewableFiles().length == 0)
                gotoStep(processSurfaces);
            else
                gotoNextStep();
        }
        else
            gotoNextStep();
    }

    public StsAbstractFile[] getSelectedFiles()
    {
        return selectProcessSurfaces.getSelectedFiles();
    }

    public StsAbstractFile[] getViewableFiles()
    {
        return selectProcessSurfaces.getViewableFiles();
    }

    public StsAbstractFile getCurrentFile()
    {
        return selectProcessSurfaces.panel.getCurrentFile();
    }

    public String[] getSelectedSurfacenames()
    {
        return selectProcessSurfaces.panel.getSelectedSurfaceNames();
    }

    public void finish()
    {
        super.finish();
    }

    public void setOrder(byte field, int col)
    {
        order[field] = col;
    }

    public void runCreateSurfaces()
    {
        Runnable createSurfacesRunnable = new Runnable()
        {
            public void run()
            {
                createSurfaces();
            }
        };

        Thread createSurfacesThread = new Thread(createSurfacesRunnable);
        createSurfacesThread.start();
    }

    /** Create surfaces */
    public boolean createSurfaces()
    {
        int nLoaded = 0;
        try
        {
            StsAbstractFile[] files = selectProcessSurfaces.getSelectedFiles();
            int nSelected = files.length;

            vUnits = selectProcessSurfaces.panel.unitsPanel.vUnits;
            hUnits = selectProcessSurfaces.panel.unitsPanel.hUnits;
            tUnits = selectProcessSurfaces.panel.unitsPanel.tUnits;
            concatLines = defineCols.panel.concatLines();

            datumShift = selectProcessSurfaces.panel.verticalUnitsPanel.datumShift;

            // turn off cursors and wait for redisplay
//             model.win3d.cursor3dPanel.setEnabled(false);
//            model.win3dDisplayAndWait();

            // don't draw 3d window until surfaces built
            model.disableDisplay();

            disablePrevious();

            // turn on the wait cursor
            StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);

            logMessage("Preparing to load " + selectProcessSurfaces.getSelectedFiles().length + " surfaces ...");

            StsProject project = model.getProject();
            StsProgressPanel panel = processSurfaces.panel;
            panel.progressBar.initializeIntervals(nSelected, 100);
            boolean reloadAscii = selectProcessSurfaces.getReloadAscii();
            String timeOrDepth = selectProcessSurfaces.panel.verticalUnitsPanel.timeDepthString;
            boolean isTrueDepth = timeOrDepth.equals(StsParameters.TD_DEPTH_STRING);
            String tvdOrSsString = selectProcessSurfaces.panel.verticalUnitsPanel.tvdOrSsString;
            boolean isTvd = tvdOrSsString.equals(StsParameters.TVD_STRING);
            for (int n = 0; n < nSelected; n++)
            {
                panel.progressBar.setInterval(n);
                StsAbstractFile file = files[n];
                String surfaceName = file.name;
                StsSurface surface = (StsSurface) model.getObjectWithName(StsSurface.class, surfaceName);
                if (surface != null && !reloadAscii)
                {
                    panel.appendLine("Surface: " + surfaceName + " already loaded...");
                    continue;
                }
                String fileType = file.fileType;
                if(fileType.equals(SEISMIC)) // could be a multiSurface file: process all surfaces
                {
                    StsImportSeismicSurfaces.initializeFile();
                    while((surface = createSurface(file, surfaceName, isTrueDepth, isTvd)) != null)
                    {
                        if (surface == null)
                            panel.appendLine("Failed to load " + surfaceName);
                        else
                            panel.appendLine("Loaded " + surfaceName + ": surface number " + ++nLoaded + " from file " + file.name + "...");
                    }
                }
                else
                {
                    surface = createSurface(file, surfaceName, isTrueDepth, isTvd);
                    if (surface == null)
                        panel.appendLine("Failed to load " + surfaceName);
                    else
                        panel.appendLine("Loaded " + surfaceName + ": surface number " + ++nLoaded + " from file " + file.name + "...");
                }
            }
            panel.finished();

            project.runCompleteLoading();
            // turn off the wait cursor
            if (cursor != null) cursor.restoreCursor();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsProcessSurfaceWizard.createSurfaces() failed.", e, StsException.WARNING);
            return false;
        }
        finally
        {
            success = nLoaded > 0;
            completeLoading(success);
        }
    }

    private StsSurface createSurface(StsAbstractFile file, String surfaceName, boolean isTrueDepth, boolean isTvd)
    {
        String fileType = file.fileType;
        if(fileType.equals(SEISMIC))
        {
            if (!StsImportSeismicSurfaces.scanSurface(model, file, surfaceName, datumShift, hUnits, vUnits, tUnits, concatLines, isTrueDepth))
                return null;
            return StsImportSeismicSurfaces.createSurface(model, file, surfaceName, concatLines, processSurfaces.panel);
        }
        else if(fileType.equals(ZMAP))
        {
            if (!StsImportAsciiSurfaces.scanSurface(model, file, surfaceName, datumShift, hUnits, vUnits, tUnits, concatLines, isTrueDepth, isTvd))
                return null;
            return StsImportAsciiSurfaces.createSurface(model, file, surfaceName, concatLines, processSurfaces.panel);
        }
        else if(fileType.equals(BINARY))
            return StsImportBinarySurfaces.createSurface(model, file, surfaceName, processSurfaces.panel);
        else if(fileType.equals(JAR))
            return StsImportBinarySurfaces.createSurface(model, file, jarFilename, processSurfaces.panel);
        else if(fileType.equals(WEBJAR))
            return StsImportBinarySurfaces.createSurface(model, file, surfaceName, processSurfaces.panel);
        else
            StsException.systemError(this, "createSurface", "Could not identify fileType: " + fileType);
            return null;
    }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsSurfaceWizard wellWizard = new StsSurfaceWizard(actionManager);
        wellWizard.start();
    }
}
