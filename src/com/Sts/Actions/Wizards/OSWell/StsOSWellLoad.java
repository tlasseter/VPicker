//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsOSWellLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsOSWellWizard wizard;
    private StsOSWell[] osWells;
    private StsWell[] wells;
    private boolean reloadOriginals = false;
    static String currentDirectory = null;
    static String binaryDataDir;
    static float logNull;
    public int nLoaded = 0;
    static boolean success = false;
    StsOSWellIO wellIO = null;
    private boolean canceled = false;
    public StsOSWellLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsOSWellWizard)wizard;
        panel = StsProgressPanel.constructor(5, 50);
        header = new StsHeaderPanel();
        super.initialize(wizard, panel, null, header);
        header.setTitle("Well Selection");
        header.setSubtitle("Load Well(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#OSWell");
        header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen.");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void run()
    {
        try
        {
            if (canceled)
            {
                success = false;
                return;
            }
            disablePrevious();

            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting well loading...");

            osWells = wizard.getSelectedWells();

            // turn on the wait cursor
            StsCursor cursor = new StsCursor(panel, Cursor.WAIT_CURSOR);

            // get the well data from OpenSpirit, we add the data to the existing osWells
            osWells = wizard.currentDatastore.getOpenSpiritImport().createStsOSWells(osWells, panel);
            panel.appendLine("Completed loading selected OpenSpirit wells. Converting to S2S Wells...");

            //wells = StsWellImport.createFromOSWells(panel, wizard.getWellDatastore(), osWells);
            wells = createFromOSWells(panel, osWells);

            success = (wells != null);
            if(!success)
                return;

            int numFailed = 0;
            for(int i=0; i<wells.length; i++)
            {
            	if(wells[i] == null)
            		numFailed++;
            }
            completeWellLoad();
            panel.appendLine("Well loading is complete. Press the Finish> button");

            if(numFailed == wells.length)
                panel.setDescriptionAndLevel("All wells failed to load.", StsProgressBar.ERROR);
            else if(numFailed > 0)
                panel.setDescriptionAndLevel("Some wells failed to load.", StsProgressBar.WARNING);
            else
                panel.setDescriptionAndLevel("All wells loaded successfully.", StsProgressBar.INFO);
            panel.finished();

            cursor.restoreCursor();
            model.win3d.cursor3d.initialize();
            model.win3d.cursor3dPanel.setSliderValues();
            model.win3d.getGlPanel3d().setDefaultView();
            model.enableDisplay();
            model.win3dDisplay(); // display the wells
            wizard.enableFinish();
        }
        catch (Exception e)
        {
        	panel.appendLine("Failed to load well.");
            panel.appendLine("Error message: " + e.getMessage());
            StsException.outputWarningException(this, "run", e);
            success = false;
        }
    }

    public StsWell[] createFromOSWells(StsProgressPanel progressPanel, StsOSWell[] osWells)
    {
        try
        {
            if (osWells == null) return null;
            int nWells = osWells.length;

            progressPanel.appendLine("Preparing to load " + nWells + " well from OpenSpirit ...");
            StsMessageFiles.logMessage("Preparing to load " + nWells + " well from OpenSpirit ...");

            StsProject project = model.getProject();
            currentDirectory = project.getRootDirString();
            binaryDataDir = project.getBinaryFullDirString();
            logNull = project.getLogNull();

            // If the user has selected both binary files and ASCII files the units may be mixed, so reload all from ASCII
            int nWithBinaries = 0;
            for (int n = 0; n < nWells; n++)
            {
                // If binaries exist for some and not for others. Load all from ASCII
                if (StsWellKeywordIO.binaryFileExist(binaryDataDir, osWells[n].getName(), StsLogVector.WELL_DEV_PREFIX, StsLogVector.X))
                {
                    nWithBinaries++;
                }
            }
            if (nWithBinaries != nWells)
            {
                progressPanel.appendLine("All binary wells will be re-loaded from ASCII to apply unit selections\n");
                StsMessageFiles.infoMessage("All binary wells will be re-loaded from ASCII to apply unit selections");
                reloadOriginals = true;
            }
            else if (nWithBinaries == nWells)
            {
                progressPanel.appendLine("Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option");
                StsMessageFiles.infoMessage("Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option"); ;
            }

            StsWell[] wells = createWells(progressPanel, osWells);

            project.adjustBoundingBoxes(true, false); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsWell.class);
            project.rangeChanged();

            progressPanel.appendLine("Loaded " + nLoaded + " of " + nWells + " well deviation files ...");

            StsSeismicVelocityModel velocityModel = project.getSeismicVelocityModel();
            if (velocityModel != null)
            {
                for (int n = 0; n < nLoaded; n++)
                {
					wells[n].adjustFromVelocityModel(velocityModel);
                    wells[n].computeMarkerTimesFromMDepth(velocityModel);
                }
            }
            //        model.win3d.resetCursorPanel();
            return wells;

        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.createWells() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsWell[] createWells(StsProgressPanel progressPanel, StsOSWell[] osWells)
    {
        StsWell well;

        try
        {
            int nWells = osWells.length;
            int nLoaded = 0;
            StsWell[] wells = new StsWell[nWells];
            progressPanel.initialize(nWells);
            for (int n = 0; n < nWells; n++)
            {
                String wellname = osWells[n].getName();

                // read and build a well (if we don't already have it)
                if ( (model.getObjectWithName(StsWell.class, wellname) == null) || reloadOriginals)
                {
                    StsWell wellTemp = (StsWell) model.getObjectWithName(StsWell.class, wellname);
                    if ((wellTemp != null) && reloadOriginals)
                    {
                        progressPanel.appendLine("Deleting well " + wellname + " from database to re-load.");
                        wellTemp.delete();
                    }
                    well = createWellFromOSWell(osWells[n], progressPanel);
                    if (well != null)
                    {
                        progressPanel.appendLine("Successfully created well " + wellname + " from database.");
                        wells[nLoaded++] = well;
                    }
                }
                else
                {
                    if (progressPanel != null)
                    {
                        progressPanel.appendLine("Well: " + wellname + " already loaded...");
                        progressPanel.setDescriptionAndLevel("Well: " + wellname + " already loaded...", StsProgressBar.WARNING);
                    }
                }
                progressPanel.setValue(n+1);
                progressPanel.setDescription("Loaded well #" + nLoaded + " of " + nWells);
            }
            return wells;

        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "createWells", e);
            progressPanel.setDescriptionAndLevel("StsWellWizard.createWells() failed.\n", StsProgressBar.WARNING);
            return null;
        }
    }

    private StsWell createWellFromOSWell(StsOSWell osWell, StsProgressPanel progressPanel)
    {
        StsOSWellIO wellIO = new StsOSWellIO(model, osWell, false);
        StsWell well = wellIO.createWell(logNull, 0.0f, progressPanel);
        success = (well != null);
        if(!success)
            return null;

        well.addToProject(); // this adds lineVertices to projectBoundingBoxes

        model.getProject().adjustBoundingBoxes(true, false); // extend displayBoundingBox as needed and set cursor3d box accordingly
        model.getProject().checkAddUnrotatedClass(StsWell.class);
        model.getProject().rangeChanged();

        well.computePoints(); // generate splined points between vertices

        wellIO.addWellMarkers(well);
        return well;
    }

    public StsWell[] getWells()
    {
        return wells;
    }

    public void completeWellLoad()
    {
        disableCancel();
        enableFinish();
    }

    public boolean end()
    {
        return success;
    }
}
