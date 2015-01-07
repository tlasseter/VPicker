//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsWellLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsWellWizard wizard;
    private String[] selectedWellnames = null;
    private StsPoint[] topPoints = null;
    int[] types = null;
    StsAbstractWellFactory wellFactory = null;
    private StsWell[] wells = null;

    private boolean canceled = false;
    public StsWellLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsWellWizard)wizard;        
        panel = StsProgressPanel.constructor(5, 50);
        header = new StsHeaderPanel();
        super.initialize(wizard, panel, null, header);    
        header.setTitle("Well Selection");
        header.setSubtitle("Load Well(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
        header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen.");
    }
/*
    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Well Selection");
        header.setSubtitle("Load Well(s)");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/WellLoad.html");
        header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen.");
    }
*/
    public void setWellFactory(StsAbstractWellFactory wellFactory)
    {
        this.wellFactory = wellFactory;
    }

    public void setSelectedWellnames(String[] selectedWellnames)
    {
        this.selectedWellnames = selectedWellnames;
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

            // turn on the wait cursor
            StsCursor cursor = new StsCursor(panel, Cursor.WAIT_CURSOR);
//            if(wizard.getFileType() == wizard.SINGLE_WELL)
            wells = StsWellImport.createWells(panel, selectedWellnames, wizard.selectedWells, wizard.wellheadXYZs, wizard.applyKbs, wellFactory, wizard.getDatumShift(), wizard.hUnits,  wizard.vUnits, model.getProject().getDateOrder());

            success = (wells != null);
            int nSucessfulWells = 0;
            for(int i=0; i<wells.length; i++)
            {
            	if(wells[i] != null)
            		nSucessfulWells++;
            }
            if(wizard.getArchiveIt())
            {
                panel.setDescription("Archiving well data.");
                String archiveDir = model.getProject().getDataFullDirString() + model.getProject().getArchiveDirString();
                if(!StsWellImport.archiveWells(selectedWellnames, archiveDir))
                    new StsMessage(model.win3d, StsMessage.ERROR, "Source file(s) does not exist or cannot be copied to database.");
                panel.appendLine("Archiving well data complete.");
            }
            completeWellLoad();

            panel.appendLine("Well loading is complete. Press the Finish> button");

            if(nSucessfulWells == 0)
                 panel.setDescriptionAndLevel("All wells failed to load.", StsProgressBar.ERROR);
             else if(nSucessfulWells < wells.length)
                 panel.setDescriptionAndLevel("Some wells failed to load.", StsProgressBar.WARNING);
             else
                 panel.setDescriptionAndLevel("All wells loaded successfully.", StsProgressBar.INFO);
            panel.finished();

            cursor.restoreCursor();
            if(nSucessfulWells > 0) model.getProject().runCompleteLoading();
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