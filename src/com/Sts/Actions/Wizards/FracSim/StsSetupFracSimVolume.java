
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSim;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.SubVolume.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;

public class StsSetupFracSimVolume extends StsWizardStep implements Runnable
{
    public StsStatusPanel panel;
    private StsHeaderPanel header;

    public StsSetupFracSimVolume(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        panel = (StsStatusPanel) getContainer();

        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sub-PostStack3d Definition");
        header.setSubtitle("Setup Sub-PostStack3d(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FracSim");
    }

    public boolean start()
    {
        panel.setTitle("Setting up Sub-PostStack3d Creation");
        run();
        return true;
    }

    public void run()
    {
        try
        {
            panel.setMaximum(100); // 10 is the scaling: max = n*scaling
            disablePrevious();
            panel.setProgress(0.0f);
            Thread.currentThread().sleep(10);
            panel.setText("Creating toolbar");
            checkAddToolbar();
            panel.setProgress(50.0f);

            Thread.currentThread().sleep(200);
            panel.setText("Setup Complete");
            panel.setProgress(100.0f);
            disableCancel();
            wizard.enableFinish();

            success = true;
        }
        catch(Exception e)
        {
            success = false;
        }
    }

    public void checkAddToolbar()
    {
        String toolbarName = null;
        StsSubVolumeWizard wizard = (StsSubVolumeWizard)this.wizard;
        StsSubVolume subVolume = null;
        toolbarName = null;
/*
        switch(wizard.getType())
        {
            case StsSubVolume.SINGLE_SURFACE:
                toolbarName = null;
//                subVolume = new StsSingleSurfaceSubVolume();
                break;
            case StsSubVolume.DUAL_SURFACE:
                toolbarName = null;
//                subVolume = new StsDualSurfaceSubVolume();
                break;
            case StsSubVolume.WELL_SET:
                toolbarName = null;
                break;
            case StsSubVolume.RESERVOIR_UNIT:
                toolbarName = null;
                break;
        }
        organizeToolbars(toolbarName);
*/
    }
/*
    private void organizeToolbars(String toolbarName)
    {
        if(model.win3d.hasToolbarNamed(StsPolyhedralTubeSubVolumeToolbar.NAME))
            model.win3d.closeToolbar(model.win3d.getToolbarNamed(StsPolyhedralTubeSubVolumeToolbar.NAME));
        if(model.win3d.hasToolbarNamed(StsPolyhedronSubVolumeToolbar.NAME))
            model.win3d.closeToolbar(model.win3d.getToolbarNamed(StsPolyhedronSubVolumeToolbar.NAME));

        if(toolbarName != null)
        {
            if (!model.win3d.hasToolbarNamed(toolbarName))
                model.win3d.addToolbar(subVolumeToolbar);
        }
    }
*/
    public boolean end()
    {
        return true;
    }

}