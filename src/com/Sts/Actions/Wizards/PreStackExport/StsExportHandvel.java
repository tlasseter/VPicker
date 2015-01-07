package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

public abstract class StsExportHandvel extends StsWizardStep {
	StsPreStackHandVelExportPanel panel;
    protected StsHeaderPanel header;
    protected StsPreStackVelocityModel velocityModel;

    public StsExportHandvel(StsWizard wizard)
	{
        super(wizard);
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet == null) {
        	System.err.println("line Set is null");
        	return;
        }
        velocityModel = lineSet.velocityModel;
        if (lineSet == null) {
        	System.err.println("velocity Model is null");
        	return;
        }
        panel = createHandVelPanel();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        //panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Export Handvel Velocities");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/PreStackDataExport.html");
        header.setInfoText(wizardDialog,"(1) Select the type of pre-stack data you wish to export.\n" +
                           " (2) Select the dimension of data.\n" +
                           "     *** Disabled if project only contains 2D or 3D but not both.\n" +
                           "(3) Press the Next > Button to proceed to the data selection step.\n");
    }

    protected abstract StsPreStackHandVelExportPanel createHandVelPanel();

	public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public void exportHandVels() {
    	StsProgressPanel progressPanel = StsProgressPanel.constructor(0, velocityModel.getNProfiles());
    	panel.exportHandVels(progressPanel);
    }
}

