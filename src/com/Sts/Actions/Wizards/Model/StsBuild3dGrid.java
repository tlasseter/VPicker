
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;

public class StsBuild3dGrid extends StsWizardStep implements Runnable
{
    // convenience flag copies
    static final int NOT_BUILDABLE = StsZone.NOT_BUILDABLE;
    static final int BUILDABLE = StsZone.BUILDABLE;
    static final int BUILT = StsZone.BUILT;

    private StsStatusPanel status;
    private StsHeaderPanel header;

    public StsBuild3dGrid(StsWizard wizard)
    {
    	super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        status = (StsStatusPanel) getContainer();

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Model Construction");
        header.setSubtitle("Build 3D Model");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Model");        
        header.setInfoText(wizardDialog,"(1) Press the Next>> Button when model construction is finished");
    }

	public boolean start()
    {
        status.setTitle("Building the 3d grid:");
        run();
		return true;
    }

    public void run()
    {
        StsBuiltModel builtModel = ((StsModelWizard)wizard).getBuiltModel();
        success = builtModel.buildModel(model, status);
		if(success)
        {
        	status.setText("3-D model built successfully.");
        	logMessage("3-D model built successfully.");
        }
        else
		{
        	status.setText("Unable to build 3-D model.");
        	logMessage("Unable to build 3-D model.");
		}

 //       status.sleep(1000);
//		wizard.finish();
    }


    public boolean end()
    {
		return success;
    }
}
