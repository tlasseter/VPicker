package com.Sts.Actions.Wizards.EditWellMarkers;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.StsMessage;

public class StsEditWellMarkersWizard extends StsWizard
{
    private StsEditWellMarkers editWellMarkers = new StsEditWellMarkers(this);
    private StsWizardStep[] wizardSteps = { editWellMarkers };

    public StsEditWellMarkersWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(wizardSteps);
   }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public void disableCancelBtn()
    {
        disableCancel();
    }

    public boolean start()
    {

        dialog.setTitle("Edit Well Markers");
        if(model.getProject().isRealtime())
        {
            new StsMessage(frame, StsMessage.WARNING, "Cannot define a marker while running real-time. Time will be stopped.");
            model.getProject().stopProjectTime();
        }
        return super.start();
    }

    public boolean checkPrerequisites()
    {
        if((model.getCreateStsClass(StsWell.class).getSize() == 0) && (model.getCreateStsClass(StsLiveWell.class).getSize() == 0))
        {
        	reasonForFailure = "No wells were found, unable to run edit markers wizard.";
    		return false;
        }
    	return true;
    }

    public boolean end()
    {
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        return super.end();
    }

}
