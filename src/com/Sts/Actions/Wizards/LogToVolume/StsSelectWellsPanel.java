package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectWellsPanel extends StsJPanel
{
    private StsLogToVolumeWizard wizard;
    private StsSelectWells wizardStep;

    private StsModel model = null;

    Object selectedWells = null;
    StsListFieldBean wellBean = new StsListFieldBean();
    StsGroupBox wellBox = new StsGroupBox("Select Wells");

    public StsSelectWellsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsLogToVolumeWizard)wizard;
    	this.wizardStep = (StsSelectWells)wizardStep;
    	model  = wizard.getModel();
        if(!buildPanel())
        	wizard.cancel();
    }

    public boolean buildPanel()
    {
        gbc.fill = gbc.HORIZONTAL;

        wellBox.gbc.fill = gbc.HORIZONTAL;
        wellBox.addEndRow(wellBean);
    	addEndRow(wellBox);

    	return true;
    }

    public void initialize()
    {
        if(model == null) return;
        StsWell[] wells =  (StsWell[])model.getCastObjectList(StsWell.class);
        StsWell[] wellsWithLogs = null;
        for(int i=0; i<wells.length; i++)
        {
            if(wells[i].getNLogCurves() == 0) continue;
            else
                wellsWithLogs = (StsWell[])StsMath.arrayAddElement(wellsWithLogs, wells[i]);
        }
        if(wellsWithLogs == null)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "No logs found in loaded wells. Exit wizard, load wells with logs and try again.");
            return;
        }
        wellBean.initialize(wizard, "selectedWells", "Wells:", wellsWithLogs);
        wizard.rebuild();
    }

    public Object[] getSelectedWells()
    {
    	return wellBean.getSelectedObjects();
    }
}