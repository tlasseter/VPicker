package com.Sts.Actions.Wizards.SubVolumeMgmt;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSubVolumeMgmtApplyPanel extends StsJPanel
{
    private StsSubVolumeMgmtWizard wizard;
    private StsSubVolumeMgmtApply wizardStep;

    StsGroupBox objectBox = new StsGroupBox("Data Objects");
    StsBooleanFieldBean seismicApplyBean = null;
    StsBooleanFieldBean virtualVolumeApplyBean = null;
    StsBooleanFieldBean crossplotApplyBean = null;

    public StsSubVolumeMgmtApplyPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsSubVolumeMgmtWizard)wizard;
        this.wizardStep = (StsSubVolumeMgmtApply)wizardStep;
        try
        {
            constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void constructBeans()
    {
        StsModel model = wizard.getModel();
        seismicApplyBean = new StsBooleanFieldBean(model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume"),
            "displayOnSubVolumes", "Seismic Volumes");
        virtualVolumeApplyBean = new StsBooleanFieldBean(model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"),
            "displayOnSubVolumes", "Virtual Volumes");
        crossplotApplyBean = new StsBooleanFieldBean(model.getCreateStsClass("com.Sts.DBTypes.StsCrossplot"),
            "displayOnSubVolumes", "Crossplots");
    }

    public void initialize()
    {
        ;
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        objectBox.addEndRow(seismicApplyBean);
        objectBox.addEndRow(virtualVolumeApplyBean);
        objectBox.addEndRow(crossplotApplyBean);
        add(objectBox);
    }
}
