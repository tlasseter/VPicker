package com.Sts.Actions.Wizards.Alarms;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineValueAlarmPanel extends StsJPanel
{
    private StsAlarmsWizard wizard;
    private StsDefineValueAlarm wizardStep;

    private StsModel model = null;

    StsFloatFieldBean valueBean;
    StsFloatFieldBean thresholdBean;
    StsBooleanFieldBean insideBean;

    public StsDefineValueAlarmPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsAlarmsWizard)wizard;
    	this.wizardStep = (StsDefineValueAlarm)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        valueBean = new StsFloatFieldBean(wizard, "value", true, "Value:");
        valueBean.setToolTipText("Set value to test sensor attribute value against.");

        thresholdBean = new StsFloatFieldBean(wizard, "threshold", true, "Threshold:");
        thresholdBean.setToolTipText("Set percentage range around value.");

        insideBean = new StsBooleanFieldBean(wizard, "inside", true, "Inside:", true);
        insideBean.setToolTipText("Alarm sounds if inside range, or outside range.");

        addEndRow(valueBean);
        addEndRow(thresholdBean);
        addEndRow(insideBean);
    }

    public void initialize()
    {
    	;
    }

}