package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>e
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineKbAndKickoffPanel extends StsFieldBeanPanel
{
    private StsDefineWellPlanWizard wizard;
    private StsDefineKbAndKickoff definePath;

    StsWellPlan planSet = null;

//    StsGroupBox wellOrPlatformBox = new StsGroupBox("Select/define platform or single well");
//    StsComboBoxFieldBean wellOrPlatformBean = new StsComboBoxFieldBean();

    StsGroupBox topHolePropertiesBox = new StsGroupBox("Well Plan Basic Description");
    StsStringFieldBean xyUnitsBean = new StsStringFieldBean("Horizontal units:");
    StsStringFieldBean zUnitsBean = new StsStringFieldBean("Vertical units:");
    StsStringFieldBean nameBean = new StsStringFieldBean();
    StsColorComboBoxFieldBean colorBean = new StsColorComboBoxFieldBean();
    StsDoubleFieldBean xOriginBean = new StsDoubleFieldBean();
    StsDoubleFieldBean yOriginBean = new StsDoubleFieldBean();
    StsFloatFieldBean zKbBean = new StsFloatFieldBean();
    StsFloatFieldBean zKickoffBean = new StsFloatFieldBean();
    StsFloatFieldBean buildRateBean = new StsFloatFieldBean();
    StsComboBoxFieldBean trajectoryTypeBean = new StsComboBoxFieldBean();
    StsFloatFieldBean intervalBean = new StsFloatFieldBean();

    public StsDefineKbAndKickoffPanel(StsWizard wizard)
    {
        this.wizard = (StsDefineWellPlanWizard)wizard;
        try
        {
            constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize(StsWellPlan planSet)
    {
        this.planSet = planSet;
        setPanelObject(planSet); // connect beans to this object
        StsProject project = wizard.getModel().getProject();
        xyUnitsBean.setValue(project.getXyUnitString());
        zUnitsBean.setValue(project.getDepthUnitString());
        xOriginBean.setValue(planSet.getXOrigin());
        yOriginBean.setValue(planSet.getYOrigin());
        zKbBean.setValue(planSet.getZKB());
        colorBean.doSetValueObject(planSet.getStsColor());
        if(planSet.getPlatform() == null)
        {
            xOriginBean.setEditable(true);
            yOriginBean.setEditable(true);
            zKbBean.setEditable(true);
        }
    }

    void constructPanel() throws Exception
    {
        nameBean.initialize(StsWellPlan.class, "name", "wellPlan", true, "Name:");
        StsSpectrum spectrum = StsObject.getCurrentModel().getSpectrum("Basic");
        colorBean.initializeColors(StsWellPlan.class, "stsColor", "Color:", spectrum);
//        xyUnitsBean.classInitialize("Horizontal units:");
//        zUnitsBean.classInitialize("Vertical units:");
        xOriginBean.classInitialize(StsWellPlan.class, "xOrigin", false, "X Origin:");
        yOriginBean.classInitialize(StsWellPlan.class, "yOrigin", false, "Y Origin:");
        zKbBean.classInitialize(StsWellPlan.class, "zKB", false, "Kelly Bushing sea level elevation:");
        zKickoffBean.classInitialize(StsWellPlan.class, "zKickoff", true, "Kickoff sealevel depth:");
        buildRateBean.classInitialize(StsWellPlan.class, "buildRate", true, "Build rate (deg/100):");
        trajectoryTypeBean.classInitialize(StsWellPlan.class, "trajectoryType", "Trajectory:", StsWellPlan.trajectoryTypeStrings);
        intervalBean.classInitialize(StsWellPlan.class, "displayInterval", 1.0f, 1000.0f, "Display interval:");

        topHolePropertiesBox.add(nameBean);
        topHolePropertiesBox.add(colorBean);
        topHolePropertiesBox.add(xyUnitsBean);
        topHolePropertiesBox.add(zUnitsBean);
        topHolePropertiesBox.add(xOriginBean);
        topHolePropertiesBox.add(yOriginBean);
        topHolePropertiesBox.add(zKbBean);
        topHolePropertiesBox.add(zKickoffBean);
        topHolePropertiesBox.add(buildRateBean);
        topHolePropertiesBox.add(trajectoryTypeBean);
        topHolePropertiesBox.add(intervalBean);
        add(topHolePropertiesBox);
    }
}
