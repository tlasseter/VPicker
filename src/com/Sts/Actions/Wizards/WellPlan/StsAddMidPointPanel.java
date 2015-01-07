package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>e
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsAddMidPointPanel extends StsFieldBeanPanel
{
    private StsDefineWellPlanWizard wizard;
    private StsAddMidPoint addMidPoint;
//    private StsWellPlan plannedWell;
    private StsWellPlan wellPlan;

    private StsPoint[] points = null;
    StsProject project = null;

    StsStringFieldBean xyUnitsBean = new StsStringFieldBean("Horizontal units:");
    StsStringFieldBean zUnitsBean = new StsStringFieldBean("Vertical units:");

    StsGroupBox midPointBox = new StsGroupBox("Mid point description");
    StsDoubleFieldBean midPointXBean = new StsDoubleFieldBean();
    StsDoubleFieldBean midPointYBean = new StsDoubleFieldBean();
    StsFloatFieldBean midPointTBean = new StsFloatFieldBean();
    StsFloatFieldBean midPointZBean = new StsFloatFieldBean(false, "midPoint Depth:");
    StsFloatFieldBean midPointRateBean = new StsFloatFieldBean();

    public StsAddMidPointPanel(StsWizard wizard, StsAddMidPoint addMidPoint)
    {
        this.wizard = (StsDefineWellPlanWizard) wizard;
        this.addMidPoint = addMidPoint;
        project = wizard.getModel().getProject();
        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize(StsWellPlan planSet)
    {
        this.wellPlan = planSet;
        midPointBox.setPanelObject(planSet); // connect beans to this object
        xyUnitsBean.setValue(project.getXyUnitString());
        zUnitsBean.setValue(project.getDepthUnitString());
    }

    void jbInit() throws Exception
    {
//        xyUnitsBean.classInitialize("Horizontal units:");
//        zUnitsBean.classInitialize("Vertical units:");
        midPointXBean.classInitialize(StsWellPlan.class, "midPointX", true, "MidPoint X:");
        midPointYBean.classInitialize(StsWellPlan.class, "midPointY", true, "MidPoint Y:");
        midPointTBean.classInitialize(StsWellPlan.class, "midPointT", true, "MidPoint Time:");
//        midPointZBean.classInitialize("midPoint Depth:");
        midPointRateBean.classInitialize(StsWellPlan.class, "midPointRate", true, "MidPoint rate (deg/100):");

        midPointBox.add(xyUnitsBean);
        midPointBox.add(zUnitsBean);
        midPointBox.add(midPointXBean);
        midPointBox.add(midPointYBean);
        midPointBox.add(midPointTBean);
        midPointBox.add(midPointZBean);
        midPointBox.add(midPointRateBean);
        add(midPointBox);
    }

    /** Called when point is manually picked; set in display box and in object */
    public void setMidPointX(double x)
    {
        midPointXBean.setValue(x);
        wellPlan.setMidPointX(x);
    }

    /** Called when point is manually picked; set in display box and in object */
    public void setMidPoint(double y)
    {
        midPointYBean.setValue(y);
        wellPlan.setMidPointY(y);
    }

    /** Called when point is manually picked; set in display box and in object */
   public void setMidPointT(float t)
    {
        midPointTBean.setValue(t);
        wellPlan.setMidPointT(t);
    }

    public void setMidPoint(float z)
    {
        midPointZBean.setValue(z);
    }
}
