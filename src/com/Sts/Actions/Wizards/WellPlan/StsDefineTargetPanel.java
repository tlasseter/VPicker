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
public class StsDefineTargetPanel extends StsFieldBeanPanel
{
    private StsDefineWellPlanWizard wizard;
    private StsDefineTarget defineTarget;
    private StsWellPlan wellPlan;

    private StsPoint[] points = null;

    StsStringFieldBean xyUnitsBean = new StsStringFieldBean(false, "Horizontal Units:");
    StsStringFieldBean zUnitsBean = new StsStringFieldBean(false, "Vertical Units:");

    StsGroupBox targetBox = new StsGroupBox("Target description");
	StsDoubleFieldBean targetXBean = new StsDoubleFieldBean(StsWellPlan.class, "targetX", "Target X:");
	StsDoubleFieldBean targetYBean = new StsDoubleFieldBean(StsWellPlan.class, "targetY", "Target Y:");
	StsFloatFieldBean targetTBean = new StsFloatFieldBean(StsWellPlan.class, "targetT", "Target Time:");
    StsFloatFieldBean targetZBean = new StsFloatFieldBean(false, "Target Depth:");

    StsGroupBox dropBox = new StsGroupBox("Drop arc description");
	StsFloatFieldBean dropRadiusBean = new StsFloatFieldBean(StsWellPlan.class, "dropRate", "Drop rate (deg/100):");
	StsFloatFieldBean dzAboveTargetBean = new StsFloatFieldBean(StsWellPlan.class, "dzAboveTarget", "Drop dz above target:");
	StsFloatFieldBean dropAngleBean = new StsFloatFieldBean(StsWellPlan.class, "dropAngle", "Drop angle:");

    StsProject project = null;

    public StsDefineTargetPanel(StsWizard wizard, StsDefineTarget defineTarget)
    {
        this.wizard = (StsDefineWellPlanWizard) wizard;
        this.defineTarget = defineTarget;
		this.wellPlan = this.wizard.getWellPlan();
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

    public void initialize(StsWellPlan wellPlan)
    {
        this.wellPlan = wellPlan;
//        planSet.initializeLegs();
        targetBox.setPanelObject(wellPlan); // connect beans to this object
		dropBox.setPanelObject(wellPlan);
        xyUnitsBean.setValue(project.getXyUnitString());
        zUnitsBean.setValue(project.getDepthUnitString());

        remove(dropBox);
        if(wellPlan.getTrajectoryType() == StsWellPlan.BUILD_HOLD_DROP_STRING)
        {
            add(dropBox);
        }
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        targetBox.add(xyUnitsBean);
        targetBox.add(zUnitsBean);
        targetBox.add(targetXBean);
        targetBox.add(targetYBean);
        targetBox.add(targetTBean);
        targetBox.add(targetZBean);
        add(targetBox);

        dropBox.add(dropRadiusBean);
        dropBox.add(dzAboveTargetBean);
        dropBox.add(dropAngleBean);
        add(dropBox);
    }

    /** Called when point is manually picked; set in display box and in object */
    public void setTargetX(double x)
    {
        targetXBean.setValue(x);
        wellPlan.setTargetX(x);
    }

    /** Called when point is manually picked; set in display box and in object */
    public void setTargetY(double y)
    {
        targetYBean.setValue(y);
        wellPlan.setTargetY(y);
    }

    /** Called when point is manually picked; set in display box and in object */
   public void setTargetT(float t)
    {
        targetTBean.setValue(t);
        wellPlan.setTargetT(t);
    }

    public void setTargetZ(float z)
    {
        targetZBean.setValue(z);
    }
}
