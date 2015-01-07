package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.DBTypes.*;

public class StsWellPlanFactory extends StsAbstractWellFactory
{
    static StsWellPlanFactory WellPlanFactory = null;
    static boolean isLive = false;
    
    public StsWellPlanFactory()
    {
    }

    static StsWellPlanFactory getWellPlanFactory()
    {
        if(WellPlanFactory == null) WellPlanFactory = new StsWellPlanFactory();
        return WellPlanFactory;
    }
    /**
     * createWellInstance
     *
     * @return StsWell
     * @todo Implement this com.Sts.DBTypes.StsAbstractWellFactory method
     */
    public StsWell createWellInstance()
    {
        return (StsWell)(new StsWellPlan(false));
    }

    public StsWell createWellInstance(String wellname)
    {
        return (StsWell)(new StsWellPlan(wellname, false));
    }
    public void setIsLive(boolean live)
    {
        isLive = live;
    }
}
