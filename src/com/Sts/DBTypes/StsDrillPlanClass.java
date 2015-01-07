package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsDrillPlanClass extends StsClass implements StsSerializable
{
	//static final Class[] subClassClasses = new Class[] { StsDrillPlan.class };
	//static StsDrillPlanClass[] subClasses = null;

	public StsDrillPlanClass()
	{
        userName = "Drill Plans";                                              
	}

    public void initializeFields()
	{
//        initializeSubClasses();
    }
/*
    private void initializeSubClasses()
    {
		int nSubClasses = subClassClasses.length;
		subClasses = new StsDrillPlanClass[nSubClasses];
		int nActualInstances = 0;
		for(int n = 0; n < nSubClasses; n++)
		{
			StsDrillPlanClass subClassInstance = (StsDrillPlanClass) currentModel.getStsClass(subClassClasses[n]);
			if(subClassInstance != null) subClasses[nActualInstances++] = subClassInstance;
		}
		subClasses = (StsDrillPlanClass[])StsMath.trimArray(subClasses, nActualInstances);
    }
*/
    public void initializeDisplayFields()
	{
//        initColors(StsDrillPlan.planDisplayFields);

		displayFields = new StsFieldBean[] { };
	}

	public void selected(StsDrillPlan drillPlan)
	{
		super.selected(drillPlan);
		setCurrentObject(drillPlan);
	}

	public StsDrillPlan getCurrentDrillPlan()
	{
		return (StsDrillPlan)currentObject;
	}

	public boolean setCurrentObject(StsObject object)
	{
		boolean changed = super.setCurrentObject(object);
		if (changed && object != null) ( (StsDrillPlan) object).treeObjectSelected();
		return changed;
	}

	public boolean setCurrentPlanName(String name)
	{
		StsDrillPlan newDrillPlan = (StsDrillPlan) getObjectWithName(name);
		return setCurrentObject(newDrillPlan);
	}

	public void close()
	{
		list.forEach("close");
	}

	static public StsDrillPlanClass getDrillPlanClass()
	{
		return (StsDrillPlanClass)currentModel.getCreateStsClass(StsDrillPlan.class);
	}

    public StsDrillPlan[] getDrillPlans()
    {
        Object drillPlanList;

        StsDrillPlan[] drillPlans = new StsDrillPlan[0];
        drillPlanList = currentModel.getCastObjectList(StsDrillPlan.class);
        drillPlans = (StsDrillPlan[])StsMath.arrayAddArray(drillPlans, drillPlanList);
        return drillPlans;
    }

}
