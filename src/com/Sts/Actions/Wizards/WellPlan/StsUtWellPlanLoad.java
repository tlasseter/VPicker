package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Well.*;
import com.Sts.DBTypes.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsUtWellPlanLoad extends StsWellLoad
{
    public StsUtWellPlanLoad(StsWizard wizard)
    {
        super(wizard);
    }

    public void completeWellLoad()
    {
        disableCancel();
        disableNext();
        enableFinish();
        StsWell[] wells = super.getWells();
        if (wells != null)
        {
            for (int n = 0; n < wells.length; n++)
            {
                StsWellPlanSet planSet = new StsWellPlanSet();
				planSet.setName(wells[n].getName());
                planSet.addWellPlan( (StsWellPlan) wells[n]);
            }
        }

        /*
         StsWellPlanSet planSet = ((StsWellPlanWizard)wizard).getWellPlanSet();
        StsWell[] wells = super.getWells();
        if(wells == null) return;
        planSet.addPlannedWell( (StsWellPlan) wells[0]);
        return;
    */
    }
}
