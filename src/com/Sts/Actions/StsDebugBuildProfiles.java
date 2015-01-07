package com.Sts.Actions;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsDebugBuildProfiles extends StsAction
{
    String title;
    public StsDebugBuildProfiles(StsActionManager actionManager)
    {
        super(actionManager);
        title = "Debug build profiles";
    }

 	public boolean start()
    {
        statusArea.setTitle(title);
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if(lineSet == null) return true;
        StsPreStackVelocityModel velocityModel = lineSet.velocityModel;
        if(velocityModel == null) return true;
        velocityModel.debugBuildProfiles();
        actionManager.endCurrentAction();
        return true;
    }
}
