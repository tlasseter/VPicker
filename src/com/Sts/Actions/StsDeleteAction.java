package com.Sts.Actions;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsDeleteAction extends StsAction
{
    private StsObject[] objects = null;

 	public StsDeleteAction(StsActionManager actionManager, StsObject[] objects)
    {
        super(actionManager, true);
        this.objects = objects;
    }

    public boolean start()
    {
        if(objects == null) return false;
        int nObjects = objects.length;

        // Ask user if they are sure
        String message = "Delete:";
        int nLoops = nObjects;
        if(nObjects > 10) nLoops = 10;
        for(int n = 0; n < nLoops; n++)
        {
            message = message + "\n   ";
            message = message + objects[n].getName();
        }
        if(nLoops == 10)
            message = message + "\n   ....\n   " + objects[nObjects-1] + "?\n";
        else
            message = message + "?\n";
        boolean deleteThem = StsYesNoDialog.questionValue(getModel().win3d, message);
        if(!deleteThem)
            return false;

        boolean viewChanged = false;
        for(int n = 0; n < nObjects; n++)
        {
            objects[n].delete();
            logMessage("Successfully deleted: " + objects[n].getName());
            viewChanged = viewChanged | model.viewObjectChanged(this, objects[n]);
        }
        actionManager.endCurrentAction();
        if(viewChanged) model.win3dDisplayAll();
        model.refreshObjectPanel();
        return true;
   	}

    public boolean end()
    {
        return true;
    }
}
