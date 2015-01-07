package com.Sts.Actions;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.MVC.*;

// Here we just want to toggle the 3d cursor
public class StsDisplaySpectrumDialog extends StsAction
{
    boolean success = false;

    public StsDisplaySpectrumDialog(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public void displayDialog() {
        StsSpectrumAction sd = new StsSpectrumAction((StsActionManager) actionManager);
        return;
    }

    public boolean start()
    {
        success = true;
        return success;
    }

    public boolean end()
    {
        success = true;
        return true;
    }
}
