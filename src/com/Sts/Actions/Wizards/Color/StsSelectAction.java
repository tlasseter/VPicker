package com.Sts.Actions.Wizards.Color;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSelectAction extends StsAction implements Runnable
{
    boolean success = false;
    StsColorscale colorscale = null;
    StsSpectrumSelect spectrumSelect = null;

    public StsSelectAction(StsActionManager actionManager, StsColorscalePanel colorscalePanel)
    {
        super(actionManager);
        spectrumSelect = new StsSpectrumSelect(colorscalePanel);
    }

    public void run()
    {
        try
        {
            spectrumSelect.setVisible(true);
            success = spectrumSelect.getSuccess();
            if(success)
                actionManager.endCurrentAction();
            return;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSelectAction.run() failed.", e, StsException.WARNING);
            return;
        }
    }


    public boolean end()
    {
        return success;
    }


}
