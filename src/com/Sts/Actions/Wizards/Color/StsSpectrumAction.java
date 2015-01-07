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
import com.Sts.Utilities.*;

public class StsSpectrumAction extends StsAction
{
    boolean success = false;
    StsColorscale colorscale = null;
    StsSpectrumDialog spectrumDialog = null;

    public StsSpectrumAction(StsActionManager actionManager)
    {
        super(actionManager);
        spectrumDialog = new StsSpectrumDialog(actionManager.getModel());
    }

    public StsSpectrumAction(StsActionManager actionManager, StsColorscale colorscale)
    {
        super(actionManager);
        this.colorscale = colorscale;
        spectrumDialog = new StsSpectrumDialog("Color Palette Editor", model, colorscale, true);
    }

    public boolean start()
    {
        try
        {
            spectrumDialog.setModal(true);
            spectrumDialog.setVisible(true);
            success = spectrumDialog.getSuccess();
            actionManager.endCurrentAction();
            return success;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrumDialog.setVisible(true) failed.", e, StsException.WARNING);
            return false;
        }
    }


    public boolean end()
    {
        if(spectrumDialog != null) spectrumDialog.dispose();
        return success;
    }


}
