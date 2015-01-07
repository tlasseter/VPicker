package com.Sts.Actions.Wizards.Well;

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

import com.Sts.Actions.Import.*;
import com.Sts.MVC.*;

public class StsWellDrillingWizard extends StsWellWizard
{
    public StsWellDrillingWizard(StsActionManager actionManager)
    {
        super(actionManager, false);
        StsWellImport.setReloadAscii(true);
    }
}
