package com.Sts.Actions.Build;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

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
public class StsWellOnCursor extends StsLineOnCursor
{
    public StsWellOnCursor(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public void initializeLine()
    {
        line = new StsWell();
    }
}
