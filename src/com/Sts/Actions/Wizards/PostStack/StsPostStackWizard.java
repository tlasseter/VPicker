package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.MVC.*;

abstract public class StsPostStackWizard extends StsSeismicWizard implements Runnable
{
    protected StsPostStackTraceDefinition traceDef;
    public StsPostStackRangeEdit rangeEdit;
    protected StsPostStackBatch filesToProcess;
    protected StsPostStackProcess volumeProcess;

    public StsPostStackWizard(StsActionManager actionManager, int width, int height, byte segyType)
    {
        super(actionManager, width, height, segyType);
    }
}
