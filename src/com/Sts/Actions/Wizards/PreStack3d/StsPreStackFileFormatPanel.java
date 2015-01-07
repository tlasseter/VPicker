package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.PreStack.*;
import com.Sts.Actions.Wizards.Seismic.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsPreStackFileFormatPanel extends StsPostStackFileFormatPanel
{
	public StsPreStackFileFormatPanel(StsSeismicWizard wizard, StsPostStackFileFormat wizardStep)
	{
        super(wizard, wizardStep);
	}

    protected void buildSubPanels()
    {
        unitsPanel = new StsPrestackUnitsPanel(wizard);
        fileFormatPanel = new StsFileFormatPanel(wizard, true);
        buildTablePanel();
    }
}
