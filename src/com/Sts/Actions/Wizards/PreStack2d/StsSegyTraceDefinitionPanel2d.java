package com.Sts.Actions.Wizards.PreStack2d;

import com.Sts.Actions.Wizards.PreStack3d.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSegyTraceDefinitionPanel2d extends StsSegyTraceDefinitionPanel
{
    public StsSegyTraceDefinitionPanel2d(StsWizard wizard, StsSegyPreStackTraceDefinition wizardStep)
    {
        super( wizard, wizardStep);
    }

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "xOrigin", "yOrigin", "colNumMin", "colNumMax", "statusString"};
        String[] columnTitles = {"Name", "Start X", "Start Y", "Start Cdp", "End Cdp", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Volumes");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }
}
