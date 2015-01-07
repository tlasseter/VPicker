package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

import javax.swing.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsPostStackTraceDefinitionPanel2d extends StsPostStackTraceDefinitionPanel implements ChangeListener, ListSelectionListener,
                                                                                    StsSelectRowNotifyListener, StsTableModelListener
{
	public StsPostStackTraceDefinitionPanel2d(StsSeismicWizard wizard, StsWizardStep wizardStep)
	{
        super(wizard, wizardStep);
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
